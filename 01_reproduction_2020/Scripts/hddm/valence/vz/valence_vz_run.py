"""Run and analyze only the valence_vz model; plotting/PPC are separate. Default mode is read-only check.

Run with WSL hddm2020 Python. Existing models/data are never overwritten.
"""
import argparse
import csv
from datetime import datetime
import json
from pathlib import Path
import pickle
import re
import sys

from valence_vz_fit import fit_one

PROJECT = Path(__file__).resolve().parents[4]


def inside(path, parent):
    return path == parent or parent in path.parents


def resolve_output(config):
    path = (PROJECT / config['output_root']).resolve()
    allowed = [(PROJECT / 'Output/hddm').resolve(),
               (PROJECT / 'Output/sensitivity').resolve()]
    if not any(inside(path, p) and path != p for p in allowed):
        raise ValueError('output_root must be a model subdirectory under Output/hddm')
    return path


def validate(config, existing=True):
    if config['task'] != 'valence' or config['model'] != 'vz':
        raise ValueError('Wrong configuration for valence_vz_run.py')
    if config['schema_version'] != 1 or config['model'] not in ['vtz', 'vz']:
        raise ValueError('Unsupported schema/model; new methods need their own implementation')
    if config['model_kwargs']['stim_col'] != 'stim' or config['model_kwargs']['split_param'] != 'v':
        raise ValueError('Unsupported coding')
    if config['sampling']['db'] != 'pickle' or config['sampling']['samples'] <= config['sampling']['burn']:
        raise ValueError('Invalid sampling settings')
    if config['sampling']['thin'] < 1:
        raise ValueError('thin must be positive')
    conditions = ['match', 'val', 'id'] if config['task'] == 'matching' else ['val', 'id']
    expected_depends = {'v': conditions}
    if config['model'] == 'vtz':
        expected_depends['t'] = conditions
    if config['model_kwargs']['depends_on'] != expected_depends or config['model_kwargs']['include'] != 'z' or config['model_kwargs']['drift_criterion'] is not False:
        raise ValueError('Model label/settings mismatch')
    runs = config['runs']
    if not runs or len({r['id'] for r in runs}) != len(runs) or len({r['seed'] for r in runs}) != len(runs):
        raise ValueError('Run IDs/seeds must be nonempty and unique')
    if any(not re.match(r'^r\d{2,}$', r['id']) for r in runs):
        raise ValueError('Invalid run ID')
    output = resolve_output(config)
    input_paths = {(PROJECT / r['input']['path']).resolve() for r in runs}
    if len(input_paths) != 1:
        raise ValueError('Runs must use the same input table')
    for run in runs:
        for role in (['input', 'model', 'database'] if existing else ['input']):
            path = (PROJECT / run[role]['path']).resolve()
            if not path.is_file():
                if role in ['model', 'database']:
                    raise FileNotFoundError('Full historical model/samples are local only. Set their paths in the config, refit, or use figure04_plot.py for bundled group plots: '+str(path))
                raise FileNotFoundError('Missing input table: '+str(path))
    with next(iter(input_paths)).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        if not set(config['data_checks']['required_columns']).issubset(reader.fieldnames):
            raise ValueError('Missing input columns')
        data = list(reader)
    if len(data) != config['data_checks']['rows'] or len({r['subj_idx'] for r in data}) != config['data_checks']['subjects']:
        raise ValueError('Input counts differ')
    return output


def split(x):
    import numpy as np
    n = x.shape[1] // 2
    return np.concatenate([x[:, :n], x[:, -n:]], axis=0)


def rank_normalize(x):
    from scipy.stats import rankdata, norm
    return norm.ppf((rankdata(x.ravel(), method='average') - .375) / (x.size + .25)).reshape(x.shape)


def rhat(x):
    import numpy as np
    w = np.var(x, axis=1, ddof=1).mean()
    if w <= 0:
        return float('inf')
    return float(np.sqrt(((x.shape[1] - 1) / x.shape[1] * w + np.var(x.mean(axis=1), ddof=1)) / w))


def ess(x):
    import numpy as np
    m, n = x.shape
    centered = x - x.mean(axis=1, keepdims=True)
    fft = np.fft.rfft(centered, n=2*n, axis=1)
    acov = np.fft.irfft(fft * fft.conjugate(), n=2*n, axis=1)[:, :n] / n
    w = np.var(x, axis=1, ddof=1).mean()
    var_plus = (n - 1) / n * w + np.var(x.mean(axis=1), ddof=1)
    if var_plus <= 0:
        return float('nan')
    rho = 1 - (w - acov.mean(axis=0)) / var_plus
    rho[0] = 1
    pairs = rho[:(n//2)*2:2] + rho[1:(n//2)*2:2]
    bad = np.flatnonzero(pairs <= 0)
    positive = np.minimum.accumulate(pairs[:int(bad[0]) if len(bad) else len(pairs)])
    tau = max(1., -1. + 2. * positive.sum())
    return float(min(m*n, m*n/tau))


def diagnostics(x):
    import numpy as np
    if x.shape[0] < 2:
        return dict(rank_split_rhat=None, bulk_ess=None, tail_ess=None)
    s = split(x)
    normalized = rank_normalize(s)
    return dict(rank_split_rhat=max(rhat(normalized), rhat(rank_normalize(abs(s - np.median(s))))),
                bulk_ess=ess(normalized), tail_ess=min(ess((s <= np.quantile(s, .05)).astype(float)),
                                                     ess((s <= np.quantile(s, .95)).astype(float))))


def load_traces(path, expected):
    import numpy as np
    # Only load trusted local historical/new PyMC databases.
    with path.open('rb') as stream:
        raw = pickle.load(stream)
    traces = {}
    for name, value in raw.items():
        if not isinstance(value, dict) or 0 not in value:
            continue
        try:
            array = np.asarray(value[0], dtype=float)
        except (ValueError, TypeError):
            continue
        if array.shape == (expected,):
            if not np.isfinite(array).all():
                raise ValueError('Nonfinite trace: ' + str(name))
            traces[str(name)] = array
    if 'z' not in traces or not any(k.startswith('v(') for k in traces):
        raise ValueError('Expected scalar traces absent or length mismatch')
    return traces


def table(path, rows):
    with path.open('x', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(config, output, source):
    import numpy as np
    expected = (config['sampling']['samples'] - config['sampling']['burn'] + config['sampling']['thin'] - 1) // config['sampling']['thin']
    chains, provenance = [], []
    for run in config['runs']:
        if source == 'existing':
            path = (PROJECT / run['database']['path']).resolve()
        else:
            folder = output / 'runs' / run['id']
            record = json.loads((folder / 'complete.json').read_text())
            saved = json.loads((folder / 'config.json').read_text())
            if record['seed'] != run['seed'] or saved['model_kwargs'] != config['model_kwargs'] or saved['sampling'] != config['sampling']:
                raise ValueError('New results and requested model settings differ')
            path = folder / 'posterior.db'
        chains.append(load_traces(path, expected))
        provenance.append(dict(run=run['id'],seed=run['seed'],database=str(path)))
    names = sorted(chains[0])
    if any(set(c) != set(names) for c in chains):
        raise ValueError('Chain parameter sets differ')
    out = output / ('analysis_' + source + '_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=False)
    (out / 'config.json').write_text(json.dumps(config, indent=2), encoding='utf-8')
    table(out / 'sources.csv', provenance)
    diag = []
    for name in names:
        row = dict(parameter=name)
        row.update(diagnostics(np.stack([c[name] for c in chains])))
        diag.append(row)
    table(out / 'diagnostics.csv', diag)
    thresholds = config['diagnostics']
    passed = len(chains) > 1 and all(
        np.isfinite(r['rank_split_rhat']) and r['rank_split_rhat'] < thresholds['rhat_max']
        and np.isfinite(r['bulk_ess']) and r['bulk_ess'] >= thresholds['bulk_ess_min']
        and np.isfinite(r['tail_ess']) and r['tail_ess'] >= thresholds['tail_ess_min'] for r in diag)
    summaries, contrasts = [], []
    scopes = [(r['id'], c) for r,c in zip(config['runs'], chains)]
    if passed:
        scopes.append(('pooled' + str(len(chains)), {n:np.concatenate([c[n] for c in chains]) for n in names}))
    group_names = [n for n in names if n in ['a','z','t','v'] or re.match(r'^[vt]\(',n)]
    for scope, values in scopes:
        for name in group_names:
            x = values[name]
            summaries.append(dict(scope=scope,parameter=name,mean=float(x.mean()),sd=float(x.std(ddof=1)),
                                  low=float(np.quantile(x,.025)),high=float(np.quantile(x,.975)),draws=len(x)))
        for left in [n for n in group_names if n.startswith('v(Self.') and n.endswith('Good)')]:
            for right in [left.replace('Self.','Other.',1),left.replace('Good)','Bad)')]:
                if right in values:
                    delta = values[left]-values[right]
                    contrasts.append(dict(scope=scope,left=left,right=right,mean=float(delta.mean()),
                        low=float(np.quantile(delta,.025)),high=float(np.quantile(delta,.975)),probability_gt_zero=float((delta>0).mean())))
        if 'z' in values:
            delta=values['z']-.5
            contrasts.append(dict(scope=scope,left='z',right='0.5',mean=float(delta.mean()),low=float(np.quantile(delta,.025)),high=float(np.quantile(delta,.975)),probability_gt_zero=float((delta>0).mean())))
    table(out / 'summary.csv', summaries)
    if contrasts:table(out / 'contrasts.csv', contrasts)
    status={'mode':'analysis_only' if source=='existing' else 'new_runs_analysis','source':source,'chains':len(chains),
            'diagnostic_gate_passed':bool(passed),'pooled_summary_written':bool(passed),'scalar_nodes':len(names),
            'thresholds':thresholds,'plots':'not_run','ppc':'not_run',
            'note':'Single-run summaries are descriptive. Failed/one-chain checks do not authorize pooled inference. This is not the historical Figure 4 report builder.'}
    (out/'status.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
    return out,status


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('valence_vz_p0_config.json'))
    parser.add_argument('--mode',choices=['check','fit','analyze','all'],default='check')
    parser.add_argument('--runs',nargs='+',help='Optional configured run IDs, e.g. r01 r02')
    parser.add_argument('--source',choices=['existing','new'],default='existing',help='Only used by analyze; all always analyzes new runs')
    parser.add_argument('--behavior-dir',type=Path,help='Completed Confirm R analysis directory; use its exported HDDM input')
    args=parser.parse_args(argv)
    config=json.loads(args.config.resolve().read_text(encoding='utf-8'))
    if args.behavior_dir:
        if args.mode == 'analyze' or (args.mode == 'check' and args.source != 'new'):
            parser.error('--behavior-dir requires fit/all, or check --source new')
        directory=args.behavior_dir.resolve()
        if not (directory/'completed.txt').is_file() or not (directory/'confirm_sample_flow.csv').is_file():
            raise ValueError('--behavior-dir must be a completed Confirm analysis, not Pilot')
        filename={'matching':'MS_match_hddm_stim.csv','identity':'MS_categ_id_hddm_stim.csv','valence':'MS_categ_val_hddm_stim.csv'}[config['task']]
        path=(directory/'hddm_input'/filename).resolve()
        for run in config['runs']:
            run['input']['path']=str(path)
    if args.runs:
        if len(set(args.runs)) != len(args.runs) or not set(args.runs).issubset({r['id'] for r in config['runs']}):
            raise ValueError('Requested run IDs must be unique and configured')
        config['runs']=[r for r in config['runs'] if r['id'] in args.runs]
    output=validate(config, existing=args.mode in ['check','analyze'] and args.source=='existing')
    if args.mode=='check':
        print(json.dumps({'config':config['name'],'runs':len(config['runs']),'files':'available','output':str(output)}));return
    if args.mode in ['fit','all']:
        for run in config['runs']:
            print('FIT_OR_RESUME',config['name'],run['id'],flush=True)
            fit_one(config,run,PROJECT,output/'runs'/run['id'])
    if args.mode in ['analyze','all']:
        out,status=analyze(config,output,'new' if args.mode=='all' else args.source)
        print(json.dumps({'output':str(out),'status':status},indent=2))


if __name__=='__main__':
    main()
