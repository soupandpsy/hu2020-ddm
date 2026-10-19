"""valence_vz model: one fit. Called by valence_vz_run.py; no sampling on import."""
import contextlib
import json
import os
from pathlib import Path
import shutil
import time


def fit_one(config, run, project, output):
    """Fit one configured run in a new directory; never overwrite an existing run."""
    if config['task'] != 'valence' or config['model'] != 'vz' or config['model_kwargs']['depends_on'] != {'v': ['val', 'id']}:
        raise ValueError('Configuration does not match this model file')
    import hddm
    import numpy as np
    if str(hddm.__version__) != config['environment']['hddm']:
        raise RuntimeError('HDDM version differs from config')
    output = Path(output).resolve()
    done = output / 'complete.json'
    if output.exists():
        raise FileExistsError('Run directory already exists. Select another configured run or output_root: '+str(output))
    data_path = (project / run['input']['path']).resolve()
    data = hddm.load_csv(str(data_path))
    if len(data) != config['data_checks']['rows'] or data.subj_idx.nunique() != config['data_checks']['subjects']:
        raise RuntimeError('Data counts differ from config')
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(str(data_path), str(output / 'input.csv'))
    (output / 'config.json').write_text(json.dumps(config, indent=2), encoding='utf-8')
    previous = Path.cwd()
    started = time.time()
    try:
        os.chdir(str(output))
        with (output / 'fit.log').open('x') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            np.random.seed(run['seed'])
            print('seed', run['seed'], 'HDDM', hddm.__version__, flush=True)
            model = hddm.HDDMStimCoding(data, include='z', stim_col='stim',
                depends_on={'v': ['val', 'id']}, split_param='v', drift_criterion=False,
                p_outlier=config['model_kwargs']['p_outlier'])
            model.find_starting_values()
            model.sample(config['sampling']['samples'], burn=config['sampling']['burn'],
                         thin=config['sampling']['thin'], db='pickle', dbname='posterior.db')
            model.save('model')
            record = {'run_id': run['id'], 'seed': run['seed'],
                      'dic': float(model.dic), 'elapsed_seconds': time.time() - started}
        with done.open('x') as stream:
            json.dump(record, stream, indent=2)
        return record
    finally:
        os.chdir(str(previous))
