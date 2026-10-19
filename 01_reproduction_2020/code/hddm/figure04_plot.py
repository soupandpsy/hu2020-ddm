"""Rebuild Figure 4 panels and related extended plots from saved r01 traces.

No fitting or PPC. Uses the historical Kabuki 50-bin/300-point spline method.
Creates a fresh output directory; historical results are never overwritten.
"""
import argparse
import csv
from datetime import datetime
import importlib
import sys
import json
from pathlib import Path
import pickle
import shlex

PROJECT = Path(__file__).resolve().parents[2]
LABELS = ['Self.Good', 'Self.Bad', 'Other.Good', 'Other.Bad']
LEGEND = ['Good-self', 'Bad-self', 'Good-other', 'Bad-other']
COLORS = ['#1f77b4', '#2ca02c', '#d62728', '#17becf']
PURPOSE = {'v': 'drift', 't': 'nondecision', 'z': 'starting_bias'}


def specifications():
    specs = []
    for task in ['matching', 'identity', 'valence']:
        for param in ['v', 't', 'z']:
            conditions = ['match', 'mismatch'] if task == 'matching' and param != 'z' else ['shared']
            for condition in conditions:
                nodes = ['z'] if param == 'z' else [param+'('+lab.replace('.', '.'+condition+'.')+')' if task == 'matching' else param+'('+lab+')' for lab in LABELS]
                panel = {'matching': 'A', 'identity': 'B', 'valence': 'C'}[task] if param == 'v' and condition != 'mismatch' else ''
                stem = task+'_vtz_p0_r01_'+(condition+'_' if condition != 'shared' else '')+PURPOSE[param]
                specs.append(dict(task=task,parameter=param,condition=condition,nodes=nodes,panel=panel,stem=stem))
    return specs


def table(path, rows):
    with path.open('x', encoding='utf-8-sig', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',choices=['all','matching','identity','valence'],default='all')
    parser.add_argument('--check',action='store_true',help='Read samples without writing figures')
    parser.add_argument('--source',choices=['bundled','existing','new'],default='bundled')
    for task in ['matching','identity','valence']:
        parser.add_argument('--'+task+'-analysis',type=Path,help='New '+task+' vtz analysis directory containing sources.csv')
    args=parser.parse_args(argv)
    if args.source != 'new' and any(getattr(args,t+'_analysis') for t in ['matching','identity','valence']):
        parser.error('Analysis directories require --source new')
    import numpy as np
    from scipy.interpolate import InterpolatedUnivariateSpline
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    specs=[s for s in specifications() if args.task=='all' or s['task']==args.task]
    tasks=sorted({s['task'] for s in specs})
    traces={};provenance=[]
    for task in tasks:
        if args.source=='bundled':
            records=json.loads((PROJECT/'data/posterior/manifest.json').read_text())
            r=next(r for r in records if r['task']==task)
            path=PROJECT/'data/posterior'/r['file']
            with np.load(str(path),allow_pickle=False) as raw:
                traces[task]={n:raw[n].copy() for n in r['parameters']}
            if any(x.shape!=(r['draws'],) or not np.isfinite(x).all() for x in traces[task].values()):raise ValueError('Invalid bundled samples')
            provenance.append(dict(task=task,source='bundled',database='data/posterior/'+r['file'],run=r['run'],seed=r['seed'],draws=r['draws']))
            continue
        model_dir=Path(__file__).resolve().parent/task/'vtz'
        sys.path.insert(0,str(model_dir))
        module=importlib.import_module(task+'_vtz_run')
        if args.source=='new':
            analysis=getattr(args,task+'_analysis')
            if analysis is None:
                parser.error('--source new requires --'+task+'-analysis for the selected task')
            analysis=analysis.resolve()
            config_path=analysis/'config.json'
            config=json.loads(config_path.read_text())
            status=json.loads((analysis/'status.json').read_text())
            if status.get('source')!='new' or config['task']!=task or config['model']!='vtz' or config['model_kwargs']['p_outlier']!=0:
                raise ValueError('Expected a new '+task+' vtz p0 analysis')
            with (analysis/'sources.csv').open(encoding='utf-8-sig',newline='') as stream:
                records=list(csv.DictReader(stream))
            record=next((r for r in records if r['run']=='r01'),None)
            if record is None:
                raise ValueError('Figure 4 requires r01 in the selected analysis')
            path=(PROJECT/record['database']).resolve()
            run=next(r for r in config['runs'] if r['id']=='r01')
            if int(record['seed'])!=run['seed']:
                raise ValueError('Analysis source seed differs from configuration')
        else:
            config_path=model_dir/(task+'_vtz_p0_config.json')
            config=json.loads(config_path.read_text())
            run=next(r for r in config['runs'] if r['id']=='r01')
            config['runs']=[run];module.validate(config)
            path=(PROJECT/run['database']['path']).resolve()
        if not path.is_file():
            raise FileNotFoundError('Selected full samples are missing: '+str(path))
        with path.open('rb') as f:raw=pickle.load(f)
        expected=(config['sampling']['samples']-config['sampling']['burn']+config['sampling']['thin']-1)//config['sampling']['thin']
        traces[task]={}
        for node in {n for s in specs if s['task']==task for n in s['nodes']}:
            x=np.asarray(raw[node][0],dtype=float)
            if x.shape!=(expected,) or not np.isfinite(x).all():raise ValueError('Invalid trace: '+node)
            traces[task][node]=x
        provenance.append(dict(task=task,source=args.source,config=str(config_path),database=str(path),run='r01',seed=run['seed'],draws=expected))
    plot_data={}
    for s in specs:
        arrays=traces[s['task']];lo=min(arrays[n].min() for n in s['nodes']);hi=max(arrays[n].max() for n in s['nodes'])
        grid=np.linspace(lo,hi,300);ys=[]
        for node in s['nodes']:
            hist=np.histogram(arrays[node],bins=50,range=(lo,hi),density=True)[0]
            y=InterpolatedUnivariateSpline(np.linspace(lo,hi,50),hist)(grid)
            ys.append(y)
        plot_data[s['stem']]=(grid,ys)
    if args.check:
        print(json.dumps(dict(figures=len(specs),parameters=sum(len(v) for v in traces.values()),samples='read successfully',writes=False)));return
    out=PROJECT/'results'/'figures'/('figure04_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True,exist_ok=False)
    def draw(ax,s):
        grid,ys=plot_data[s['stem']]
        for i,y in enumerate(ys):ax.plot(grid,y,color=COLORS[i],lw=2,label='z' if s['parameter']=='z' else LEGEND[i])
        if s['parameter']=='z':ax.axvline(.5,color='gray',linestyle='--',lw=1)
        prefix='Figure 4'+s['panel'] if s['panel'] else 'Extended'
        condition=' / '+s['condition'] if s['condition']!='shared' else ''
        ax.set_title(prefix+' / '+s['task'].title()+condition+'\nvtz, p_outlier=0, r01')
        ax.set_xlabel({'v':'Drift rate (v)','t':'Nondecision time (t, seconds)','z':'Starting-point proportion (z)'}[s['parameter']])
        ax.set_ylabel('Posterior density');ax.legend(loc='best',fontsize=8)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    mapping=[]
    def entry_for(task):
        command=['python','-B','code/hddm/figure04_plot.py','--task',task,'--source',args.source]
        if args.source=='new':
            command+=['--'+task+'-analysis',str(getattr(args,task+'_analysis').resolve())]
        return ' '.join(shlex.quote(x) for x in command)
    for s in specs:
        fig,ax=plt.subplots(figsize=(7,4.8));draw(ax,s);fig.tight_layout()
        for ext in ['pdf','png']:fig.savefig(str(out/(s['stem']+'.'+ext)),dpi=180)
        plt.close(fig)
        mapping.append(dict(paper_panel='Figure 4'+s['panel'] if s['panel'] else 'Extended (not Figure 4)',task=s['task'],model='vtz',p_outlier=0,run='r01',condition=s['condition'],parameter=s['parameter'],nodes='; '.join(s['nodes']),entry=entry_for(s['task']),pdf=s['stem']+'.pdf',png=s['stem']+'.png'))
    if args.task=='all':
        fig,axes=plt.subplots(1,3,figsize=(17,4.8))
        for ax,s in zip(axes,[s for s in specs if s['panel']]):draw(ax,s)
        fig.tight_layout()
        for ext in ['pdf','png']:fig.savefig(str(out/('figure04_vtz_p0_r01_drift.'+ext)),dpi=180)
        plt.close(fig)
    table(out/'figure_index.csv',mapping)
    (out/'sources.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print(json.dumps(dict(output=str(out),figures=len(specs))))



if __name__=='__main__':main()
