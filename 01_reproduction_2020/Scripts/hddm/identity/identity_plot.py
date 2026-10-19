"""Plot identity results only; no fitting, diagnostics recomputation or PPC."""
import argparse,csv,json,pickle,re
from pathlib import Path
from datetime import datetime
PROJECT=Path(__file__).resolve().parents[3]
def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-dir',required=True,type=Path)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args(argv);folder=args.analysis_dir.resolve()
    roots=[(PROJECT/'Output/hddm').resolve(),(PROJECT/'Output/sensitivity').resolve()]
    if not any(root in folder.parents for root in roots):raise ValueError('Analysis must be inside this project')
    config=json.loads((folder/'config.json').read_text());status=json.loads((folder/'status.json').read_text())
    if config['task']!='identity':raise ValueError('Wrong task for plotting file')
    with (folder/'sources.csv').open(encoding='utf-8-sig',newline='') as f:sources=list(csv.DictReader(f))
    for row in sources:
        if not (PROJECT/row['database']).is_file():raise FileNotFoundError('Full samples are needed for diagnostic plots: '+row['database'])
    if args.check:print('PLOT_REFERENCES_OK identity');return
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    expected=(config['sampling']['samples']-config['sampling']['burn']+config['sampling']['thin']-1)//config['sampling']['thin']
    chains=[]
    for row in sources:
        with (PROJECT/row['database']).open('rb') as f:raw=pickle.load(f)
        values={}
        for name,value in raw.items():
            if (name in ['a','z','t','v'] or re.match(r'^[vt]\(',str(name))) and isinstance(value,dict) and 0 in value:
                x=np.asarray(value[0],dtype=float)
                if x.shape!=(expected,) or not np.isfinite(x).all():raise ValueError('Invalid trace')
                values[name]=x
        chains.append(values)
    names=sorted(chains[0])
    if not names or any(set(c)!=set(names) for c in chains):raise ValueError('Parameter sets differ')
    out=folder/('plots_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));out.mkdir(exist_ok=False)
    for kind in ['traceplots','posterior']:
        with PdfPages(str(out/(kind+'.pdf'))) as pdf:
            for start in range(0,len(names),6):
                subset=names[start:start+6];fig,axes=plt.subplots(len(subset),1,figsize=(10,2.3*len(subset)),squeeze=False)
                for ax,name in zip(axes[:,0],subset):
                    if kind=='traceplots':
                        for row,c in zip(sources,chains):ax.plot(np.arange(len(c[name]))[::10],c[name][::10],label=row['run'],linewidth=.4)
                        ax.set_xlabel('Retained draw')
                    elif status['pooled_summary_written']:
                        ax.hist(np.concatenate([c[name] for c in chains]),bins=50,density=True,histtype='step',label='pooled'+str(len(chains)))
                        ax.set_ylabel('Density')
                    else:
                        for row,c in zip(sources,chains):ax.hist(c[name],bins=50,density=True,histtype='step',label=row['run'])
                        ax.set_ylabel('Density')
                    ax.set_title(name)
                axes[0,0].legend();fig.tight_layout();pdf.savefig(fig);plt.close(fig)
    print(json.dumps({'plots':str(out),'parameters':len(names)}))
if __name__=='__main__':main()
