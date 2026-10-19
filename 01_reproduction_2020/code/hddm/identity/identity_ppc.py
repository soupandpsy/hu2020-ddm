"""Standalone identity PPC. Default checks references; --samples starts simulation."""
import argparse,json,os,importlib,sys
from pathlib import Path
from datetime import datetime
PROJECT=Path(__file__).resolve().parents[3]
def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--run',default='r01')
    parser.add_argument('--source',choices=['existing','new'],default='existing')
    parser.add_argument('--samples',type=int)
    args=parser.parse_args(argv);config=json.loads(args.config.resolve().read_text())
    if config['task']!='identity' or config['model'] not in ['vtz','vz']:raise ValueError('Wrong task/model')
    sys.path.insert(0,str(Path(__file__).resolve().parent/config['model']))
    workflow=importlib.import_module('identity_'+config['model']+'_run')
    run=next((r for r in config['runs'] if r['id']==args.run),None)
    if run is None:raise ValueError('Run ID is not configured')
    config['runs']=[run]
    output=workflow.resolve_output(config)
    if args.source=='existing':
        model=(PROJECT/run['model']['path']).resolve();data=(PROJECT/run['input']['path']).resolve()
        database=(PROJECT/run['database']['path']).resolve()
    else:
        d=output/'runs'/run['id'];record=json.loads((d/'complete.json').read_text())
        saved=json.loads((d/'config.json').read_text())
        if record['seed']!=run['seed'] or saved['model_kwargs']!=config['model_kwargs'] or saved['sampling']!=config['sampling']:
            raise ValueError('New results and requested model settings differ')
        model=d/'model';data=d/'input.csv';database=d/'posterior.db'
        run['input']['path']=str(data)
    for label,path in [('model',model),('input CSV',data),('sample database',database)]:
        if not path.is_file():
            raise FileNotFoundError('Missing PPC '+label+' for '+args.source+' '+run['id']+': '+str(path))
    workflow.validate(config,existing=args.source=='existing')
    if args.samples is None:print('PPC_REFERENCES_OK identity '+run['id']);return
    if args.samples<=0:raise ValueError('samples must be positive')
    import hddm,numpy as np
    if str(hddm.__version__)!=config['environment']['hddm']:raise ValueError('HDDM version mismatch')
    out=output/('ppc_'+args.source+'_'+run['id']+'_s'+str(args.samples)+'_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True,exist_ok=False);previous=Path.cwd()
    try:
        os.chdir(str(model.parent));loaded=hddm.load(str(model));observed=hddm.load_csv(str(data))
        np.random.seed(run['seed'])
        simulated=hddm.utils.post_pred_gen(loaded,samples=args.samples,progress_bar=True)
        comparison=hddm.utils.post_pred_stats(observed,simulated)
        comparison.to_csv(str(out/'ppc.csv'))
        finite=bool(np.isfinite(comparison.select_dtypes(include=[np.number]).values).all())
        (out/'status.json').write_text(json.dumps({'config':config['name'],'source':args.source,'run':run['id'],'samples':args.samples,'seed':run['seed'],'finite_numeric_summary':finite,'model':str(model),'input':str(data)},indent=2))
    finally:os.chdir(str(previous))
    print(json.dumps({'ppc':str(out),'finite_numeric_summary':finite}))
if __name__=='__main__':main()
