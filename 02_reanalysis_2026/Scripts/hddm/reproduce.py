"""Reproduce the declared public outputs in a new directory."""
import argparse,os,sys,json,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('mode',choices=['plan','run','audit-data','reports','verify'],default='plan',nargs='?')
p.add_argument('--output',type=Path,required=True)
p.add_argument('--workers',type=int,choices=range(1,5),default=1)
p.add_argument('--source-repository',type=Path)
a=p.parse_args(); out=a.output.resolve()
os.environ['HU2026_OUTPUT']=str(out);os.environ['HU2026_TASK']='matching'
policy=json.loads((HERE/'report_config.json').read_text())
if any(v['model']!='m4' for v in policy['selection'].values()):raise ValueError('This reporting pipeline requires the declared M4 candidates')
commands=[]
def add(task,op,*args):
 if op in ('compare','diagnostics','ppc'):
  args=(op,*args);op='analyze'
 commands.append([sys.executable,str(HERE/task/(task+'_'+op+'.py')),*args,'--output',str(out)])
add('matching','run','preflight')
for task in policy['selection']:
 for phase in ['smoke','screening']:add(task,'run',phase,'--workers',str(a.workers))
 add(task,'compare')
 add(task,'run','final','--model',policy['selection'][task]['model'],'--workers',str(a.workers))
 add(task,'diagnostics')
add('matching','ppc');add('matching','diagnostics','--review')
for kind in ['summary','traces']:add('matching','plot','--kind',kind)
add('matching','sensitivity','--workers',str(a.workers));add('matching','sensitivity','--analyze');add('matching','recovery')
for task in ['identity','valence']:
 add(task,'run','long-refit','--workers',str(a.workers));add(task,'diagnostics','--long-refit')
if a.mode=='plan':
 print('Data audit -> commands below -> reporting -> manifest verification')
 for cmd in commands:print(subprocess.list2cmdline(cmd))
 print('Public output files:',len(json.loads((HERE/'output_manifest.json').read_text())))
 sys.exit()
import common as c
if a.mode=='run':
 if out.exists():raise FileExistsError('Use a new output directory: '+str(out))
 out.mkdir(parents=True)
 c.regenerate_data_audit(a.source_repository)
 for cmd in commands:subprocess.run(cmd,check=True)
 c.numerical_libraries();c.regenerate_reports();c.verify_output_manifest()
elif a.mode=='audit-data':c.regenerate_data_audit(a.source_repository)
elif a.mode=='verify':c.verify_output_manifest()
else:
 c.numerical_libraries();c.regenerate_reports()
