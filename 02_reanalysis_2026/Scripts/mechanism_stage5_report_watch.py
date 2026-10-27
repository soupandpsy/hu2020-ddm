import os,time,json,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1];base=root/'Output/hddm/matching/continuation_stage5_20261027';previous=None
while True:
    states={}
    for task in ['kfold5_repair','support_refits']:
        p=base/task/'queue_status.json'
        if p.exists():
            try:states[task]=json.loads(p.read_text())
            except json.JSONDecodeError:pass
    cv=states.get('kfold5_repair',{}).get('status');support=states.get('support_refits',{}).get('status')
    key=(cv,states.get('kfold5_repair',{}).get('phase'),len(states.get('kfold5_repair',{}).get('finished',[])),support)
    if key!=previous:
        if cv in ['complete','pilot_or_fit_diagnostics_failed','failed']:subprocess.run([sys.executable,str(root/'Scripts/mechanism_stage5_compare.py')],check=True)
        subprocess.run([sys.executable,str(root/'Scripts/mechanism_stage5_report.py')],check=True);previous=key
        (base/'report_watcher_status.json').write_text(json.dumps(dict(pid=os.getpid(),status='watching',cv=cv,support=support)))
    if cv in ['complete','pilot_or_fit_diagnostics_failed','failed'] and support in ['complete','failed']:
        (base/'report_watcher_status.json').write_text(json.dumps(dict(status='finished',cv=cv,support=support)));break
    time.sleep(30)
