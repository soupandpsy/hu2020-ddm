"""Bounded, restartable process queue. Four chains per fit; at most two fitting processes."""
import os,sys,json,time,subprocess,traceback,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT)
OUT=ROOT/'Output/hddm/matching/continuation_20261026';OUT.mkdir(parents=True,exist_ok=True)
STATE=OUT/'queue_status.json'
env=os.environ.copy()
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
history=[]

def save(phase,**extra):
    obj=dict(status=phase,pid=os.getpid(),updated_epoch=time.time(),max_concurrent_fits=2,
        chains_per_fit=4,blas_threads=1,history=history,**extra)
    tmp=STATE.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2));tmp.replace(STATE)

def status(stage):
    p=OUT/(stage+'_status.json')
    if not p.exists():return {}
    return json.loads(p.read_text())

def alive(pid):
    try:os.kill(pid,0);return True
    except ProcessLookupError:return False

def fits(stages):
    pending=list(stages);running={}
    while pending or running:
        for stage in list(pending):
            s=status(stage)
            if s.get('status')=='complete':pending.remove(stage);continue
            if s.get('status')=='failed':raise RuntimeError(stage+' failed; inspect '+str(OUT/(stage+'_status.json')))
            if s.get('status')=='running' and alive(s.get('pid',-1)):
                running[stage]=None;pending.remove(stage)
        while pending and len(running)<2:
            stage=pending.pop(0);log=OUT/(stage+'_console.log')
            handle=open(log,'a',buffering=1)
            p=subprocess.Popen([sys.executable,'-u','Scripts/mechanism_continue.py',stage],stdout=handle,stderr=subprocess.STDOUT,env=env)
            running[stage]=p;history.append(dict(stage=stage,pid=p.pid,started=time.time(),log=str(log)))
        for stage,p in list(running.items()):
            s=status(stage)
            if s.get('status')=='complete':
                history.append(dict(stage=stage,finished=time.time(),exit_code=0));del running[stage]
            elif s.get('status')=='failed' or (p is not None and p.poll() is not None and p.returncode!=0):
                raise RuntimeError(stage+' failed; independent existing processes are left running')
            elif p is not None and p.poll() is not None and s.get('status')!='complete':
                raise RuntimeError(stage+' exited without a completion receipt')
            elif p is None and s.get('status')=='running' and not alive(s.get('pid',-1)):
                raise RuntimeError(stage+' process interrupted; checkpoints preserved')
        save('fitting',active=list(running),pending=pending)
        if pending or running:time.sleep(20)

def command(label,args):
    log=OUT/(label+'_console.log');save('postprocessing',active=[label])
    with log.open('a') as f:
        result=subprocess.run([sys.executable,'-u']+args,stdout=f,stderr=subprocess.STDOUT,env=env)
    history.append(dict(stage=label,finished=time.time(),exit_code=result.returncode,log=str(log)))
    if result.returncode:raise RuntimeError(label+' failed; see '+str(log))

def run():
    # Existing exact and H0 jobs are adopted by PID; never duplicated.
    fits(['exact','fit_H0','fit_H1','fit_H2','fit_H3'])
    command('corrected_recovery',['Scripts/mechanism_compare_outputs.py'])
    for model in ['H0','H1','H2','H3']:
        command('score_'+model,['Scripts/mechanism_continue.py','score_'+model])
    refits=[]
    for model in ['H0','H1','H2','H3']:
        s=json.loads((OUT/('real_'+model)/'predictive_status.json').read_text())
        refits.extend('reloo_'+model+'_'+str(i) for i in s['bad_indices'])
    # Exact held-out refits preserve every event; their cost is explicit in queue state.
    fits(refits)
    command('comparisons',['Scripts/mechanism_compare_outputs.py'])
    for model in ['H0','H1','H2','H3']:
        folder=OUT/('real_'+model)
        command('ppc_'+model,['Scripts/mechanism_real_outputs.py','--fit',str(folder/('real_'+model+'_fit.npz')),'--output',str(folder/'posterior_predictive')])
    fits(['sensitivity_0.01','sensitivity_0.025'])
    for width in ['0.01','0.025']:
        folder=OUT/('sensitivity_H3_width_'+width)
        command('sensitivity_ppc_'+width,['Scripts/mechanism_real_outputs.py','--fit',str(folder/'real_H3_fit.npz'),'--output',str(folder/'posterior_predictive'),'--width',width])
    save('complete',active=[],pending=[])
    subprocess.run([sys.executable,'Scripts/mechanism_report.py'],check=True,env=env)

if __name__=='__main__':
    # Exclusive supervisor lock, automatically recoverable after a dead PID.
    lock=OUT/'supervisor.pid'
    if lock.exists() and alive(int(lock.read_text())):raise RuntimeError('Supervisor already running')
    lock.write_text(str(os.getpid()))
    try:run()
    except Exception:
        save('failed',error=traceback.format_exc())
        subprocess.run([sys.executable,'Scripts/mechanism_report.py'],env=env)
        raise
    finally:lock.unlink(missing_ok=True)
