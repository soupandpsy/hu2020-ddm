"""True contamination-support sensitivity refits, guarded four-chain queue."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import sys,json,time,hashlib,traceback,subprocess,contextlib
from pathlib import Path
from multiprocessing import get_context
import numpy as np
import pandas as pd
import mechanism_stage3_sampler as s
m=s.m;diag=s.diag
BASE=m.ROOT/'Output/hddm/matching/continuation_stage5_20261027/support_refits'
FITS={'H2':s.OLD/'H2_repair_v2_precision/formal_fit.npz','H3':s.BASE/'H3/formal_fit.npz'}
ORIGINAL_LL=m.subject_ll
MODEL=None

def setup_target():
    def target(theta,lam,eps,data,width=.001,nquad=7,upper=1.5):
        assert upper==1.5,'Support sensitivity requires upper=1.5 on every likelihood path'
        return ORIGINAL_LL(theta,lam,eps,data,width=width,nquad=nquad,upper=1.5)
    m.subject_ll=target

def validation(model,folder):
    setup_target();frame,audit=m.prepare_real_matching();m.data=m.original_prepare(frame);fit=FITS[model];saved=np.load(fit);assert bool(saved['diagnostics_pass']);warm=saved['draws']
    assert frame.loc[frame.choice.ge(0),'rt_new'].max()<1.5
    rows=[]
    for ch in range(4):
        th,lam,eps,eta=m.real_fit_parameters(warm[ch,-1],model)
        for i in [0,m.subjects.index(7338),41]:
            item=m.data[i];via=m.subject_ll(th[i],lam,eps,item);direct=ORIGINAL_LL(th[i],lam,eps,item,upper=1.5);old=ORIGINAL_LL(th[i],lam,eps,item,upper=2.)
            key,rt,resp,nt=next(g for g in item['groups'] if len(g[1]));a,t,z,mu=m.unpack(th[i],*key,lam)
            dnew=m.mix_density(rt,resp,a,t,z,mu,eps,upper=1.5);dold=m.mix_density(rt,resp,a,t,z,mu,eps,upper=2.)
            sn=m.mix_survival(np.array([item['proxy']]),a,t,z,mu,eps,upper=1.5)[0];so=m.mix_survival(np.array([item['proxy']]),a,t,z,mu,eps,upper=2.)[0]
            rows.append(dict(chain=ch,subject=item['sub'],wrapper_error=via-direct,changed_joint_ll=via-old,max_density_change=float(np.max(abs(dnew-dold))),survival_change=float(sn-so)))
    tab=pd.DataFrame(rows);assert tab.wrapper_error.abs().max()<1e-10;assert (tab.max_density_change>0).all();assert (tab.survival_change.abs()>0).all()
    folder.mkdir(parents=True,exist_ok=True);tab.to_csv(folder/'target_support_validation.csv',index=False);m.write_json(folder/'data_audit.json',audit)
    m.write_json(folder/'target_definition.json',dict(model=model,change='contamination Uniform RT [.1,2] to [.1,1.5]; density denominator 2*(upper-.1) and survival (upper-D)/(upper-.1) both changed',all_other_priors_and_likelihood_same=True,no_model_ranking=True,source_fit=str(fit),source_fit_sha256=hashlib.sha256(fit.read_bytes()).hexdigest(),max_observed_rt=float(frame.loc[frame.choice.ge(0),'rt_new'].max()),validation_max_wrapper_error=float(tab.wrapper_error.abs().max()),likelihood_resolution='All sampler, hyperparameter MH, frozen quadratic proposal and cross-curvature paths resolve m.subject_ll, whose wrapper calls the original with explicit upper=1.5. Original subject_ll explicitly passes upper to BOTH mix_density and mix_survival.',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),imported_source_sha256={str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(s.__file__),Path(m.__file__),Path(diag.__file__),m.NB]}))
    return warm

def kernel(warm,model,folder):
    folder.mkdir(parents=True,exist_ok=True);s.install(warm,model,folder)
    active=m.allowed(model);k=len(active);center=warm.mean(axis=(0,1));starts=np.zeros((42,11));starts[:,active]=center[2*k+4:].reshape(42,k)
    m.lc_reference=center[2*k:2*k+2];m.eps_reference=center[2*k+2]
    m.proposals=m.build_proposals(m.data,active,starts,m.lc_reference,m.eps_reference);m.proposal_cross=m.proposal_cross_curvature(m.data,active,m.proposals,m.lc_reference,m.eps_reference)
    np.savez_compressed(folder/'quadratic_proposals.npz',H=[x[0] for x in m.proposals],b=[x[1] for x in m.proposals],center=[x[2] for x in m.proposals],cross=m.proposal_cross,lc_reference=m.lc_reference,eps_reference=m.eps_reference)

def worker(args):
    ch,phase,blocks,burn,state,folder=args;folder=Path(folder);collected=[];accept=[]
    for block in range(blocks):
        seed=2026500000+int(MODEL[1])*100000+phase*10000+ch*100+block
        path=folder/f'chain{ch}_block{block}.npz'
        if path.exists():
            old=np.load(path);assert int(old['seed'])==seed;draws=old['draws'];state=draws[-1];ar=old['subject_acceptance']
        else:
            m.write_json(folder/f'chain{ch}_status.json',dict(status='sampling',pid=os.getpid(),seed=seed,block=block,blocks=blocks))
            with open(folder/f'chain{ch}.log','a',buffering=1) as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
                print('stage',phase,'seed',seed,'block',block,flush=True)
                result=m.run_laplace_chain(seed,model=MODEL,iterations=1000,burn=0,global_correction=True,dispersed=False,collapsed=True,initial_state=state)
            draws=result[0];state=draws[-1];ar=result[1]
            tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,draws=draws,seed=seed,subject_acceptance=ar,global_acceptance=result[2],collapsed_acceptance=result[3],empirical_jump_counts=m.empirical_jump_counts.copy(),ta_update_counts=m.mode_shift_counts.copy());tmp.replace(path)
        collected.append(draws);accept.append(ar)
    dr=np.concatenate(collected)[burn:];m.write_json(folder/f'chain{ch}_status.json',dict(status='complete',retained=len(dr),blocks=blocks,exit_code=0));return dr

def batch(model,phase,blocks,burn,states,folder):
    global MODEL
    MODEL=model;folder.mkdir(parents=True,exist_ok=True)
    with get_context('fork').Pool(4) as pool:result=pool.map(worker,[(ch,phase,blocks,burn,states[ch],str(folder)) for ch in range(4)])
    return np.stack(result)

def compare(fit,baseline,out,model):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=[]
    for label,path in [('original_upper_2',baseline),('refitted_upper_1.5',fit)]:
        v=np.load(path);assert bool(v['diagnostics_pass']);names=v['labels'].tolist()
        for name in ['mean_beta_v','mean_beta_dc','contamination']:
            x=v['draws'][:,:,names.index(name)].ravel() if name in names else np.zeros(v['draws'].shape[0]*v['draws'].shape[1]);lo,med,hi=np.quantile(x,[.025,.5,.975]);rows.append(dict(target=label,parameter=name,mean=float(x.mean()),median=med,lower=lo,upper=hi,fixed_zero=name not in names))
    t=pd.DataFrame(rows);t.to_csv(out/'baseline_vs_support_posterior.csv',index=False);fig,axs=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for ax,name in zip(axs,t.parameter.unique()):
        sub=t[t.parameter.eq(name)].reset_index(drop=True)
        for i,r in sub.iterrows():ax.errorbar(r['median'],i,xerr=[[r['median']-r.lower],[r.upper-r['median']]],fmt='o',capsize=4)
        ax.set_yticks([0,1]);ax.set_yticklabels(sub.target);ax.set_ylim(-.5,1.5);ax.set_title(name);ax.axvline(0,color='gray',ls='--')
    fig.suptitle(model+' true refits: within-model sensitivity; no model ranking')
    for ext in ['png','pdf']:fig.savefig(out/('baseline_vs_support_forest.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)

def run(model):
    folder=BASE/model;folder.mkdir(parents=True,exist_ok=True);status=folder/'status.json';lock=folder/'RUNNING.lock'
    if status.exists() and json.loads(status.read_text()).get('status') in ['complete','short_diagnostics_failed','formal_diagnostics_failed']:return
    if lock.exists():
        pid=int(lock.read_text())
        if Path('/proc/'+str(pid)).exists():raise RuntimeError('Existing live fitting process')
    lock.write_text(str(os.getpid()));started=time.time()
    try:
        warm=validation(model,folder);m.write_json(status,dict(status='preparing_calibration',pid=os.getpid()))
        m.write_json(folder/'specification.json',dict(model=model,chains=4,support_upper=1.5,deadline_width=.001,calibration_iterations=2000,calibration_burn=500,diagnostic_iterations=4000,diagnostic_burn=1000,formal_iterations=10000,formal_burn=2000,seed_formula='2026500000+model_digit*100000+phase*10000+chain*100+block',phases={'calibration':0,'diagnostic':1,'formal':2},checkpoint='Each 1000-iteration block saves all samples and end-state atomically. Unique deterministic block seed permits exact completed-block restart; incomplete block is rerun with same seed. No samples from previous target enter inference.',adaptation='Initial proposal uses old target only for calibration. Frozen final proposal is retrained on changed-target calibration plus 5% historical subsample for mode coverage. Formal kernel is fixed.',gate='All raw and derived Rhat<1.01, bulk/tail ESS>=400; no automatic extension on failure',prior_unchanged=True))
        kernel(warm,model,folder/'calibration_kernel');cal=batch(model,0,2,500,list(warm[:,-1]),folder/'calibration');np.savez_compressed(folder/'calibration_combined.npz',draws=cal,model=model,for_proposals_only=True)
        train=np.concatenate([cal,warm[:,::max(1,warm.shape[1]//75)]],axis=1);kernel(train,model,folder/'frozen_kernel');diag.MODEL=model
        m.write_json(status,dict(status='short_diagnostic_sampling',pid=os.getpid()))
        dr=batch(model,1,4,1000,list(cal[:,-1]),folder/'diagnostic');report=diag.summarize(dr,folder,'diagnostic');np.savez_compressed(folder/'diagnostic_combined.npz',draws=dr,labels=m.recovery_labels(model),model=model,diagnostics_pass=report['all_parameters_pass'])
        if not report['all_parameters_pass']:m.write_json(status,dict(status='short_diagnostics_failed',diagnostics=report,elapsed_seconds=time.time()-started));return
        m.write_json(status,dict(status='formal_sampling',pid=os.getpid(),short_diagnostics=report))
        dr=batch(model,2,10,2000,list(dr[:,-1]),folder/'formal');report=diag.summarize(dr,folder,'formal');fit=folder/'formal_fit.npz';np.savez_compressed(fit,draws=dr,labels=m.recovery_labels(model),model=model,diagnostics_pass=report['all_parameters_pass'],contamination_upper=1.5,deadline_width=.001,data_kind='observed')
        if report['all_parameters_pass']:compare(fit,FITS[model],folder,model)
        m.write_json(status,dict(status='complete' if report['all_parameters_pass'] else 'formal_diagnostics_failed',diagnostics=report,elapsed_seconds=time.time()-started,ppc='Not generated: must use upper1.5 in BOTH contamination generation and predictive deadlines. Original-upper PPC cannot be reused.',model_ranking=False))
    except Exception:
        m.write_json(status,dict(status='failed',error=traceback.format_exc(),elapsed_seconds=time.time()-started));raise
    finally:lock.unlink(missing_ok=True)

def wait_and_run():
    BASE.mkdir(parents=True,exist_ok=True);lock=BASE/'waiter.lock'
    if lock.exists() and Path('/proc/'+lock.read_text().strip()).exists():raise RuntimeError('Existing live waiter')
    lock.write_text(str(os.getpid()));state=BASE/'queue_status.json'
    try:
        while True:
            grant=BASE/'slot_grant.json'
            allowed=grant.exists() and json.loads(grant.read_text()).get('granted') is True
            m.write_json(state,dict(status='waiting_for_explicit_four_chain_slot',pid=os.getpid(),grant=str(grant),models=['H2','H3'],max_concurrent_models=1))
            if allowed:break
            time.sleep(30)
        finished=[]
        for model in ['H2','H3']:
            m.write_json(state,dict(status='running',pid=os.getpid(),model=model,finished=finished,max_concurrent_models=1))
            with open(BASE/(model+'_console.log'),'a',buffering=1) as log:r=subprocess.run([sys.executable,str(Path(__file__).resolve()),'fit',model],stdout=log,stderr=subprocess.STDOUT)
            finished.append(dict(model=model,exit_code=r.returncode,status=json.loads((BASE/model/'status.json').read_text()).get('status') if (BASE/model/'status.json').exists() else 'missing'))
        m.write_json(state,dict(status='complete',review_required=True,finished=finished,slot_released=True,no_model_ranking=True))
    finally:lock.unlink(missing_ok=True)

if __name__=='__main__':
    BASE.mkdir(parents=True,exist_ok=True)
    if sys.argv[1]=='validate':
        for model in ['H2','H3']:validation(model,BASE/model)
        m.write_json(BASE/'ready.json',dict(status='validated_ready',command='python Scripts/mechanism_stage5_support_refit.py wait',slot_requirement='Explicit slot_grant.json granted=true: one model with four chains, H2 then H3 sequentially',target='true refits with upper1.5',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    elif sys.argv[1]=='wait':wait_and_run()
    elif sys.argv[1]=='fit':run(sys.argv[2])

