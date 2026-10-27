"""Stratified five-fold same-participant prediction after all full-data fits pass.
No PSIS weights are used for K-fold held-out predictions. Historical files read-only.
"""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import sys,json,time,traceback,subprocess,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import mechanism_stage3_sampler as s
m=s.m;diag=s.diag
BASE=m.ROOT/'Output/hddm/matching/continuation_stage4_20261027'
CV=BASE/'kfold5'
FITS={'H0':BASE/'H0/formal_fit.npz','H1':BASE/'H1/formal_fit.npz','H2':s.OLD/'H2_repair_v2_precision/formal_fit.npz','H3':s.BASE/'H3/formal_fit.npz'}

def make_design(frame):
    CV.mkdir(parents=True,exist_ok=True);path=CV/'fold_assignments.csv'
    rng=np.random.default_rng(2026278000);fold=np.full(len(frame),-1);offsets={}
    for key,g in frame.groupby(['Subject','condition','stim','relation'],sort=True):
        ii=rng.permutation(g.index.to_numpy());offset=offsets.get(key[0],0);fold[ii]=(np.arange(len(ii))+offset)%5;offsets[key[0]]=offset+len(ii)
    design=frame[['Subject','source_row','condition','stim','relation']].copy();design['fold']=fold
    if path.exists():pd.testing.assert_frame_equal(pd.read_csv(path),design.reset_index(drop=True),check_dtype=False)
    else:design.to_csv(path,index=False)
    assert (fold>=0).all() and all((fold==j).sum()==5040 for j in range(5))
    if not (CV/'specification.json').exists():m.write_json(CV/'specification.json',dict(folds=5,seed=2026278000,strata=['Subject','condition','stim','relation'],observations=len(frame),holdout='120 events per participant in each fold; event-marginal prediction conditional only on the other four training folds and fixed deadline proxy',same_folds_all_models=True,no_response_used_to_assign_fold=True,deadline_caveat='Proxy itself is derived from full raw records and fixed. This remains conditional prediction, not leakage-free prospective validation.',likelihood='Shared participant deadline integrated with training observations only; each held-out event added singly, never conditioned on other held-out responses.',uncertainty='Paired model score differences aggregated over 42 participant clusters. Five folds are not five independent datasets.',sampler='Same priors and likelihood; new four-chain short gate then independent formal segment',short_iterations=4000,short_burn=1000,formal_iterations=6000,formal_burn=1000,gate='All base and derived Rhat<1.01 and bulk/tail ESS>=400',reason='Full-data H2 has 17 and H3 has19 high Pareto-k events; uniform K-fold replaces unreliable PSIS for all models, without deleting observations.'))
    return fold

def heldout_loglik(theta,lam,eps,eta,items,test):
    result=np.full(len(test),np.nan);gx,gw=m.quadrature_nodes(7)
    for si,(th,item) in enumerate(zip(theta,items)):
        sub=m.subjects[si];gsub=test.loc[test.Subject.eq(sub)];L=item['proxy']-.001;U=item['proxy']+.001
        def integral(lower,extra=None):
            lo=max(L,lower)
            if lo>=U:return -np.inf
            ds=lo+(gx+1)*(U-lo)/2;lw=np.log(gw/2)+np.log((U-lo)/(U-L))
            for key,rt,resp,nt in item['groups']:
                if nt:
                    a,t,z,mu=m.unpack(th,*key,lam);lw+=nt*np.log(m.mix_survival(ds,a,t,z,mu,eps))
            if extra is not None:
                a,t,z,mu=m.unpack(th,*extra,lam);lw+=np.log(m.mix_survival(ds,a,t,z,mu,eps))
            return m.logsumexp(lw)
        base=integral(item['maxrt'])
        for key,g in gsub.groupby(['condition','stim','relation'],sort=True):
            legal=g.choice.ge(0);rt=g.loc[legal,'rt_new'].to_numpy();idx=g.loc[legal,'test_position'].to_numpy()
            a,t,z,mu=m.unpack(th,*key,lam)
            value=np.log1p(-eta)+np.log(m.mix_density(rt,g.loc[legal,'choice'].to_numpy(),a,t,z,mu,eps))
            for j,rr in enumerate(rt):
                if rr>max(L,item['maxrt']):value[j]+=integral(rr)-base
            result[idx]=value
            result[g.loc[g.choice.eq(-1),'test_position'].to_numpy()]=np.log1p(-eta)+integral(item['maxrt'],key)-base
            result[g.loc[g.choice.eq(-2),'test_position'].to_numpy()]=np.log(eta)
    assert np.isfinite(result).all()
    return result

def score(fit,model,items,test,out):
    saved=np.load(fit);assert bool(saved['diagnostics_pass']);dr=saved['draws'];ids=np.linspace(0,dr.shape[1]-1,800,dtype=int);selected=dr[:,ids]
    path=out/'heldout_pointwise_loglik.npy';ll=np.lib.format.open_memmap(path,mode='w+',dtype='float64',shape=(4,800,len(test)))
    for ch in range(4):
        for it,row in enumerate(selected[ch]):
            theta,lam,eps,eta=m.real_fit_parameters(row,model);ll[ch,it]=heldout_loglik(theta,lam,eps,eta,items,test)
        ll.flush();m.write_json(out/'prediction_progress.json',dict(completed_chains=ch+1,total_chains=4))
    tab=test[['observation','Subject','source_row','condition','stim','relation','choice','rt_new']].copy();tab['elpd']=m.logsumexp(ll.reshape(3200,-1),axis=0)-np.log(3200)
    tab.to_csv(out/'heldout_scores.csv',index=False);np.savez_compressed(out/'posterior_selection.npz',indices=ids)

def fit_fold(model,fold_id):
    out=CV/f'{model}_fold{fold_id}';out.mkdir(parents=True,exist_ok=True);status=out/'status.json';start=time.time()
    if status.exists():
        old=json.loads(status.read_text())
        if old['status']=='complete':return
        raise RuntimeError('Prior attempt exists; inspect rather than duplicate')
    m.write_json(status,dict(status='preparing',pid=os.getpid(),model=model,fold=fold_id))
    try:
        frame,_=m.prepare_real_matching();fold=make_design(frame);train=frame.loc[fold!=fold_id].copy();test=frame.loc[fold==fold_id].copy();test['observation']=test.index;test['test_position']=np.arange(len(test))
        assert train.groupby('Subject').size().eq(480).all() and test.groupby('Subject').size().eq(120).all()
        m.data=m.original_prepare(train);items=m.data;warm=np.load(FITS[model])['draws'];active=m.allowed(model);k=len(active)
        assert bool(np.load(FITS[model])['diagnostics_pass'])
        if model in ['H0','H1']:
            import mechanism_stage4_sampler as s4
            s4.install(warm,model,out)
        else:s.install(warm,model,out)
        center=warm.mean(axis=(0,1));starts=np.zeros((42,11));starts[:,active]=center[2*k+4:].reshape(42,k)
        m.lc_reference=center[2*k:2*k+2];m.eps_reference=center[2*k+2]
        m.proposals=m.build_proposals(m.data,active,starts,m.lc_reference,m.eps_reference);m.proposal_cross=m.proposal_cross_curvature(m.data,active,m.proposals,m.lc_reference,m.eps_reference)
        # Validate this NEW holdout calculation at a few training+single-event ratios, not historical grids.
        checks=[];picked=list(test.groupby('choice',sort=True).head(2).index)
        for si,item in enumerate(items):
            picked+=test.index[test.Subject.eq(m.subjects[si]) & test.choice.ge(0) & test.rt_new.gt(max(item['proxy']-.001,item['maxrt']))].tolist()
        picked=list(dict.fromkeys(picked))
        for state_id,row in enumerate([warm[0,-1],warm[3,0]]):
            th,lam,eps,eta=m.real_fit_parameters(row,model);v=heldout_loglik(th,lam,eps,eta,items,test)
            for idx in picked:
                ev=test.loc[idx];si=m.subjects.index(int(ev.Subject));plus=m.original_prepare(pd.concat([train,frame.loc[[idx]]]).sort_index())
                direct=m.subject_ll(th[si],lam,eps,plus[si])-m.subject_ll(th[si],lam,eps,items[si])+(np.log(eta) if ev.choice==-2 else np.log1p(-eta))
                checks.append(dict(state=state_id,observation=int(idx),choice=int(ev.choice),deadline_lower_bound_changes=bool(ev.choice>=0 and ev.rt_new>max(items[si]['proxy']-.001,items[si]['maxrt'])),conditional=float(v[int(ev.test_position)]),joint_ratio=float(direct),error=float(v[int(ev.test_position)]-direct)))
        pd.DataFrame(checks).to_csv(out/'new_holdout_ratio_validation.csv',index=False);assert max(abs(x['error']) for x in checks)<1e-8
        train[['Subject','source_row']].to_csv(out/'training_rows.csv',index=True);test[['observation','Subject','source_row']].to_csv(out/'heldout_rows.csv',index=False)
        diag.MODEL=model;seed=2026280000+int(model[1])*10000+fold_id*100
        m.write_json(out/'fit_specification.json',dict(model=model,fold=fold_id,seed=seed,training_events=len(train),heldout_events=len(test),source_fit=str(FITS[model]),source_fit_sha256=m.hashlib.sha256(FITS[model].read_bytes()).hexdigest(),model_changed=False,initialization='Full posterior used only for initial states and fixed proposal calibration, not as a prior or retained training posterior. Exact MH uses training rows only.'))
        m.write_json(status,dict(status='short_sampling',pid=os.getpid(),started=start))
        dr=diag.batch('diagnostic',seed,4000,1000,list(warm[:,-1]),out);rep=diag.summarize(dr,out,'diagnostic');np.savez_compressed(out/'diagnostic_combined.npz',draws=dr,model=model,labels=m.recovery_labels(model),diagnostics_pass=rep['all_parameters_pass'])
        if not rep['all_parameters_pass']:m.write_json(status,dict(status='short_diagnostics_failed',diagnostics=rep));return
        m.write_json(status,dict(status='formal_sampling',pid=os.getpid(),short_diagnostics=rep))
        dr=diag.batch('formal',seed+20,6000,1000,list(dr[:,-1]),out);rep=diag.summarize(dr,out,'formal');fit=out/'formal_fit.npz';np.savez_compressed(fit,draws=dr,model=model,labels=m.recovery_labels(model),diagnostics_pass=rep['all_parameters_pass'])
        if not rep['all_parameters_pass']:m.write_json(status,dict(status='formal_diagnostics_failed',diagnostics=rep));return
        m.write_json(status,dict(status='predicting',pid=os.getpid(),diagnostics=rep));score(fit,model,items,test,out)
        m.write_json(status,dict(status='complete',diagnostics=rep,elapsed_seconds=time.time()-start))
    except Exception:m.write_json(status,dict(status='failed',error=traceback.format_exc(),elapsed_seconds=time.time()-start));raise

def queue():
    CV.mkdir(parents=True,exist_ok=True);state=CV/'queue_status.json'
    m.write_json(state,dict(status='waiting_for_full_data_convergence',pid=os.getpid()))
    while True:
        st=[json.loads((BASE/model/'status.json').read_text()) for model in ['H0','H1']]
        if any(x['status'] in ['short_diagnostics_failed','formal_diagnostics_failed','failed'] for x in st):
            m.write_json(state,dict(status='blocked_by_full_data_convergence',models=st,ranking_permitted=False));return
        if all(x['status']=='complete' for x in st):break
        time.sleep(20)
    for model,path in FITS.items():assert path.exists() and bool(np.load(path)['diagnostics_pass']),model
    frame,_=m.prepare_real_matching();make_design(frame)
    jobs=[(model,fold) for fold in range(5) for model in m.MODELS];running=[];finished=[]
    while jobs or running:
        for obj in list(running):
            proc,model,fold,log=obj
            if proc.poll() is not None:
                log.close();running.remove(obj);finished.append(dict(model=model,fold=fold,exit_code=proc.returncode))
        while jobs and len(running)<2:
            model,fold=jobs.pop(0);log=open(CV/f'{model}_fold{fold}_console.log','a',buffering=1);p=subprocess.Popen([sys.executable,__file__,model,str(fold)],stdout=log,stderr=subprocess.STDOUT);running.append((p,model,fold,log))
        m.write_json(state,dict(status='running' if running else 'fits_finished_review_required',pid=os.getpid(),active=[dict(pid=p.pid,model=mo,fold=fo) for p,mo,fo,_ in running],remaining=len(jobs),finished=finished,ranking_permitted=False))
        if running:time.sleep(20)
if __name__=='__main__':
    if len(sys.argv)==2 and sys.argv[1]=='queue':queue()
    elif len(sys.argv)==3:fit_fold(sys.argv[1],int(sys.argv[2]))
    else:raise SystemExit('usage: queue OR MODEL FOLD')



