"""Fold-specific frozen conditional-block proposals; exact MH, unchanged scientific model."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import sys,ast,time,json,traceback,shutil,subprocess
from pathlib import Path
import numpy as np,pandas as pd
from scipy.special import gammaln,logsumexp
from scipy.stats import multivariate_t
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import mechanism_stage4_kfold as cv
m=cv.m;s=cv.s;jp=s.jp;diag=s.diag
ROOT=m.ROOT/'Output/hddm/matching/continuation_stage5_20261027';BASE=ROOT/'kfold5_repair';BLOCK=4

def components(x):
    p=x.shape[1];t=.1+np.exp(x[:,1]);lo,hi=np.quantile(t,[.001,.999]);hist,edges=np.histogram(t,bins=160,range=(lo,hi));smooth=gaussian_filter1d(hist.astype(float),2)
    peaks,_=find_peaks(smooth,prominence=.035*smooth.max(),distance=10)
    if len(peaks)==0:peaks=np.array([np.argmax(smooth)])
    cuts=[float((edges[z]+edges[z+1])/2) for z in [a+int(np.argmin(smooth[a:b+1])) for a,b in zip(peaks[:-1],peaks[1:])]]
    # Merge a sparse proposed mode rather than estimate an unstable local covariance.
    while cuts:
        groups=np.digitize(t,cuts);sizes=np.bincount(groups,minlength=len(cuts)+1)
        small=np.where(sizes<max(100,10*p))[0]
        if len(small)==0:break
        j=int(small[0]);del cuts[min(j,len(cuts)-1)]
    groups=np.digitize(t,cuts);parts=[x[groups==i] for i in range(len(cuts)+1)];result=[]
    for vals,w,scale in [(vals,.9/len(parts),1.2) for vals in parts]+[(x,.1,3.)]:
        cov=np.cov(vals,rowvar=False)*scale**2;e,U=np.linalg.eigh(cov);cov=(U*np.maximum(e,1e-6))@U.T
        invb=np.linalg.inv(cov[BLOCK:,BLOCK:]);reg=cov[:BLOCK,BLOCK:]@invb;cond=cov[:BLOCK,:BLOCK]-reg@cov[BLOCK:,:BLOCK]
        result.append(dict(center=vals.mean(axis=0),cov=cov,weight=w,nu=5.,invb=invb,reg=reg,cond=cond))
    return result,dict(cuts_t_seconds=cuts,mode_sample_sizes=[len(x) for x in parts],weights=[c['weight'] for c in result])

def conditional(rest,comps):
    result=[];lw=[]
    for c in comps:
        delta=rest-c['center'][BLOCK:];db=delta@c['invb']@delta;d=len(rest);nu=c['nu']+d
        mean=c['center'][:BLOCK]+c['reg']@delta;cov=(c['nu']+db)/nu*c['cond'];L=np.linalg.cholesky(cov)
        lw.append(np.log(c['weight'])+multivariate_t.logpdf(rest,loc=c['center'][BLOCK:],shape=c['cov'][BLOCK:,BLOCK:],df=c['nu']))
        result.append(dict(center=mean,cov=cov,L=L,inv=np.linalg.inv(cov),nu=nu,normalizer=gammaln((nu+BLOCK)/2)-gammaln(nu/2)-.5*BLOCK*np.log(nu*np.pi)-np.log(np.diag(L)).sum()))
    weights=np.exp(lw-logsumexp(lw))
    for c,w in zip(result,weights):c['weight']=w
    return result

def install(warm,model,folder,diagnostics):
    s.install(warm,model,folder);original=m.empirical_joint_update;active=m.allowed(model);k=len(active)
    failed=diagnostics.loc[(diagnostics.r_hat>=1.01)|(diagnostics.ess_bulk<400)|(diagnostics.ess_tail<400)]
    ids=sorted(set(int(name.split('_')[1]) for name in failed.index if name.startswith(('subject_','derived_'))));si_set=[m.subjects.index(x) for x in ids]
    # Keep proven difficult 7338 in the conditioned block even if this pilot happened to pass it.
    if m.subjects.index(7338) not in si_set:si_set.append(m.subjects.index(7338))
    fitted={};summaries=[];errors=[];balance=[];rng=np.random.default_rng(2026299911)
    for si in si_set:
        x=warm[:,:,2*k+4+si*k:2*k+4+(si+1)*k].reshape(-1,k);comps,summary=components(x);fitted[si]=comps;summaries.append(dict(subject=m.subjects[si],**summary))
        for row in x[rng.choice(len(x),min(12,len(x)),replace=False)]:
            cs=conditional(row[BLOCK:],comps)
            joint=logsumexp([np.log(c['weight'])+multivariate_t.logpdf(row,loc=c['center'],shape=c['cov'],df=5) for c in comps]);marginal=logsumexp([np.log(c['weight'])+multivariate_t.logpdf(row[BLOCK:],loc=c['center'][BLOCK:],shape=c['cov'][BLOCK:,BLOCK:],df=5) for c in comps])
            errors.append(abs(jp.logq(row[:BLOCK],cs)-(joint-marginal)));balance.append(jp.verify(cs)['max_detailed_balance_error'])
        np.savez_compressed(folder/f'conditional_proposal_subject_{m.subjects[si]}.npz',centers=[c['center'] for c in comps],covariances=[c['cov'] for c in comps],weights=[c['weight'] for c in comps],df=5,block_size=BLOCK,cuts=summary['cuts_t_seconds'])
    assert max(errors)<1e-8 and max(balance)<1e-8
    m.targeted_acceptance=np.zeros((42,3),dtype=int)
    def update(i,theta,ll,mean,sd,lc,eps,rng):
        theta,ll=original(i,theta,ll,mean,sd,lc,eps,rng)
        if i not in fitted:return theta,ll
        cs=conditional(theta[active][BLOCK:],fitted[i])
        for _ in range(4):
            cand=theta.copy();cand[active[:BLOCK]]=jp.sample(rng,cs)
            try:newll=m.subject_ll(cand,m.C@lc,eps,m.data[i])
            except (ArithmeticError,OverflowError):newll=-np.inf
            ratio=newll-ll-.5*np.sum(((cand[active]-mean[active])**2-(theta[active]-mean[active])**2)/sd[active]**2)+jp.logq(theta[active[:BLOCK]],cs)-jp.logq(cand[active[:BLOCK]],cs)
            m.targeted_acceptance[i,1]+=1
            if not np.isnan(ratio) and np.log(rng.random())<ratio:
                m.targeted_acceptance[i,0]+=1
                # Occupancy transitions, not an assumption about the number of true modes.
                cuts=next(z['cuts_t_seconds'] for z in summaries if z['subject']==m.subjects[i])
                m.targeted_acceptance[i,2]+=int(np.digitize(.1+np.exp(theta[1]),cuts)!=np.digitize(.1+np.exp(cand[1]),cuts));theta,ll=cand,newll
        return theta,ll
    m.empirical_joint_update=update
    m.write_json(folder/'conditional_kernel_validation.json',dict(block_original_coordinates=active[:BLOCK].tolist(),block_names=[m.recovery_labels(model)[2*k+4+j] for j in range(BLOCK)],subjects=summaries,max_conditional_factorization_error=max(errors),max_detailed_balance_error=max(balance),updates_per_selected_subject=4,mode_rule='Smoothed t histogram, prominence3.5% maximum, separated peaks >=10/160 bins; valleys split; components <max(100,10*p) samples merged; equal proposal mode weights totaling0.9 plus broad0.1; exact MH corrects all proposal weighting.',model_changed=False))
    source=''.join(m.nb['cells'][15]['source']);fn=ast.get_source_segment(source,next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='run_laplace_chain'))
    anchor='if np.log(rng.random())<ratio:theta[i]=cand;lls[i]=newll;accepted[i]+=1';assert fn.count(anchor)==1
    fn=fn.replace(anchor,anchor+'\n            theta[i],lls[i]=empirical_joint_update(i,theta[i],lls[i],m,sig,lc,eps,rng)')
    fn=fn.replace("print('chain',seed,'iteration',it+1,flush=True)","print('chain',seed,'iteration',it+1,'targeted_accepted_attempts_crossings',targeted_acceptance.sum(axis=0).tolist(),flush=True)\n            stage5_checkpoint(seed,it+1,samples,np.r_[m[active],sig[active],lc,eps,eta,theta[:,active].ravel()],rng.bit_generator.state)")
    def checkpoint(seed,iteration,samples,state,rngstate):
        target=folder/f'checkpoint_seed{seed}.npz';temp=folder/f'checkpoint_seed{seed}.tmp.npz';np.savez_compressed(temp,iteration=iteration,draws=np.asarray(samples),state=state,rng_state=json.dumps(rngstate),targeted_acceptance=m.targeted_acceptance);temp.replace(target)
    m.stage5_checkpoint=checkpoint;exec(compile(fn,'stage5_fold_specific_conditional_block','exec'),m.__dict__)

def fit(model,fold_id):
    folder=BASE/f'{model}_fold{fold_id}';folder.mkdir(parents=True,exist_ok=True);status=folder/'status.json';start=time.time()
    if status.exists():
        if json.loads(status.read_text())['status']=='complete':return
        raise RuntimeError('Existing attempt: inspect before restart')
    m.write_json(status,dict(status='preparing',model=model,fold=fold_id,pid=os.getpid()))
    try:
        old=cv.CV/f'{model}_fold{fold_id}';warm=np.load(old/'diagnostic_combined.npz')['draws'];diagnostics=pd.read_csv(old/'diagnostic_diagnostics.csv',index_col=0);diagnostics.to_csv(folder/'previous_diagnostics.csv');failed=diagnostics.loc[(diagnostics.r_hat>=1.01)|(diagnostics.ess_bulk<400)|(diagnostics.ess_tail<400)];failed.to_csv(folder/'previous_failed_parameters.csv')
        frame,_=m.prepare_real_matching();design=pd.read_csv(cv.CV/'fold_assignments.csv');assert np.array_equal(frame.source_row,design.source_row);fold=design.fold.to_numpy();train=frame.loc[fold!=fold_id].copy();test=frame.loc[fold==fold_id].copy();test['observation']=test.index;test['test_position']=np.arange(len(test));m.data=m.original_prepare(train);assert all(x['n']==480 for x in m.data)
        s.inspect(model,warm,folder,train);install(warm,model,folder,diagnostics)
        active=m.allowed(model);k=len(active);center=warm.mean(axis=(0,1));starts=np.zeros((42,11));starts[:,active]=center[2*k+4:].reshape(42,k);m.lc_reference=center[2*k:2*k+2];m.eps_reference=center[2*k+2]
        m.proposals=m.build_proposals(m.data,active,starts,m.lc_reference,m.eps_reference);m.proposal_cross=m.proposal_cross_curvature(m.data,active,m.proposals,m.lc_reference,m.eps_reference)
        seed=2026300000+int(model[1])*10000+fold_id*100;diag.MODEL=model
        shutil.copy2(__file__,folder/'sampler_source_snapshot.py');shutil.copy2(cv.CV/'fold_assignments.csv',folder/'fold_assignments.csv')
        m.write_json(folder/'specification.json',dict(model=model,fold=fold_id,seed=seed,source=str(old/'diagnostic_combined.npz'),source_sha256=m.hashlib.sha256((old/'diagnostic_combined.npz').read_bytes()).hexdigest(),sampler_source_sha256=m.hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),model_changed=False,training_events=len(train),heldout_events=len(test),short_iterations=4000,short_burn=1000,formal_iterations=6000,formal_burn=1000,criterion='All base and derived Rhat<1.01, bulk/tail ESS>=400',proposal_calibration='Frozen training-fold pilot samples only; no posterior concatenation; new independent seeds and fresh burn-in',deadline_caveat='Same fixed full-record-derived proxy as earlier fits; conditional validation only'))
        m.write_json(status,dict(status='short_sampling',pid=os.getpid(),started=start));dr=diag.batch('diagnostic',seed,4000,1000,list(warm[:,-1]),folder);rep=diag.summarize(dr,folder,'diagnostic');np.savez_compressed(folder/'diagnostic_combined.npz',draws=dr,labels=m.recovery_labels(model),model=model,diagnostics_pass=rep['all_parameters_pass'])
        if not rep['all_parameters_pass']:m.write_json(status,dict(status='short_diagnostics_failed',diagnostics=rep,elapsed_seconds=time.time()-start));return
        m.write_json(status,dict(status='formal_sampling',pid=os.getpid(),short_diagnostics=rep));dr=diag.batch('formal',seed+20,6000,1000,list(dr[:,-1]),folder);rep=diag.summarize(dr,folder,'formal');fitpath=folder/'formal_fit.npz';np.savez_compressed(fitpath,draws=dr,labels=m.recovery_labels(model),model=model,diagnostics_pass=rep['all_parameters_pass'])
        if not rep['all_parameters_pass']:m.write_json(status,dict(status='formal_diagnostics_failed',diagnostics=rep,elapsed_seconds=time.time()-start));return
        m.write_json(status,dict(status='predicting',pid=os.getpid(),diagnostics=rep));cv.score(fitpath,model,m.data,test,folder);m.write_json(status,dict(status='complete',diagnostics=rep,elapsed_seconds=time.time()-start))
    except Exception:m.write_json(status,dict(status='failed',error=traceback.format_exc(),elapsed_seconds=time.time()-start));raise

def run_queue():
    BASE.mkdir(parents=True,exist_ok=True);state=BASE/'queue_status.json';finished=[]
    pilots=[('H0',4),('H2',4)];representatives=[('H0',0),('H2',0)];rest=[(mo,fo) for fo in range(5) for mo in m.MODELS if (mo,fo) not in pilots+representatives]
    for phase,jobs in [('severe_pilots',pilots),('representative_validation',representatives),('remaining_folds',rest)]:
        jobs=list(jobs);running=[]
        while jobs or running:
            for proc,mo,fo,log in list(running):
                if proc.poll() is not None:
                    log.close();running.remove((proc,mo,fo,log));status=json.loads((BASE/f'{mo}_fold{fo}'/'status.json').read_text());finished.append(dict(model=mo,fold=fo,exit_code=proc.returncode,status=status['status']))
            support_status=ROOT/'support_refits/queue_status.json';support_done=support_status.exists() and json.loads(support_status.read_text()).get('status')=='complete'
            capacity=2 if phase=='severe_pilots' or support_done else 1
            while jobs and len(running)<capacity:
                mo,fo=jobs.pop(0);log=open(BASE/f'{mo}_fold{fo}_console.log','a',buffering=1);proc=subprocess.Popen([sys.executable,__file__,mo,str(fo)],stdout=log,stderr=subprocess.STDOUT);running.append((proc,mo,fo,log))
            m.write_json(state,dict(status='running',phase=phase,pid=os.getpid(),capacity_models=capacity,active=[dict(pid=p.pid,model=mo,fold=fo) for p,mo,fo,_ in running],remaining_this_phase=len(jobs),finished=finished,ranking_permitted=False))
            if running:time.sleep(20)
        if phase=='severe_pilots':
            grant=ROOT/'support_refits/slot_grant.json';grant.parent.mkdir(parents=True,exist_ok=True);m.write_json(grant,dict(granted=True,max_models=1,reason='Both CV severe pilot processes exited; CV reserves at most one model until support queue complete'))
        if any(x['status']!='complete' for x in finished):m.write_json(state,dict(status='pilot_or_fit_diagnostics_failed',phase=phase,finished=finished,ranking_permitted=False));return
    m.write_json(state,dict(status='complete',finished=finished,ranking_permitted=True))
if __name__=='__main__':
    if sys.argv[1]=='queue':run_queue()
    else:fit(sys.argv[1],int(sys.argv[2]))
