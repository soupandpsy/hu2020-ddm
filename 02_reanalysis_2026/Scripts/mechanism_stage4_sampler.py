"""Targeted exact conditional mixture update for the remaining subject-7338 modes."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import sys,time,json,ast,traceback
from pathlib import Path
import numpy as np
from scipy.special import logsumexp,gammaln
from scipy.stats import multivariate_t
import mechanism_stage3_sampler as s
m=s.m;jp=s.jp;diag=s.diag
BASE=m.ROOT/'Output/hddm/matching/continuation_stage4_20261027'

def component(x,weight,scale):
    p=x.shape[1];cov=np.cov(x,rowvar=False)*scale**2;e,U=np.linalg.eigh(cov);cov=(U*np.maximum(e,1e-6))@U.T
    invb=np.linalg.inv(cov[2:,2:]);reg=cov[:2,2:]@invb;cond=cov[:2,:2]-reg@cov[2:,:2]
    return dict(center=x.mean(axis=0),cov=cov,weight=weight,invb=invb,reg=reg,cond=cond,nu=5.)

def conditional(rest,comps):
    results=[];weights=[]
    for c in comps:
        delta=rest-c['center'][2:];db=delta@c['invb']@delta;d=len(rest);nu=c['nu']+d
        center=c['center'][:2]+c['reg']@delta;cov=(c['nu']+db)/nu*c['cond'];L=np.linalg.cholesky(cov)
        weights.append(np.log(c['weight'])+multivariate_t.logpdf(rest,loc=c['center'][2:],shape=c['cov'][2:,2:],df=c['nu']))
        results.append(dict(center=center,cov=cov,L=L,inv=np.linalg.inv(cov),nu=nu,normalizer=gammaln((nu+2)/2)-gammaln(nu/2)-np.log(nu*np.pi)-np.log(np.diag(L)).sum()))
    weights=np.exp(weights-logsumexp(weights))
    for c,w in zip(results,weights):c['weight']=w
    return results

def install(warm,model,folder):
    s.install(warm,model,folder);original=m.empirical_joint_update;active=m.allowed(model);k=len(active);si=m.subjects.index(7338)
    x=warm[:,:,2*k+4+si*k:2*k+4+(si+1)*k].reshape(-1,k);mask=(.1+np.exp(x[:,1]))>.465
    assert min(mask.sum(),(~mask).sum())>=100
    comps=[component(x[~mask],.45,1.2),component(x[mask],.45,1.2),component(x,.1,3.)]
    m.targeted_acceptance=np.zeros(3,dtype=int)
    def update(i,theta,ll,mean,sd,lc,eps,rng):
        theta,ll=original(i,theta,ll,mean,sd,lc,eps,rng)
        if i!=si:return theta,ll
        cs=conditional(theta[active][2:],comps)
        for _ in range(4):
            cand=theta.copy();cand[:2]=jp.sample(rng,cs)
            try:newll=m.subject_ll(cand,m.C@lc,eps,m.data[i])
            except (ArithmeticError,OverflowError):newll=-np.inf
            ratio=newll-ll-.5*np.sum(((cand[active]-mean[active])**2-(theta[active]-mean[active])**2)/sd[active]**2)+jp.logq(theta[:2],cs)-jp.logq(cand[:2],cs)
            m.targeted_acceptance[1]+=1
            if np.isfinite(ratio) and np.log(rng.random())<ratio:
                m.targeted_acceptance[0]+=1;m.targeted_acceptance[2]+=int((theta[1]>np.log(.365))!=(cand[1]>np.log(.365)));theta,ll=cand,newll
        return theta,ll
    m.empirical_joint_update=update
    rng=np.random.default_rng(2026279911);errors=[];balance=[]
    for row in x[rng.choice(len(x),100,replace=False)]:
        cs=conditional(row[2:],comps)
        joint=logsumexp([np.log(c['weight'])+multivariate_t.logpdf(row,loc=c['center'],shape=c['cov'],df=5) for c in comps])
        marginal=logsumexp([np.log(c['weight'])+multivariate_t.logpdf(row[2:],loc=c['center'][2:],shape=c['cov'][2:,2:],df=5) for c in comps])
        errors.append(abs(jp.logq(row[:2],cs)-(joint-marginal)))
        balance.append(jp.verify(cs)['max_detailed_balance_error'])
    assert max(errors)<1e-8 and max(balance)<1e-8
    m.write_json(folder/'targeted_proposal_validation.json',dict(subject=7338,split_t_seconds=.465,training_low=int((~mask).sum()),training_high=int(mask.sum()),weights=[.45,.45,.1],max_joint_minus_marginal_error=max(errors),max_detailed_balance_error=max(balance),updates_per_iteration=4,justification='Stage3 cross-mode transitions only 0-5 per chain despite >35% acceptance: proposal now conditions only t/a on the unchanged other individual parameters. Proposal weights do not set posterior weights.'))
    np.savez_compressed(folder/'targeted_proposal_components.npz',centers=[c['center'] for c in comps],covariances=[c['cov'] for c in comps],weights=[.45,.45,.1],df=5)
    src=''.join(m.nb['cells'][15]['source']);fn=ast.get_source_segment(src,next(n for n in ast.parse(src).body if isinstance(n,ast.FunctionDef) and n.name=='run_laplace_chain'))
    anchor='if np.log(rng.random())<ratio:theta[i]=cand;lls[i]=newll;accepted[i]+=1'
    fn=fn.replace(anchor,anchor+'\n            theta[i],lls[i]=empirical_joint_update(i,theta[i],lls[i],m,sig,lc,eps,rng)')
    fn=fn.replace("print('chain',seed,'iteration',it+1,flush=True)","print('chain',seed,'iteration',it+1,'targeted_acceptance_crossings',targeted_acceptance.tolist(),flush=True)\n            stage4_checkpoint(seed,it+1,samples,np.r_[m[active],sig[active],lc,eps,eta,theta[:,active].ravel()],rng.bit_generator.state)")
    def checkpoint(seed,iteration,samples,state,rngstate):
        target=folder/f'checkpoint_seed{seed}.npz';tmp=folder/f'checkpoint_seed{seed}.tmp.npz'
        np.savez_compressed(tmp,iteration=iteration,draws=np.asarray(samples),state=state,rng_state=json.dumps(rngstate),targeted_acceptance=m.targeted_acceptance);tmp.replace(target)
    m.stage4_checkpoint=checkpoint
    exec(compile(fn,'stage4_conditional_mixture_kernel','exec'),m.__dict__)

def run(model):
    folder=BASE/model;folder.mkdir(parents=True,exist_ok=True);status=folder/'status.json';started=time.time()
    if status.exists():raise RuntimeError('Existing run; refusing duplicate')
    m.write_json(status,dict(status='preparing',pid=os.getpid(),started=started))
    try:
        recent=np.load(s.BASE/model/'diagnostic_combined.npz')['draws'];old=np.load(s.source(model))['draws'];warm=np.concatenate([recent,old[:,::16]],axis=1)
        frame,audit=m.prepare_real_matching();m.data=m.original_prepare(frame);m.write_json(folder/'data_audit.json',audit)
        s.inspect(model,recent,folder,frame);install(warm,model,folder)
        active=m.allowed(model);k=len(active);center=recent.mean(axis=(0,1));starts=np.zeros((42,11));starts[:,active]=center[2*k+4:].reshape(42,k)
        m.lc_reference=center[2*k:2*k+2];m.eps_reference=center[2*k+2]
        m.proposals=m.build_proposals(m.data,active,starts,m.lc_reference,m.eps_reference);m.proposal_cross=m.proposal_cross_curvature(m.data,active,m.proposals,m.lc_reference,m.eps_reference)
        diag.MODEL=model;seed=2026270000+int(model[1])*100
        m.write_json(folder/'specification.json',dict(model=model,seed=seed,model_changed=False,sampler='Stage3 fixed kernel, refreshed calibration; subject7338 adds four conditional two-mode t/a MH updates',hashes=m.hashes,short_iterations=4000,short_burn=1000,formal_iterations=10000,formal_burn=2000,source=str(s.BASE/model/'diagnostic_combined.npz'),source_sha256=m.hashlib.sha256((s.BASE/model/'diagnostic_combined.npz').read_bytes()).hexdigest(),gate='All Rhat<1.01, bulk/tail ESS>=400',checkpoint='Every1000 iterations, per chain state, retained draws and RNG saved; completed-chain resume in batch'))
        m.write_json(status,dict(status='short_diagnostic_sampling',pid=os.getpid(),started=started))
        dr=diag.batch('diagnostic',seed,4000,1000,list(recent[:,-1]),folder);report=diag.summarize(dr,folder,'diagnostic');np.savez_compressed(folder/'diagnostic_combined.npz',draws=dr,labels=m.recovery_labels(model),model=model,diagnostics_pass=report['all_parameters_pass'])
        if not report['all_parameters_pass']:m.write_json(status,dict(status='short_diagnostics_failed',diagnostics=report,elapsed=time.time()-started));return
        m.write_json(status,dict(status='formal_sampling',pid=os.getpid(),started=started,short_diagnostics=report))
        dr=diag.batch('formal',seed+20,10000,2000,list(dr[:,-1]),folder);report=diag.summarize(dr,folder,'formal');np.savez_compressed(folder/'formal_fit.npz',draws=dr,labels=m.recovery_labels(model),model=model,diagnostics_pass=report['all_parameters_pass'])
        m.write_json(status,dict(status='complete' if report['all_parameters_pass'] else 'formal_diagnostics_failed',diagnostics=report,elapsed=time.time()-started))
    except Exception:m.write_json(status,dict(status='failed',error=traceback.format_exc(),elapsed=time.time()-started));raise
if __name__=='__main__':run(sys.argv[1])
