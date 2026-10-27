import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import sys,ast,json,time,contextlib,traceback
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import mechanism_continue as m
import mechanism_joint_proposals as jp
import mechanism_h2h3_diagnostic as diag
BASE=m.ROOT/'Output/hddm/matching/continuation_stage3_20261026'
OLD=m.ROOT/'Output/hddm/matching/diagnostic_H2H3_20261026_1425'

def source(model):
    return OLD/'H3_repair_v2/formal_fit.npz' if model=='H3' else m.BASE/f'real_{model}/real_{model}_fit.npz'

def inspect(model,warm,folder,frame):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    table=m.recovery_diagnostics(warm,model);table.to_csv(folder/'historical_diagnostics.csv',index_label='parameter')
    labels=m.recovery_labels(model); selected=table.loc[labels].sort_values('r_hat',ascending=False).head(6).index.tolist()
    fig,axs=plt.subplots(len(selected),3,figsize=(15,2*len(selected)),constrained_layout=True)
    for r,name in enumerate(selected):
        j=labels.index(name);x=warm[:,:,j];rank=rankdata(x).reshape(x.shape)
        for c in range(4):
            axs[r,0].plot(x[c],lw=.3,alpha=.6,label=str(c))
            axs[r,1].hist(x[c],bins=50,histtype='step',density=True)
            axs[r,2].hist(rank[c],bins=20,histtype='step',density=True)
        axs[r,0].set_ylabel(name)
    fig.suptitle(model+' failed historical chains: trace / density / pooled ranks; diagnostic only')
    fig.savefig(folder/'historical_trace_density_rank.png',dpi=150);plt.close(fig)
    subs=list(dict.fromkeys(int(p.split('_')[1]) for p in selected if p.startswith('subject_')))
    if model=='H3' and 7338 not in subs:subs.append(7338)
    fig,axs=plt.subplots(len(subs),2,figsize=(12,4*len(subs)),squeeze=False,constrained_layout=True)
    rows=[]
    for r,sub in enumerate(subs):
        a=np.exp(warm[:,:,labels.index(f'subject_{sub}_log_a')]);t=.1+np.exp(warm[:,:,labels.index(f'subject_{sub}_log_t_excess')])
        for c in range(4):axs[r,0].scatter(a[c,::20],t[c,::20],s=2,alpha=.25,label=f'chain {c}')
        axs[r,0].set(xlabel='a',ylabel='t (s)',title=f'{sub} joint sampled states')
        f=frame.loc[frame.Subject.eq(sub)];legal=f.choice.ge(0)
        for acc in [0,1]:axs[r,1].hist(f.loc[legal & f.ACC.eq(acc),'rt_new'],bins=50,alpha=.5,label=f'ACC {acc}')
        axs[r,1].set(xlabel='Observed legal RT (s)',title=f'{sub} RT');axs[r,1].legend()
        rows.append(dict(model=model,subject=sub,n=len(f),legal=int(legal.sum()),timeout=int(f.choice.eq(-1).sum()),illegal=int(f.choice.eq(-2).sum()),accuracy_legal=float(f.loc[legal,'ACC'].mean()),min_legal_rt=float(f.loc[legal,'rt_new'].min()),rt_q01=float(f.loc[legal,'rt_new'].quantile(.01)),t_chain_means=t.mean(axis=1).tolist(),a_t_correlation=np.corrcoef(a.ravel(),t.ravel())[0,1]))
        f.to_csv(folder/f'subject_{sub}_observations.csv',index=False)
    fig.savefig(folder/'historical_t_a_and_observed_rt.png',dpi=150);plt.close(fig)
    pd.DataFrame(rows).to_csv(folder/'problem_subject_summary.csv',index=False)
    m.write_json(folder/'diagnostic_interpretation.json',dict(subjects=rows,confirmed='Historical chain locations disagree in t/a. Original t=.1+exp(log_t_excess) is not hard-bounded by minimum RT: contamination mixture allows responses faster than t. Different t regions can assign short responses to contamination, giving separated likelihood regions.',hypothesis='Sparse frozen proposal mixture and high-dimensional independent proposals may undersample narrow or newly visited t regions. New proposal trains on longest historical chains and adds t/a-only multiscale symmetric updates.',constraint='No participant removed; no likelihood, prior, or scientific parameter changed. Historical chains used only for proposal calibration and diagnostics.'))

def install(warm,model,folder):
    active=m.allowed(model);k=len(active);flat=warm.reshape(-1,warm.shape[-1]);mixtures=[];covta=[]
    for i in range(42):
        x=flat[:,2*k+4+i*k:2*k+4+(i+1)*k];groups=np.array_split(np.argsort(x[:,list(active).index(1)]),12); comps=[]
        for idx,weight,scale in [(g,.9/12,1.15) for g in groups]+[(np.arange(len(x)),.1,3.)]:
            z=x[idx];cov=np.cov(z,rowvar=False)*scale**2;e,U=np.linalg.eigh(cov);cov=(U*np.maximum(e,1e-6))@U.T;L=np.linalg.cholesky(cov);nu=5
            comps.append(dict(weight=weight,center=z.mean(axis=0),cov=cov,L=L,inv=np.linalg.inv(cov),nu=nu,normalizer=jp.gammaln((nu+k)/2)-jp.gammaln(nu/2)-.5*k*np.log(nu*np.pi)-np.log(np.diag(L)).sum()))
        mixtures.append(comps);covta.append(np.linalg.cholesky(np.cov(x[:,:2],rowvar=False)+np.eye(2)*1e-7))
    m.empirical_jump_counts=np.zeros((42,2),int);m.mode_shift_counts=np.zeros((42,2),int)
    def update(i,theta,ll,mean,sd,lc,eps,rng):
        cand=theta.copy();cand[active]=jp.sample(rng,mixtures[i])
        try:newll=m.subject_ll(cand,m.C@lc,eps,m.data[i])
        except (ArithmeticError,OverflowError):newll=-np.inf
        ratio=newll-ll-.5*np.sum(((cand[active]-mean[active])**2-(theta[active]-mean[active])**2)/sd[active]**2)+jp.logq(theta[active],mixtures[i])-jp.logq(cand[active],mixtures[i])
        m.empirical_jump_counts[i,1]+=1
        if np.log(rng.random())<ratio:m.empirical_jump_counts[i,0]+=1;theta,ll=cand,newll
        cand=theta.copy();cand[:2]+=rng.choice([.15,.5,1.5,4.])*covta[i]@rng.normal(size=2)
        try:newll=m.subject_ll(cand,m.C@lc,eps,m.data[i])
        except (ArithmeticError,OverflowError):newll=-np.inf
        ratio=newll-ll-.5*np.sum(((cand[active]-mean[active])**2-(theta[active]-mean[active])**2)/sd[active]**2)
        m.mode_shift_counts[i,1]+=1
        if np.log(rng.random())<ratio:m.mode_shift_counts[i,0]+=1;theta,ll=cand,newll
        return theta,ll
    m.empirical_joint_update=update
    src=''.join(m.nb['cells'][15]['source']);fn=ast.get_source_segment(src,next(n for n in ast.parse(src).body if isinstance(n,ast.FunctionDef) and n.name=='run_laplace_chain'))
    anchor='if np.log(rng.random())<ratio:theta[i]=cand;lls[i]=newll;accepted[i]+=1';assert fn.count(anchor)==1
    fn=fn.replace(anchor,anchor+'\n            theta[i],lls[i]=empirical_joint_update(i,theta[i],lls[i],m,sig,lc,eps,rng)')
    exec(compile(fn,'stage3_fixed_kernel','exec'),m.__dict__)
    checks=[jp.verify(mixtures[i]) for i in [0,m.subjects.index(7338),41]]
    m.write_json(folder/'proposal_validation.json',dict(checks=checks,random_walk='Fixed covariance and equiprobable scales, zero-centered Gaussian: q(delta)=q(-delta); full original posterior MH ratio.'))
    np.savez_compressed(folder/'proposal_components.npz',centers=np.array([[c['center'] for c in cs] for cs in mixtures]),covariances=np.array([[c['cov'] for c in cs] for cs in mixtures]),weights=[.9/12]*12+[.1],ta_cholesky=covta)

def run(model):
    folder=BASE/model;folder.mkdir(parents=True,exist_ok=True);status=folder/'status.json';started=time.time()
    if status.exists():
        prior=json.loads(status.read_text())
        if prior['status']!='failed' or list(folder.glob('*chain*.npz')):raise RuntimeError('Output already exists; do not duplicate sampling')
        status.rename(folder/('preparation_failure_'+str(int(time.time()))+'.json'))
    m.write_json(status,dict(status='diagnostic_preparation',pid=os.getpid(),started=started))
    try:
        warm=np.load(source(model))['draws'];frame,audit=m.prepare_real_matching();m.data=m.original_prepare(frame);m.write_json(folder/'data_audit.json',audit)
        inspect(model,warm,folder,frame);install(warm,model,folder)
        active=m.allowed(model);k=len(active);center=warm.mean(axis=(0,1));starts=np.zeros((42,11));starts[:,active]=center[2*k+4:].reshape(42,k)
        m.lc_reference=center[2*k:2*k+2];m.eps_reference=center[2*k+2]
        m.proposals=m.build_proposals(m.data,active,starts,m.lc_reference,m.eps_reference);m.proposal_cross=m.proposal_cross_curvature(m.data,active,m.proposals,m.lc_reference,m.eps_reference)
        diag.MODEL=model;seed=2026260000+int(model[1])*100
        m.write_json(folder/'specification.json',dict(model=model,seed=seed,source=str(source(model)),source_sha256=m.hashlib.sha256(source(model).read_bytes()).hexdigest(),hashes=m.hashes,model_changed=False,sampler='Original kernel plus frozen 12-quantile Student-t mixture and symmetric t/a multiscale update',short_iterations=4000,short_burn=1000,formal_iterations=10000,formal_burn=2000,gate='All parameters Rhat<1.01 and bulk/tail ESS>=400',initialization='Last state of each independent historical chain; new seeds, full burn-in'))
        m.write_json(status,dict(status='short_diagnostic_sampling',pid=os.getpid(),started=started))
        dr=diag.batch('diagnostic',seed,4000,1000,list(warm[:,-1]),folder);summary=diag.summarize(dr,folder,'diagnostic');np.savez_compressed(folder/'diagnostic_combined.npz',draws=dr,labels=m.recovery_labels(model),model=model,diagnostics_pass=summary['all_parameters_pass'])
        if not summary['all_parameters_pass']:
            m.write_json(status,dict(status='short_diagnostics_failed',diagnostics=summary,elapsed=time.time()-started));return
        m.write_json(status,dict(status='formal_sampling',pid=os.getpid(),started=started,short_diagnostics=summary))
        formal=diag.batch('formal',seed+20,10000,2000,list(dr[:,-1]),folder);summary=diag.summarize(formal,folder,'formal');np.savez_compressed(folder/'formal_fit.npz',draws=formal,labels=m.recovery_labels(model),model=model,diagnostics_pass=summary['all_parameters_pass'])
        m.write_json(status,dict(status='complete' if summary['all_parameters_pass'] else 'formal_diagnostics_failed',diagnostics=summary,elapsed=time.time()-started))
    except Exception:
        m.write_json(status,dict(status='failed',error=traceback.format_exc(),elapsed=time.time()-started));raise
if __name__=='__main__':run(sys.argv[1])

