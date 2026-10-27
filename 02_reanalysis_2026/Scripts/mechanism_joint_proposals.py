"""Frozen multimodal Student-t proposals; exact MH correction preserves posterior."""
import numpy as np
from scipy.special import gammaln,logsumexp
from scipy.stats import multivariate_t

def train(samples,active):
    # Fit proposals only, never alter priors. Distinct t ranges retain separated modes.
    tpos=list(active).index(1);order=np.argsort(samples[:,tpos]);groups=np.array_split(order,3)
    components=[];k=len(active)
    for indices,weight,inflation in [(g,.3,1.2) for g in groups]+[(np.arange(len(samples)),.1,3.)]:
        x=samples[indices];center=np.mean(x,axis=0);cov=np.cov(x,rowvar=False)*inflation**2
        e,U=np.linalg.eigh(cov);cov=(U*np.maximum(e,1e-5))@U.T
        L=np.linalg.cholesky(cov);nu=5.
        normalizer=gammaln((nu+k)/2)-gammaln(nu/2)-.5*k*np.log(nu*np.pi)-np.log(np.diag(L)).sum()
        components.append(dict(weight=weight,center=center,cov=cov,L=L,inv=np.linalg.inv(cov),nu=nu,normalizer=normalizer))
    return components

def logq(x,components):
    values=[]
    for c in components:
        z=x-c['center'];k=len(z)
        values.append(np.log(c['weight'])+c['normalizer']-.5*(c['nu']+k)*np.log1p((z@c['inv']@z)/c['nu']))
    return logsumexp(values)

def sample(rng,components):
    c=components[rng.choice(len(components),p=[x['weight'] for x in components])]
    return c['center']+c['L']@rng.normal(size=len(c['center']))/np.sqrt(rng.chisquare(c['nu'])/c['nu'])

def verify(components):
    rng=np.random.default_rng(2026209900);errors=[];balance=[]
    for _ in range(50):
        x=sample(rng,components);y=sample(rng,components)
        independent=logsumexp([np.log(c['weight'])+multivariate_t.logpdf(x,loc=c['center'],shape=c['cov'],df=c['nu']) for c in components])
        errors.append(abs(logq(x,components)-independent))
        # Arbitrary non-Gaussian target; verifies the two-direction MH balance algebra.
        tx=-np.sum(x*x)/2+.2*np.cos(x).sum();ty=-np.sum(y*y)/2+.2*np.cos(y).sum()
        ratio=ty-tx+logq(x,components)-logq(y,components)
        balance.append(abs(tx+logq(y,components)+min(0,ratio)-(ty+logq(x,components)+min(0,-ratio))))
    assert max(errors)<1e-8 and max(balance)<1e-8
    return dict(max_scipy_logq_error=max(errors),max_detailed_balance_error=max(balance))

def install(m,reference_draws,model,mode_shifts=False):
    import ast
    active=m.allowed(model);k=len(active);flat=reference_draws.reshape(-1,reference_draws.shape[-1])
    mixtures=[train(flat[:,2*k+4+i*k:2*k+4+(i+1)*k],active) for i in range(42)]
    validation=verify(mixtures[0])
    m.empirical_jump_counts=np.zeros((42,2),dtype=int)
    m.mode_shift_counts=np.zeros((42,2),dtype=int)
    def update(i,theta,ll,mean,sd,lc,eps,rng):
        cand=theta.copy();cand[active]=sample(rng,mixtures[i])
        try:newll=m.subject_ll(cand,m.C@lc,eps,m.data[i])
        except (ArithmeticError,OverflowError):newll=-np.inf
        ratio=newll-ll-.5*np.sum(((cand[active]-mean[active])**2-(theta[active]-mean[active])**2)/sd[active]**2)
        ratio+=logq(theta[active],mixtures[i])-logq(cand[active],mixtures[i])
        m.empirical_jump_counts[i,1]+=1
        if np.log(rng.random())<ratio:
            m.empirical_jump_counts[i,0]+=1;theta,ll=cand,newll
        if mode_shifts:
            # Pair, sign, scale and noise are independent of current state.
            # q(delta)=q(-delta), hence MH needs only the exact target ratio.
            a,b=rng.choice(3,size=2,replace=False)
            shift=mixtures[i][b]['center']-mixtures[i][a]['center']
            scale=rng.choice([.5,1.,1.5],p=[.15,.7,.15])
            cand=theta.copy();cand[active]+=scale*shift+.15*mixtures[i][3]['L']@rng.normal(size=k)
            try:newll=m.subject_ll(cand,m.C@lc,eps,m.data[i])
            except (ArithmeticError,OverflowError):newll=-np.inf
            ratio=newll-ll-.5*np.sum(((cand[active]-mean[active])**2-(theta[active]-mean[active])**2)/sd[active]**2)
            m.mode_shift_counts[i,1]+=1
            if np.log(rng.random())<ratio:m.mode_shift_counts[i,0]+=1;theta,ll=cand,newll
        return theta,ll
    m.empirical_joint_update=update
    source=''.join(m.nb['cells'][15]['source']);tree=ast.parse(source)
    fn=ast.get_source_segment(source,next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_laplace_chain'))
    anchor='if np.log(rng.random())<ratio:theta[i]=cand;lls[i]=newll;accepted[i]+=1'
    assert fn.count(anchor)==1
    replacement=anchor+'\n            theta[i],lls[i]=empirical_joint_update(i,theta[i],lls[i],m,sig,lc,eps,rng)'
    exec(compile(fn.replace(anchor,replacement),'frozen_joint_mixture_kernel','exec'),m.__dict__)
    return validation,mixtures
