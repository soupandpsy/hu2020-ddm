from pathlib import Path
import json,itertools
import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal
BASE=Path(__file__).resolve().parents[1]/'Output/hddm/matching/diagnostic_H2H3_20261026_1425'
for model in ['H2','H3']:
    folder=BASE/(model+'_repair_v2');s=np.load(folder/'proposal_components.npz');rng=np.random.default_rng(2026209911)
    errors=[]
    for i in range(42):
        centers=s['centers'][i];cov=.15**2*s['covariances'][i,3]
        distributions=[]
        for a,b in itertools.permutations(range(3),2):
            for scale,weight in [(.5,.15),(1.,.7),(1.5,.15)]:
                distributions.append((np.log(weight/6),multivariate_normal(mean=scale*(centers[b]-centers[a]),cov=cov)))
        for _ in range(5):
            delta=rng.multivariate_normal(np.zeros(len(centers[0])),s['covariances'][i,3])
            pos=logsumexp([w+d.logpdf(delta) for w,d in distributions]);neg=logsumexp([w+d.logpdf(-delta) for w,d in distributions])
            errors.append(abs(pos-neg))
    assert max(errors)<1e-8
    result=dict(max_logq_delta_minus_logq_negative_delta=max(errors),checked_subjects=42,checks_per_subject=5,
        weights_sum=1.,note='Ordered peak pairs are equiprobable; scale law and centered Gaussian noise independent of state. q(delta)=q(-delta), so exact MH uses original target difference only.')
    (folder/'symmetric_shift_validation.json').write_text(json.dumps(result,indent=2));print(model,result)
