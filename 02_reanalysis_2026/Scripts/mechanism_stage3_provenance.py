import json,inspect,hashlib,platform
from pathlib import Path
import numpy as np
import mechanism_continue as m
base=m.ROOT/'Output/hddm/matching/continuation_stage3_20261026'
f,a=m.prepare_real_matching()
p=m.ROOT/'Output/hddm/matching/diagnostic_H2H3_20261026_1425/H2_repair_v2_precision/real_data.npz'
s=np.load(p)
assert np.array_equal(s['choices'],f.choice)
assert np.array_equal(s['trial_subject'],f.Subject)
assert np.allclose(s['rt'],f.rt_new,equal_nan=True,rtol=0,atol=0)
assert np.array_equal(s['source_row'],f.source_row)
obj=dict(data_audit=a,H2_observations_exactly_match=True,notebook_sha256=hashlib.sha256(m.NB.read_bytes()).hexdigest(),python=platform.python_version(),hddm=m.hddm.__version__,numpy=np.__version__,arviz=m.az.__version__,parameter_names=m.names,models={name:dict(active_indices=list(map(int,m.allowed(name))),active_names=[m.names[j] for j in m.allowed(name)]) for name in m.MODELS},prior_mu=m.prior_mu.tolist(),prior_sd=m.prior_sd.tolist() if hasattr(m,'prior_sd') else 'see frozen notebook',relation_reference='sum log multipliers =0',deadline_halfwidth=.001,model_change=False,retained_subjects=42,retained_events=25200)
m.write_json(base/'comparability_manifest.json',obj)
print(json.dumps(obj,indent=2,default=str))
