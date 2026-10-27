"""Persist final receipts if the interactive chat is interrupted; no new sampling."""
from pathlib import Path
import os,sys,time,json,hashlib,subprocess
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT)
BASE=ROOT/'Output/hddm/matching/diagnostic_H2H3_20261026_1425'
names=['H2_repair_v2_precision','H3_repair_v2'];terminal={'complete','formal_diagnostics_failed','needs_sampler_investigation','failed'}
for iteration in range(180):
    states={name:json.loads((BASE/name/'status.json').read_text()) if (BASE/name/'status.json').exists() else {} for name in names}
    if all(s.get('status') in terminal for s in states.values()):break
    time.sleep(30)
results=[]
for name in names:
    folder=BASE/name;fit=folder/'formal_fit.npz'
    if not fit.exists():continue
    saved=np.load(fit);diag=pd.read_csv(folder/'formal_diagnostics.csv')
    passed=bool(np.isfinite(diag[['r_hat','ess_bulk','ess_tail']]).all().all() and (diag.r_hat<1.01).all() and (diag.ess_bulk>=400).all() and (diag.ess_tail>=400).all())
    assert passed==bool(saved['diagnostics_pass'])
    finite=bool(np.isfinite(saved['draws']).all());assert finite
    metadata=dict(model=str(saved['model']),status='valid_fit' if passed else 'diagnostics_failed',chains=4,retained_per_chain=saved['draws'].shape[1],
        all_samples_finite=finite,diagnostics_pass=passed,max_rhat=float(diag.r_hat.max()),min_bulk=float(diag.ess_bulk.min()),min_tail=float(diag.ess_tail.min()),
        seeds=saved['seeds'].tolist(),files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [fit,folder/'real_data.npz',folder/'formal_diagnostics.csv',folder/'proposal_components.npz']},
        interpretation='H2/H3 real-data fits only. No reliable four-model ranking while H0/H1 remain unconverged. dc depends on zero-sum log relation reference.')
    (folder/'formal_fit_metadata.json').write_text(json.dumps(metadata,indent=2));results.append(metadata)
    if passed:
        labels=saved['labels'].tolist();draws=saved['draws'];rows=[]
        for label in ['mean_beta_v','mean_beta_dc','mean_log_a','sd_log_a','mean_log_t_excess','sd_log_t_excess']:
            v=draws[:,:,labels.index(label)].ravel() if label in labels else np.zeros(1)
            lo,med,hi=np.quantile(v,[.025,.5,.975]);rows.append(dict(parameter=label,mean=v.mean(),median=med,lower=lo,upper=hi,fixed_zero=label not in labels))
        pd.DataFrame(rows).to_csv(folder/'formal_key_parameters.csv',index=False)
(BASE/'formal_results.json').write_text(json.dumps(results,indent=2))
subprocess.run([sys.executable,'Scripts/mechanism_h2h3_report.py'],check=True)
print(json.dumps(results,indent=2),flush=True)
