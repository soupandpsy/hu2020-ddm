"""Four independent chains/model, short gate BEFORE any long production run."""
import os,sys,time,json,traceback,contextlib
from pathlib import Path
import numpy as np
import pandas as pd
import mechanism_continue as m
from multiprocessing import get_context
BASE=m.ROOT/'Output/hddm/matching/diagnostic_H2H3_20261026_1425'
MODEL=sys.argv[1] if __name__=='__main__' else None
PRECISION='--precision' in sys.argv
REPAIR='--repair' in sys.argv or '--repair-v2' in sys.argv or PRECISION
REPAIR_V2='--repair-v2' in sys.argv or PRECISION

def job(args):
    chain,seed,iterations,burn,state,stage,folder=args
    folder=Path(folder);target=folder/f'{stage}_chain{chain}.npz'
    if target.exists():
        old=np.load(target);assert int(old['seed'])==seed
        return old['draws'],old['subject_acceptance'],old['global_acceptance'],old['collapsed_acceptance']
    status=folder/f'{stage}_chain{chain}_status.json';started=time.time()
    m.write_json(status,dict(status='running',pid=os.getpid(),seed=seed,iterations=iterations,burn=burn,started_epoch=started))
    try:
        with open(folder/f'{stage}_chain{chain}.log','a',buffering=1) as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
            print('seed',seed,'pid',os.getpid(),'stage',stage,flush=True)
            result=m.run_laplace_chain(seed,model=MODEL,iterations=iterations,burn=burn,
                global_correction=True,dispersed=state is None,collapsed=True,initial_state=state,
                progress_path=str(folder/f'{stage}_chain{chain}_progress.log'))
        np.savez_compressed(target,draws=result[0],subject_acceptance=result[1],global_acceptance=result[2],collapsed_acceptance=result[3],seed=seed,model=MODEL,
            empirical_jump_counts=getattr(m,'empirical_jump_counts',np.zeros((42,2))),mode_shift_counts=getattr(m,'mode_shift_counts',np.zeros((42,2))))
        m.write_json(status,dict(status='complete',pid=os.getpid(),seed=seed,exit_code=0,elapsed_seconds=time.time()-started))
        return result
    except Exception:
        m.write_json(status,dict(status='failed',seed=seed,exit_code=1,error=traceback.format_exc()));raise

def batch(stage,base,iterations,burn,states,folder):
    jobs=[(i,base+i,iterations,burn,states[i],stage,str(folder)) for i in range(4)]
    with get_context('fork').Pool(4) as pool:result=pool.map(job,jobs)
    return np.stack([r[0] for r in result])

def summarize(draws,folder,stage):
    table=m.recovery_diagnostics(draws,MODEL);table.to_csv(folder/(stage+'_diagnostics.csv'),index_label='parameter')
    target=table.loc[[p for p in table.index if p.endswith(('log_a','log_t_excess'))]]
    target.to_csv(folder/(stage+'_a_t_diagnostics.csv'),index_label='parameter')
    report=dict(max_rhat=float(table.r_hat.max()),min_bulk_ess=float(table.ess_bulk.min()),min_tail_ess=float(table.ess_tail.min()),
        all_parameters_pass=m.recovery_gate(table),a_t_max_rhat=float(target.r_hat.max()),a_t_min_bulk=float(target.ess_bulk.min()),
        a_t_min_tail=float(target.ess_tail.min()),severe_mixing=bool((target.r_hat>1.05).any() or (target.ess_bulk<100).any() or (target.ess_tail<100).any()),
        retained_draws_per_chain=draws.shape[1],criterion='Formal sampling allowed only if all checked quantities Rhat<1.01 and bulk/tail ESS>=400. Severe t/a flag: Rhat>1.05 or ESS<100. Short-run failure is not proof of model failure.')
    m.write_json(folder/(stage+'_summary.json'),report)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels=m.recovery_labels(MODEL)
    selected=['mean_log_a','sd_log_a','mean_log_t_excess','sd_log_t_excess']+target.sort_values('r_hat',ascending=False).head(6).index.tolist()
    selected=list(dict.fromkeys(selected));fig,axes=plt.subplots(len(selected),2,figsize=(13,2*len(selected)),constrained_layout=True)
    for row,name in enumerate(selected):
        j=labels.index(name)
        for ch in range(4):
            axes[row,0].plot(draws[ch,:,j],lw=.4,alpha=.7,label=f'Chain {ch+1}')
            axes[row,1].hist(draws[ch,:,j],bins=40,density=True,histtype='step',label=f'Chain {ch+1}')
        axes[row,0].set_ylabel(name)
    axes[0,1].legend();fig.suptitle(f'{MODEL}: {stage} sampler diagnostics, not inference')
    for ext in ['png','pdf']:fig.savefig(folder/(stage+'_a_t_chains.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    return report

def run():
    assert MODEL in ['H2','H3']
    folder=BASE/(MODEL+'_repair_v2_precision' if PRECISION else MODEL+'_repair_v2' if REPAIR_V2 else MODEL+'_repair_v1' if REPAIR else MODEL);folder.mkdir(parents=True,exist_ok=True);status=folder/'status.json';start=time.time()
    seed=2026200000+(2 if MODEL=='H2' else 3)*1000+(300 if PRECISION else 200 if REPAIR_V2 else 100 if REPAIR else 0)
    spec=dict(model=MODEL,chains=4,calibration_iterations=0 if REPAIR else 1000,calibration_burn=0 if REPAIR else 500,diagnostic_iterations=2000,diagnostic_burn=0 if PRECISION else 500,
        seed_base=seed,source_hashes=m.hashes,sampler='empirical mixture plus symmetric mode shifts' if REPAIR_V2 else 'frozen empirical multimodal MH supplement' if REPAIR else 'original frozen quadratic proposal',model_and_likelihood='Frozen original H0-H3; zero-sum log relation difficulty; shared deadline +/-1 ms',
        production_rule='Only after complete short-run gate; otherwise investigate mixing without blindly extending',max_concurrent_models=2,blas_threads=1)
    manifest=folder/'specification.json'
    if manifest.exists():
        previous=json.loads(manifest.read_text());previous.setdefault('sampler','original frozen quadratic proposal')
        assert previous==spec
    else:m.write_json(manifest,spec)
    try:
        m.write_json(status,dict(status='calibration',pid=os.getpid(),started_epoch=start))
        frame,audit=m.prepare_real_matching();m.write_json(folder/'data_audit.json',audit)
        np.savez_compressed(folder/'real_data.npz',choices=frame.choice.to_numpy(),rt=frame.rt_new.to_numpy(),trial_subject=frame.Subject.to_numpy(),source_row=frame.source_row.to_numpy(),data_kind='observed')
        m.data=m.original_prepare(frame);active=m.allowed(MODEL);k=len(active)
        m.lc_reference=np.zeros(2);m.eps_reference=.05
        if REPAIR:
            source=BASE/MODEL/'diagnostic_combined.npz'
            warm=np.load(source)['draws']
            if REPAIR_V2:
                newer=np.load(BASE/(MODEL+'_repair_v1')/'diagnostic_combined.npz')['draws']
                warm=np.concatenate([warm,newer],axis=1)
            from mechanism_joint_proposals import install
            validation,mixtures=install(m,warm,MODEL,mode_shifts=REPAIR_V2)
            m.write_json(folder/'proposal_validation.json',dict(**validation,reference_sha256=m.hashlib.sha256(source.read_bytes()).hexdigest(),
                note='Reference draws train fixed proposals ONLY, not priors or reported posteriors; new burn-in and independent seeds. Full original likelihood and hyperprior used in MH correction.'))
            np.savez_compressed(folder/'proposal_components.npz',centers=np.array([[c['center'] for c in cs] for cs in mixtures]),covariances=np.array([[c['cov'] for c in cs] for cs in mixtures]),weights=[.3,.3,.3,.1],df=5)
        else:
            m.proposals=m.build_proposals(m.data,active,np.tile(m.prior_mu,(42,1)),m.lc_reference,m.eps_reference)
            m.proposal_cross=m.proposal_cross_curvature(m.data,active,m.proposals,m.lc_reference,m.eps_reference)
            warm=batch('calibration',seed,1000,500,[None]*4,folder)
        center=warm.mean(axis=(0,1));starts=np.zeros((42,11));starts[:,active]=center[2*k+4:].reshape(42,k)
        m.lc_reference=center[2*k:2*k+2];m.eps_reference=center[2*k+2]
        m.proposals=m.build_proposals(m.data,active,starts,m.lc_reference,m.eps_reference)
        m.proposal_cross=m.proposal_cross_curvature(m.data,active,m.proposals,m.lc_reference,m.eps_reference)
        m.write_json(status,dict(status='short_diagnostic_sampling',pid=os.getpid(),started_epoch=start))
        if PRECISION:
            previous=np.load(BASE/(MODEL+'_repair_v2')/'diagnostic_combined.npz')['draws']
            extra=batch('diagnostic',seed+10,2000,0,list(previous[:,-1]),folder)
            draws=np.concatenate([previous,extra],axis=1)
            m.write_json(folder/'precision_continuation.json',dict(previous=str(BASE/(MODEL+'_repair_v2')/'diagnostic_combined.npz'),
                previous_draws_per_chain=previous.shape[1],additional_draws_per_chain=2000,proposal_adaptation=False,
                rationale='Serious t/a mixing resolved; limited precision check for remaining Rhat 1.01475 and tail ESS 364. Same frozen transition kernel.'))
        else:draws=batch('diagnostic',seed+10,2000,500,list(warm[:,-1]),folder)
        summary=summarize(draws,folder,'diagnostic')
        np.savez_compressed(folder/'diagnostic_combined.npz',draws=draws,labels=m.recovery_labels(MODEL),model=MODEL,seeds=np.arange(4)+seed+10,diagnostics_pass=summary['all_parameters_pass'])
        if not summary['all_parameters_pass']:
            m.write_json(status,dict(status='needs_sampler_investigation',pid=os.getpid(),exit_code=0,elapsed_seconds=time.time()-start,diagnostics=summary,formal_sampling_started=False));return
        m.write_json(status,dict(status='formal_sampling',pid=os.getpid(),started_epoch=start))
        production=batch('formal',seed+20,10000,2000,list(draws[:,-1]),folder)
        final=summarize(production,folder,'formal')
        np.savez_compressed(folder/'formal_fit.npz',draws=production,labels=m.recovery_labels(MODEL),model=MODEL,seeds=np.arange(4)+seed+20,diagnostics_pass=final['all_parameters_pass'])
        m.write_json(status,dict(status='complete' if final['all_parameters_pass'] else 'formal_diagnostics_failed',exit_code=0,elapsed_seconds=time.time()-start,diagnostics=final))
    except Exception:
        m.write_json(status,dict(status='failed',exit_code=1,error=traceback.format_exc(),elapsed_seconds=time.time()-start));raise
if __name__=='__main__':run()
