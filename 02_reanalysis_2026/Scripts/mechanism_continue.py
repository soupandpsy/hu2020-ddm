"""Resumable continuation of the frozen notebook; never executes historical experiments."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import ast, json, sys, time, traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT)
NB=ROOT/'matching_observation_recovery.ipynb'
nb=json.loads(NB.read_text(encoding='utf-8'))
# Execute setup, then only declarations from historical experimental cells.
exec(compile(''.join(nb['cells'][1]['source']),str(NB)+':setup','exec'))
for cell in [2,6,15,18,21,26,29]:
    tree=ast.parse(''.join(nb['cells'][cell]['source']))
    if cell==2:
        nodes=[]
        for node in tree.body:
            if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='validation' for t in node.targets):break
            nodes.append(node)
        nodes += [n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='plot_pilot_trace']
    else:
        nodes=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.Import,ast.ImportFrom))]
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(NB)+':cell'+str(cell),'exec'))
BASE=ROOT/'Output/hddm/matching/continuation_20261026'
OLD=ROOT/'Output/hddm/matching/mechanism_recovery_20261025'
BASE.mkdir(exist_ok=True,parents=True)
original_prepare=prepare
# Resume completed production segments without overwriting them. The transition kernel is unchanged.
fit_source=ast.get_source_segment(''.join(nb['cells'][26]['source']),next(n for n in ast.parse(''.join(nb['cells'][26]['source'])).body if isinstance(n,ast.FunctionDef) and n.name=='recovery_fit'))
old_block="""production=recovery_chain_batch(model,seeds+10,10000,2000,list(warm[:,-1,:]),log)
    draws=np.stack([x[0] for x in production]);table=recovery_diagnostics(draws,model)
    stage=stem.with_name(stem.name+'_segment1.npz')
    if stage.exists():raise FileExistsError(stage)
    np.savez_compressed(stage,draws=draws,seeds=seeds+10,labels=recovery_labels(model),
                        acceptance=np.stack([x[1] for x in production]))"""
new_block="""stage=stem.with_name(stem.name+'_segment1.npz')
    if stage.exists():
        draws=np.load(stage)['draws']
    else:
        production=recovery_chain_batch(model,seeds+10,10000,2000,list(warm[:,-1,:]),log)
        draws=np.stack([x[0] for x in production])
        np.savez_compressed(stage,draws=draws,seeds=seeds+10,labels=recovery_labels(model),
                            acceptance=np.stack([x[1] for x in production]),global_acceptance=np.stack([x[2] for x in production]),
                            collapsed_acceptance=np.stack([x[3] for x in production]))
    table=recovery_diagnostics(draws,model)"""
assert old_block in fit_source
exec(compile(fit_source.replace(old_block,new_block),str(NB)+':resumable_fit','exec'))
_frozen_fit=recovery_fit
OMITTED_INDEX=None
DEADLINE_WIDTH=.001

def fit_manifest(dataset,model,fit_id):
    return dict(model=model,dataset_sha256=hashlib.sha256(Path(dataset).read_bytes()).hexdigest(),
        seed_base=2026120000+fit_id*100,omitted_index=OMITTED_INDEX,deadline_halfwidth=DEADLINE_WIDTH,
        frozen_notebook_cells_sha256=hashlib.sha256(json.dumps([nb['cells'][i]['source'] for i in [1,2,6,15,18,21,26,29]],sort_keys=True).encode()).hexdigest(),
        chains=4,calibration_iterations=2000,production_iterations=10000,production_burn=2000,extension_iterations=8000)

def recovery_fit(dataset,model,fit_id):
    path=RECOVERY_DIR/(Path(dataset).stem.replace('_data','')+'_'+model+'_manifest.json')
    spec=fit_manifest(dataset,model,fit_id)
    if path.exists():
        assert json.loads(path.read_text())==spec,'Checkpoint manifest differs; use a new output directory'
    else:write_json(path,spec)
    return _frozen_fit(dataset,model,fit_id)

def write_json(path,obj):
    tmp=path.with_suffix('.tmp')
    tmp.write_text(json.dumps(obj,indent=2,default=str),encoding='utf-8');tmp.replace(path)

def frame_from(dataset):
    s=np.load(dataset);f=d.copy().reset_index(drop=True)
    assert np.array_equal(f.Subject,s['trial_subject'])
    f['choice']=s['choices'];f['rt_new']=s['rt'];return f

def audit_point():
    dataset=OLD/'H3_r02_data.npz';f=frame_from(dataset);idx=8158
    s=np.load(OLD/'H3_r02_H1_fit.npz');dr=s['draws'];context=recovery_loo_context(f)
    sub=int(f.loc[idx,'Subject']);si=subjects.index(sub)
    full=original_prepare(f);minus=original_prepare(f.drop(index=idx))
    rows=[]
    for ch in range(4):
        for it in np.linspace(0,dr.shape[1]-1,20,dtype=int):
            theta,lam,eps,eta=real_fit_parameters(dr[ch,it],'H1')
            conditional=recovery_conditional_loglik(context,theta,lam,eps,eta)[idx]
            ratio=subject_ll(theta[si],lam,eps,full[si])-subject_ll(theta[si],lam,eps,minus[si])
            ratio+=np.log(eta) if f.loc[idx,'choice']==-2 else np.log1p(-eta)
            rows.append(dict(chain=ch,draw=int(it),conditional=conditional,joint_difference=ratio,error=conditional-ratio))
    table=pd.DataFrame(rows);table.to_csv(BASE/'problem_point_loglik_audit.csv',index=False)
    assert table.error.abs().max()<1e-8
    write_json(BASE/'problem_point.json',dict(index_zero_based=idx,observation=f.loc[idx].to_dict(),
        candidate='H1',generator='H3',replicate=2,max_abs_loglik_error=float(table.error.abs().max()),
        note='No observation discarded. Leave-one-out removes this event only from training and scores it as held out.'))
    print(table.describe(),flush=True)

def exact():
    global RECOVERY_DIR,prepare,OMITTED_INDEX
    OMITTED_INDEX=8158
    RECOVERY_DIR=BASE/'exact_loo_H3_r02_H1_8158';RECOVERY_DIR.mkdir(exist_ok=True)
    # Index retained until prepare; all 42 subjects and remaining 25199 events enter training.
    prepare=lambda frame: original_prepare(frame.drop(index=8158))
    t=time.time();fit,frame=recovery_fit(OLD/'H3_r02_data.npz','H1',200)
    prepare=original_prepare
    saved=np.load(fit);draws=saved['draws'];context=recovery_loo_context(frame)
    # Use every tenth retained draw, separately by chain, for held-out predictive integration.
    ll=[]
    for ch in range(4):
        values=[]
        for row in draws[ch,::10]:
            theta,lam,eps,eta=real_fit_parameters(row,'H1')
            values.append(recovery_conditional_loglik(context,theta,lam,eps,eta)[8158])
        ll.append(values)
    ll=np.asarray(ll);np.savez_compressed(RECOVERY_DIR/'heldout_loglik.npz',loglik=ll)
    score=float(logsumexp(ll)-np.log(ll.size))
    scaled=np.exp(ll-ll.max());mcse=float(az.mcse(az.from_dict(posterior={'density':scaled}),method='mean').density)/scaled.mean()
    old=pd.read_csv(OLD/'H3_r02_H1_loo.csv');old_value=float(old.loc[8158,'elpd_loo'])
    old.loc[8158,'elpd_loo']=score;old['method']='PSIS';old.loc[8158,'method']='exact_refit'
    old.to_csv(RECOVERY_DIR/'corrected_loo.csv',index=False)
    write_json(RECOVERY_DIR/'result.json',dict(status='complete',index=8158,old_psis=old_value,
        exact_elpd=score,change=score-old_value,log_predictive_mcse=mcse,elapsed_seconds=time.time()-t,
        heldout_draws=ll.size,fit=str(fit),training_trials=25199,seed_base=2026140000,
        deadline_proxy='fixed exogenous simulation proxy',diagnostics_pass=True))

def set_width(width):
    global subject_ll,DEADLINE_WIDTH
    DEADLINE_WIDTH=width
    base_ll=subject_ll
    def widened(theta,lam,eps,data,width_override=None,nquad=7,upper=2.,**kwargs):
        return base_ll(theta,lam,eps,data,width=width,nquad=max(21,nquad) if width>.001 else nquad,upper=upper)
    subject_ll=widened

def real_fit_stage(model,width=.001):
    global RECOVERY_DIR
    tag='real_'+model if width==.001 else 'sensitivity_'+model+'_width_'+str(width)
    RECOVERY_DIR=BASE/tag;RECOVERY_DIR.mkdir(exist_ok=True)
    frame,audit=prepare_real_matching()
    dataset=RECOVERY_DIR/'real_data.npz'
    if not dataset.exists():
        np.savez_compressed(dataset,choices=frame.choice.to_numpy(),rt=frame.rt_new.to_numpy(),
            trial_subject=frame.Subject.to_numpy(),source_row=frame.source_row.to_numpy(),
            subject_ids=np.asarray(subjects),data_kind='observed',audit=json.dumps(audit))
    review=dict(supports_exploratory_real_fitting=True,formal_identification_claim=False,
        rationale='All 32 independent-recovery fits pass full diagnostics; pure v versus pure dc discriminate in the pilot. Two replicates per generator do not validate precise recovery rates. Exact LOO repair runs before any final recovery ranking. Real parameter inference and PPC remain conditional on stated structural assumptions.',
        difficulty_reference='sum(log relation multipliers)=0; dc is reference dependent',
        deadline_halfwidth=width,predictive_target='same-subject held-out event conditional on record-derived deadline proxy; not unconditional prospective prediction',
        audit=audit)
    write_json(RECOVERY_DIR/'analysis_specification.json',review)
    if width!=.001:set_width(width)
    fit_id=300+MODELS.index(model) if width==.001 else 400+int(width*1000)
    fit,_=recovery_fit(dataset,model,fit_id)
    if width!=.001:
        # Targeted check for the NEW wider integration range; historical validation is not rerun.
        saved=np.load(fit);draws=saved['draws'];items=original_prepare(frame);checks=[]
        for ch in range(4):
            for it in np.linspace(0,draws.shape[1]-1,3,dtype=int):
                theta,lam,eps,eta=real_fit_parameters(draws[ch,it],model)
                for si in range(42):
                    a=subject_ll(theta[si],lam,eps,items[si],nquad=21)
                    b=subject_ll(theta[si],lam,eps,items[si],nquad=41)
                    checks.append(dict(chain=ch,draw=int(it),subject=subjects[si],loglik21=a,loglik41=b,difference=a-b))
        table=pd.DataFrame(checks);table.to_csv(RECOVERY_DIR/'new_width_quadrature_check.csv',index=False)
        assert table.difference.abs().max()<1e-4,'Wider deadline quadrature insufficient; do not interpret this sensitivity fit'

def real_score_stage(model):
    out=BASE/('real_'+model);fit=out/('real_'+model+'_fit.npz')
    frame,_=prepare_real_matching()
    table=recovery_score(fit,frame)
    bad=table.index[table.pareto_k>table.pareto_threshold].tolist()
    write_json(out/'predictive_status.json',dict(status='needs_exact_refits' if bad else 'conditional_psis_acceptable',
        max_pareto_k=float(table.pareto_k.max()),bad_indices=bad,
        caveat='Conditional on response-derived deadline proxy. Scores are not leakage-free prospective CV. No unequivocal mechanism claim from a marginal maximum.'))

def real_reloo(model,idx):
    global RECOVERY_DIR,prepare,OMITTED_INDEX
    OMITTED_INDEX=idx
    base=BASE/('real_'+model);RECOVERY_DIR=base/('reloo_'+str(idx));RECOVERY_DIR.mkdir(exist_ok=True)
    prepare=lambda frame: original_prepare(frame.drop(index=idx))
    fit,frame=recovery_fit(base/'real_data.npz',model,500+MODELS.index(model)*30000+idx)
    prepare=original_prepare;saved=np.load(fit);dr=saved['draws'];ctx=recovery_loo_context(frame);ll=[]
    for ch in range(4):
        vals=[]
        for row in dr[ch,::10]:
            th,lam,eps,eta=real_fit_parameters(row,model)
            vals.append(recovery_conditional_loglik(ctx,th,lam,eps,eta)[idx])
        ll.append(vals)
    ll=np.asarray(ll);np.savez_compressed(RECOVERY_DIR/'heldout_loglik.npz',loglik=ll)
    v=np.exp(ll-ll.max())
    write_json(RECOVERY_DIR/'result.json',dict(index=idx,exact_elpd=float(logsumexp(ll)-np.log(ll.size)),
        log_predictive_mcse=float(az.mcse(az.from_dict(posterior={'density':v}),method='mean').density)/v.mean(),
        diagnostics_pass=True,conditional_on_record_derived_proxy=True))

def main():
    stage=sys.argv[1] if len(sys.argv)>1 else 'audit'
    state=BASE/(stage+'_status.json');start=time.time()
    write_json(state,dict(stage=stage,status='running',pid=os.getpid(),started=time.time(),notebook_sha256=hashlib.sha256(NB.read_bytes()).hexdigest()))
    try:
        if stage.startswith('fit_'):real_fit_stage(stage.split('_')[1])
        elif stage.startswith('score_'):real_score_stage(stage.split('_')[1])
        elif stage.startswith('sensitivity_'):real_fit_stage('H3',float(stage.split('_')[1]))
        elif stage.startswith('reloo_'):
            _,model,idx=stage.split('_');real_reloo(model,int(idx))
        else:{'audit':audit_point,'exact':exact}[stage]()
        write_json(state,dict(stage=stage,status='complete',elapsed_seconds=time.time()-start,exit_code=0))
    except Exception:
        write_json(state,dict(stage=stage,status='failed',elapsed_seconds=time.time()-start,exit_code=1,error=traceback.format_exc()))
        raise
if __name__=='__main__':main()
