"""Independent read-only validation of the new five-fold conditional scoring."""
import json
import numpy as np,pandas as pd
import mechanism_stage4_kfold as k
m=k.m;out=k.BASE/'kfold_independent_review';out.mkdir(exist_ok=True)
frame,_=m.prepare_real_matching();rng=np.random.default_rng(2026278000);fold=np.full(len(frame),-1)
for _,subject_frame in frame.groupby('Subject',sort=True):
    offset=0
    for _,g in subject_frame.groupby(['condition','stim','relation'],sort=True):
        ii=rng.permutation(g.index.to_numpy());fold[ii]=(offset+np.arange(len(ii)))%5;offset+=len(ii)
counts=pd.crosstab(frame.Subject,fold);assert (counts.to_numpy()==120).all();counts.to_csv(out/'participant_fold_counts.csv')
checks=[];coverage=[]
for fid in range(5):
    train=frame.loc[fold!=fid].copy();test=frame.loc[fold==fid].copy();test['test_position']=np.arange(len(test));items=m.original_prepare(train)
    assert set(train.index).isdisjoint(test.index) and len(train)==20160 and len(test)==5040
    chosen=list(test.groupby('choice',sort=True).head(1).index);extra=[]
    for si,item in enumerate(items):
        assert item['sub']==m.subjects[si] and item['n']==480
        ss=test.loc[test.Subject.eq(item['sub']) & test.choice.ge(0)]
        indices=ss.index[ss.rt_new>max(item['proxy']-.001,item['maxrt'])].tolist();extra+=indices
    chosen=list(dict.fromkeys(chosen+extra));coverage.append(dict(fold=fid,train=len(train),test=len(test),legal_deadline_lower_bound_updates=len(extra),validated_events=len(chosen)))
    for model in ['H0','H1','H2','H3']:
        a=np.load(k.FITS[model]);warm=a['draws']
        for ch,it in [(0,0),(2,-1)]:
            th,lam,eps,eta=m.real_fit_parameters(warm[ch,it],model);conditional=k.heldout_loglik(th,lam,eps,eta,items,test)
            for idx in chosen:
                ev=test.loc[idx];si=m.subjects.index(int(ev.Subject));base=items[si]
                small=pd.concat([train.loc[train.Subject.eq(ev.Subject)],frame.loc[[idx]]]).sort_index();plus=m.original_prepare(small)[0]
                direct=m.subject_ll(th[si],lam,eps,plus)-m.subject_ll(th[si],lam,eps,base)+(np.log(eta) if ev.choice==-2 else np.log1p(-eta))
                val=conditional[int(ev.test_position)]
                checks.append(dict(fold=fid,model=model,chain=ch,draw=it,observation=idx,choice=ev.choice,deadline_update=idx in extra,conditional=val,joint_ratio=direct,error=val-direct))
    print('validated fold',fid,flush=True)
pd.DataFrame(checks).to_csv(out/'joint_ratio_checks.csv',index=False);maxerr=max(abs(x['error']) for x in checks);assert maxerr<1e-8
m.write_json(out/'summary.json',dict(status='passed',max_abs_error=maxerr,checks=len(checks),coverage=coverage,training_only_likelihood=True,per_participant_training_events=480,per_participant_test_events=120,conditional_fixed_full_record_deadline_caveat=True,review='Event-marginal p(event|training) ratios agree with direct shared-deadline joint likelihood differences. This validates scoring algebra at selected full-posterior states; it does not validate future fold convergence or eliminate record-derived-proxy leakage.'))

