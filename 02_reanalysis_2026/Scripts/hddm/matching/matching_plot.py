"""Matching figures and Wiener paths simulated at fitted group parameter means."""
import os,sys
from pathlib import Path
os.environ['HU2026_TASK']='matching'
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common import entry

SEED=2026102301
DT=0.0001
N_PATHS=2

def simulate_paths(a,z,mu,rng,n=N_PATHS,dt=DT,max_time=10.):
    """dX=mu*dt+dW; stop at first grid hit, interpolate the final segment.

    No endpoint conditioning, outliers, posterior uncertainty, or trial windows.
    """
    import numpy as np
    if not (a>0 and 0<z<1 and dt>0 and n>0 and np.isfinite([a,z,mu,dt]).all()):
        raise ValueError('Invalid diffusion parameters')
    paths=[]
    for _ in range(n):
        x=a*z; times=[0.]; evidence=[x]
        for step in range(1,int(max_time/dt)+1):
            nxt=x+mu*dt+np.sqrt(dt)*rng.normal()
            if nxt<=0 or nxt>=a:
                boundary=0. if nxt<=0 else a
                times.append((step-1+(boundary-x)/(nxt-x))*dt); evidence.append(boundary)
                paths.append((np.array(times),np.array(evidence))); break
            times.append(step*dt); evidence.append(nxt); x=nxt
        else:raise RuntimeError('No boundary hit; path was not discarded or replaced')
    return paths

def draw_schematic(output_root,figures_output=None):
    """Rebuild from public summaries, without HDDM or posterior databases."""
    import csv,hashlib,json
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    root=Path(output_root); source=root/'hddm/matching/m4/final/diagnostics'
    def records(name):
        with (source/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
    precise={r['']:r for r in records('matching_precise_diagnostics.csv')}
    estimates={r['parameter']:r for r in records('matching_parameter_contrasts.csv')}
    effective=records('matching_effective_drift.csv')
    conditions=['Self.Good','Self.Bad','Other.Good','Other.Bad']
    labels=['Good-self','Bad-self','Good-other','Bad-other']
    a=float(precise['a']['mean']); t=float(precise['t']['mean']); z=float(estimates['z']['mean'])
    v=float(estimates['v(Self.Good)']['mean']); dc=float(estimates['dc(Self.Good)']['mean'])
    if len(effective)!=8 or len({(r['condition'],r['stim']) for r in effective})!=8:raise ValueError('Expected eight effective drift rows')
    for r in effective:
        expected=(2*int(r['stim'])-1)*float(estimates['v('+r['condition']+')']['mean'])+float(estimates['dc('+r['condition']+')']['mean'])
        if not np.isclose(expected,float(r['mean']),rtol=0,atol=1e-10):raise ValueError('Inconsistent posterior summaries')
    rng=np.random.RandomState(SEED)
    groups={s:simulate_paths(a,z,(2*s-1)*v+dc,rng) for s in [1,0]}
    xmax=max(ts[-1] for paths in groups.values() for ts,_ in paths)*1.09
    blue,orange,ink,grey='#236A95','#B85C22','#26313A','#8A939B'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'mathtext.fontset':'dejavusans'})
    fig,(ax,bx)=plt.subplots(1,2,figsize=(14,8),gridspec_kw={'width_ratios':[1.08,1]})
    fig.subplots_adjust(left=.12,right=.97,top=.77,bottom=.34,wspace=.32)
    fig.suptitle('Matching M4: simulated decision paths and effective drift',x=.075,ha='left',y=.97,fontsize=18)
    fig.text(.075,.915,'Fitted group parameters; simulated paths are not measured evidence trajectories',color='#58636D',fontsize=11)
    ax.set_title('A   Good-self: paths at posterior means',loc='left',fontweight='bold',pad=22)
    ax.axhline(a,color=ink,lw=1.5); ax.axhline(0,color=ink,lw=1.5); ax.axhline(a*.5,color=grey,lw=.8,ls=':')
    ax.scatter([0],[a*z],color=ink,zorder=5,s=30)
    ax.text(xmax,a+.045,'Match (response = 1)',ha='right',fontsize=10)
    ax.text(xmax,-.085,'Mismatch (response = 0)',ha='right',fontsize=10)
    endpoints=[]
    for stim,color,style in [(1,blue,'-'),(0,orange,'--')]:
        for index,(times,evidence) in enumerate(groups[stim]):
            ax.plot(times,evidence,color=color,ls=style,lw=1,alpha=.8)
            ax.scatter(times[-1],evidence[-1],color=color,s=25,zorder=5)
            endpoints.append(dict(stim=stim,path=index+1,response=int(evidence[-1]>0),decision_time=float(times[-1]),rt=float(t+times[-1])))
    ax.set_xlim(-xmax*.07,xmax); ax.set_ylim(-.14,a+.13)
    ax.set_xlabel('Time since accumulation begins (s)',labelpad=12); ax.set_ylabel('Accumulated evidence')
    ax.set_yticks([0,a*z,a]); ax.set_yticklabels(['0',f'az = {a*z:.3f}',f'a = {a:.3f}'])
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.legend(handles=[Line2D([0],[0],color=blue,label='Match stimulus'),Line2D([0],[0],color=orange,ls='--',label='Mismatch stimulus')],loc='upper left',bbox_to_anchor=(0,-.19),frameon=False,ncol=2,fontsize=10)
    bx.set_title('B   Effective drift: posterior estimates',loc='left',fontweight='bold',pad=22)
    bx.axvline(0,color=grey,lw=1)
    for i,c in enumerate(conditions):
        for stim,offset,color,marker in [(1,-.13,blue,'o'),(0,.13,orange,'s')]:
            r=next(r for r in effective if r['condition']==c and int(r['stim'])==stim)
            mean,lo,hi=[float(r[k]) for k in ['mean','q025','q975']]
            bx.errorbar(mean,i+offset,xerr=[[mean-lo],[hi-mean]],fmt=marker,color=color,capsize=3,ms=5)
    bx.set_yticks(range(4)); bx.set_yticklabels(labels); bx.invert_yaxis()
    bx.set_xlim(-3,3); bx.set_xlabel('Signed drift (positive toward Match)',labelpad=12)
    bx.spines['top'].set_visible(False); bx.spines['right'].set_visible(False)
    bx.grid(axis='x',color='#E9ECEF',lw=.7); bx.set_axisbelow(True)
    bx.legend(handles=[Line2D([0],[0],marker='o',color=blue,label=r'Match: $+v+dc$'),Line2D([0],[0],marker='s',color=orange,label=r'Mismatch: $-v+dc$')],loc='upper left',bbox_to_anchor=(0,-.19),frameon=False,fontsize=10,ncol=2)
    fig.text(.075,.12,f'a = {a:.3f}    z = {z:.3f}    t = {t:.3f} s    v = {v:.3f}    dc = {dc:.3f}',fontsize=11)
    fig.text(.075,.075,r'$RT=t+T_{decision}$; only the DDM component is simulated.',fontsize=11)
    fig.text(.075,.03,'First two paths per stimulus; fixed seed; no selection by outcome or appearance.',fontsize=10,color='#58636D')
    fig.text(.56,.12,'Points: posterior means; bars: 95% equal-tailed credible intervals.',fontsize=10)
    fig.text(.56,.075,'Intervals use paired v and dc draws from all four chains.',fontsize=10)
    fig.text(.56,.03,'Neither panel establishes predictive adequacy or parameter recovery.',fontsize=10,color='#58636D')
    out=Path(figures_output) if figures_output else root/'hddm/matching/figures'; out.mkdir(parents=True,exist_ok=True)
    image=out/'matching_ddm_schematic.png'; fig.savefig(image,dpi=180,facecolor='white'); plt.close(fig)
    inputs=['matching_precise_diagnostics.csv','matching_parameter_contrasts.csv','matching_effective_drift.csv']
    metadata=dict(task='matching',model='m4',condition_panel_a='Self.Good',parameters=dict(a=a,z=z,t=t,v=v,dc=dc,absolute_start=a*z),
        parameter_policy='Separate group posterior means plugged into DDM; a_mean*z_mean is not E[a*z].',
        method='Euler-Maruyama Wiener increments; grid boundary detection with linear final-segment interpolation',
        seed=SEED,dt_seconds=DT,diffusion_sd=1.,paths_per_stimulus=N_PATHS,max_time_seconds=10.,
        outlier_component=False,trial_response_window=False,posterior_uncertainty_in_paths=False,subject_variability_in_paths=False,
        path_selection='First two paths per stimulus, Match then Mismatch; no filtering or resampling',endpoints=endpoints,
        inputs={(source/f).relative_to(root).as_posix():hashlib.sha256((source/f).read_bytes()).hexdigest() for f in inputs},
        numpy_version=np.__version__,matplotlib_version=matplotlib.__version__)
    (out/'matching_ddm_simulation.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    return image

if __name__=='__main__':entry('matching','plot')
