"""Matching plot: task-scoped executable entry."""
import os,sys
from pathlib import Path
os.environ['HU2026_TASK']='matching'
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common import entry


def draw_schematic(output):
    """Draw a conceptual diagram; no observations, posterior draws or sampling."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'mathtext.fontset': 'dejavusans'})
    blue, orange, grey, ink = '#236A95', '#B85C22', '#8A939B', '#26313A'
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(14, 7.4))
    fig.subplots_adjust(left=.07, right=.98, top=.84, bottom=.24, wspace=.27)
    fig.suptitle('Stimulus-coded drift diffusion model for the Matching task',
                 fontsize=18, x=.07, ha='left', y=.97)
    fig.text(.07, .91, 'Conceptual illustration only — no observed trajectories or fitted parameter values',
             fontsize=11, color='#58636D')

    ax.set_title('A   Decision process', loc='left', fontweight='bold', pad=20)
    ax.set_xlim(-.17, 1.12); ax.set_ylim(-.12, 1.2)
    ax.hlines([0, 1], 0, 1.05, color=ink, lw=1.7)
    ax.text(1.04, 1.035, 'Match  (response = 1)', ha='right', color=ink)
    ax.text(1.04, -.065, 'Mismatch  (response = 0)', ha='right', color=ink)
    ax.hlines(.5, 0, 1.05, colors=grey, lw=1, linestyles=':')
    ax.text(.99, .515, 'center: z = 0.5', ha='right', fontsize=10, color=grey)
    start = .58
    # Hand-authored conceptual paths, each ending at its first boundary contact.
    xs = [0,.04,.08,.12,.16,.20,.24,.28,.32,.36,.40,.44,.48,.52,.56,.60,.64,.68,.72,.76,.80,.84,.88]
    upper = [.58,.60,.56,.64,.62,.68,.65,.70,.67,.73,.75,.71,.78,.76,.81,.79,.85,.83,.89,.87,.92,.95,1.0]
    lower = [.58,.55,.59,.51,.53,.46,.49,.43,.45,.37,.40,.33,.36,.29,.32,.24,.27,.19,.22,.13,.16,.07,0.0]
    ax.plot(xs, upper, color=blue, lw=1.8)
    ax.plot(xs, lower, color=grey, lw=1.6, ls='--')
    ax.scatter([0,.88,.88], [start,1,0], color=[ink,blue,grey], s=32, zorder=5)
    ax.annotate('start = az', xy=(0,start), xytext=(.02,.84),
                arrowprops={'arrowstyle':'->','color':ink}, color=ink)
    for x, y0, y1, label in [(-.12,0,1,'a'),(-.045,0,start,'az')]:
        ax.annotate('', xy=(x,y1), xytext=(x,y0),
                    arrowprops={'arrowstyle':'<->','color':ink})
        ax.text(x-.018,(y0+y1)/2,label,ha='right',va='center',color=ink)
    ax.annotate('', xy=(.55,.90), xytext=(.22,.69),
                arrowprops={'arrowstyle':'->','color':blue,'lw':2})
    ax.text(.43,.96,r'mean slope $\mu=+v+dc$',color=blue,ha='center',fontsize=11)
    ax.set_xlabel('Time since accumulation begins (schematic)', labelpad=12)
    ax.set_ylabel('Accumulated evidence', labelpad=17)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():spine.set_visible(False)

    bx.set_title('B   Stimulus coding and drift bias', loc='left', fontweight='bold', pad=20)
    bx.set_xlim(-.75,2.15); bx.set_ylim(-1.45,1.85)
    bx.axhline(0,color=grey,lw=1)
    bx.set_yticks([0]); bx.set_yticklabels(['0'])
    bx.set_xticks([]); bx.set_ylabel('Signed drift (toward Match is positive)')
    for spine in ['top','right','bottom']:bx.spines[spine].set_visible(False)
    bx.spines['left'].set_color(grey)
    for x, sign, name in [(0,1,'Match stimulus'),(1.4,-1,'Mismatch stimulus')]:
        base=sign*1.0; shifted=base+.35
        bx.annotate('',xy=(x,base),xytext=(x,0),
                    arrowprops={'arrowstyle':'->','color':grey,'lw':2.2})
        bx.text(x-.07,base/2, r'$+v$' if sign==1 else r'$-v$',ha='right',va='center',color='#606B74')
        bx.annotate('',xy=(x,shifted),xytext=(x,base),
                    arrowprops={'arrowstyle':'->','color':orange,'lw':2.8})
        bx.text(x-.08,base+.175,r'$+dc$',ha='right',va='center',color=orange)
        bx.annotate('',xy=(x+.27,shifted),xytext=(x+.27,0),
                    arrowprops={'arrowstyle':'->','color':blue,'lw':2.5})
        bx.text(x+.27,shifted+(.12 if sign==1 else -.18),
                r'$+v+dc$' if sign==1 else r'$-v+dc$',ha='center',color=blue)
        bx.text(x+.10,-1.35,name,ha='center',fontweight='bold',fontsize=10.5)
    bx.text(.65,1.65,r'Example: $v>0$, $0<dc<v$',ha='center',fontsize=11)
    # A separate RT decomposition avoids locating all nondecision time before accumulation.
    fig.text(.07,.135,r'$RT=t+T_{decision}$',fontsize=15,color=ink)
    fig.patches.extend([
        Rectangle((.07,.075),.085,.027,transform=fig.transFigure,facecolor='#D6DADF',edgecolor='none'),
        Rectangle((.155,.075),.285,.027,transform=fig.transFigure,facecolor=blue,edgecolor='none')])
    fig.text(.1125,.056,'t',ha='center',fontsize=11)
    fig.text(.2975,.056,'decision time',ha='center',fontsize=11)
    fig.text(.07,.015,'t combines encoding and response processes; this strip is a sum, not a stage timeline.',fontsize=9,color='#58636D')
    fig.text(.56,.135,'Gray: stimulus-related drift     Orange: same +dc shift',fontsize=10.5,color=ink)
    fig.text(.56,.10,'Blue: effective drift after adding dc',fontsize=10.5,color=blue)
    fig.text(.56,.055,'z changes the starting position; dc changes the drift.',fontsize=11,color=ink)
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    path=output/'matching_ddm_schematic.png'
    fig.savefig(path,dpi=180,facecolor='white'); plt.close(fig)
    return path

if __name__=='__main__':entry('matching','plot')
