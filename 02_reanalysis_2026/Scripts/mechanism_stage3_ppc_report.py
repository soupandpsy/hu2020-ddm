"""Write the H2 PPC report/notebook from completed stage-three artifacts."""
import json
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'Output/hddm/matching/continuation_stage3_20261026/H2_ppc'
status=json.loads((OUT/'output_status.json').read_text());assert status['status']=='complete'
df=pd.read_csv(OUT/'ppc_summary.csv')
key=df.loc[(df.group.str.endswith('|pooled')&df.metric.isin(['correct_mean_rt','error_mean_rt','timeout_rate','accuracy_given_response']))|df.metric.isin(['dprime_good_minus_bad','correct_rt_good_minus_bad','accuracy_good_minus_bad','fast_accuracy_good_minus_bad','fast_match_choice_good_minus_bad','fast_correct_rt_good_minus_bad'])].copy()
key.to_csv(OUT/'key_ppc_results.csv',index=False)
def md(t):
    lines=['|Group|Statistic|Observed|Predictive median [95% interval]|Outside interval|','|---|---|---:|---:|---|']
    for r in t.itertuples():lines.append(f'|{r.group.replace("|"," / ")}|{r.metric}|{r.observed:.5f}|{r.median:.5f} [{r.lower:.5f}, {r.upper:.5f}]|{not r.observed_in_95_predictive_interval}|')
    return '\n'.join(lines)
intro='''# H2 已通过拟合的真实数据后验预测检验

复用 `diagnostic_H2H3_20261026_1425/H2_repair_v2_precision/formal_fit.npz`，四链每链8000保留样本；未重拟合或修改任何原模型。数据逐数组核对25200真实事件及来源行号。最大Rhat=1.002674、最低bulk ESS=2351.01、tail ESS=1620.84。

每链独立选50个后验样本，共200个完整实验设计复制。每次复制保留42名被试×600事件，通过观测条件下的共同截止时间后验抽取每个被试一个截止值。沿用WFPT密度、均匀污染、超时及非法反应生成规则。图中区间为95%逐统计量预测区间；不是同时置信带，也不是留出预测或模型排名。200复本的尾部分位数仍有Monte Carlo不确定性。多指标中的越界计数只作诊断描述，不解释为全局显著性检验。

快速反应规则在生成前固定：两种Self效价合并，按Match状态使用观测RT第20百分位阈值；另固定观测五分位边界检验整个反应速度形状。统计量采用被试内Good-self减Bad-self再跨被试平均；缺少某一效价的分箱被试不进入该配对值，逐复本保留实际配对人数。d-prime以Match为信号、Match选择为肯定反应，在合法反应内用(count+0.5)/(n+1)校正，然后跨被试平均；它不把超时视为某个选择，不能解释为包括超时的无条件敏感性。

H2的β_log(v)由模型固定为零。β_dc后验中位数0.50131，95%可信区间[0.36249, 0.64104]，只在当前零和log关系难度参照、先验、污染和截止代理模型下成立。该区间不能证明H2优于H1/H3或推翻原论文。参照变换图使用等预测重参数化及推前先验；H2的v目标效应固定零，因而这里的dc效应在所示变换下保持不变，也不能据此推断一般H3中的dc解释不依赖参照。未完成截止/污染先验重拟合稳健性。

## 主要检验值

'''
images=['accuracy_timeout_ppc','rt_means_quantiles_ppc','rt_distribution_ppc','fast_good_self_ppc','fast_shape_ppc','dprime_ppc','target_forest','difficulty_reference_sensitivity']
report=intro+md(key)+'\n\n## 图形\n\n'+'\n\n'.join(f'![{name}]({name}.png)' for name in images)+'\n\n## 可复查文件\n\n`ppc_summary.csv`保留全部分组/统计量的观测、预测均值/中位数/95%区间、偏差、越界距离及有限复本数。`ppc_misfit_flags.csv`保存非CDF指标越界项；`ppc_all.csv`及`ppc_checkpoints/`保留全部复本统计量。分组细化到条件×Match×错配关系，并额外保存跨关系汇总。`ppc_settings.json`记录随机种子、后验抽样索引和固定快速阈值；`extended_statistic_definitions.json`记录定义与源文件哈希。所有图片同时提供PNG/PDF。\n'
(OUT/'REPORT.md').write_text(report,encoding='utf-8')
cells=[dict(cell_type='markdown',metadata={},source=intro.splitlines(True)),dict(cell_type='code',execution_count=None,metadata={},outputs=[],source=['from pathlib import Path\n','import pandas as pd\n','from IPython.display import display, Image\n',"root = Path.cwd()\n","if not (root/'ppc_summary.csv').exists():\n","    root = root/'Output/hddm/matching/continuation_stage3_20261026/H2_ppc'\n","display(pd.read_csv(root/'key_ppc_results.csv'))\n"])]
for name in images:cells.append(dict(cell_type='code',execution_count=None,metadata={},outputs=[],source=[f"display(Image(filename=str(root/'{name}.png')))\n"]))
nb=dict(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3',language='python',name='python3')),nbformat=4,nbformat_minor=5)
(OUT/'H2_ppc_review.ipynb').write_text(json.dumps(nb,ensure_ascii=False,indent=2),encoding='utf-8')
print(md(key))
