# 道德自我关联与知觉决策：HDDM 再分析

[English](README.md) · [完整结果与计算过程](results.ipynb)

## 研究目的

使用 Hu 等（2020）的验证性研究数据，考察 Good-self（好自我）关联是否比 Bad-self（坏自我）和 Good-other（好他人）表现出更快的证据积累，并区分刺激相关证据与反应偏置。

## 数据与实验任务

数据来自[原研究公开材料](https://github.com/hcp4715/moralSelf_ddm)。身份（Self、Other）与道德属性（Good、Bad）构成四种关联条件，三个任务分别拟合。

|任务|判断内容|被试数|有效试次|
|---|---|---:|---:|
|Matching|形状与标签是否匹配|42|24527|
|Identity|刺激属于 Self 还是 Other|41|14493|
|Valence|刺激属于 Good 还是 Bad|41|14431|

输入文件位于 `Data/`。主分析保留正反应时，单位为秒；Matching 排除一个非法按键试次。

## 模型设定

### 刺激编码

模型使用 `HDDMStimCoding`，设置 `split_param="v"`。`response` 编码实际选择，`stim` 编码刺激类别，不以正确或错误重新标记反应边界。

|任务|上边界：response=1；stim=1 对应类别|下边界：response=0；stim=0 对应类别|
|---|---|---|
|Matching|Match（匹配）|Mismatch（不匹配）|
|Identity|Self（自我）|Other（他人）|
|Valence|Good（好）|Bad（坏）|

以朝向上边界为正方向，模型中的有符号漂移为：

$$
\mu=(2\,\mathrm{stim}-1)v+dc.
$$

因此，两个刺激方向共享同一条件下的 $v$，但符号相反；$dc$ 在两个方向上均以相同符号加入。$v$ 是刺激相关的漂移分量，不等同于某一反应方向的总体倾向。

### 模型参数及心理学含义

|参数|计算含义|在本研究任务中的解释|
|---|---|---|
|$v$|单位时间的刺激相关平均证据变化；其符号随刺激类别翻转|Matching 中反映匹配性辨别的证据积累，Identity 中反映 Self/Other 分类，Valence 中反映 Good/Bad 分类。正向 $v$ 增大表示更强的正确方向证据；条件差异不直接等同于偏好。|
|$a$|上下决策边界之间的距离|三项任务中作出选择所需证据量的模型指标。其他参数不变时，较大的 $a$ 通常对应更谨慎的决策。|
|$t$|证据积累以外的反应时间成分，单位为秒|三项任务中感觉编码、反应执行等过程的合计时间；该模型不能进一步区分这些阶段，对应原研究的非决策时间 $t_0$。|
|$z$|相对起点，范围为 0–1；绝对起点为 $az$，0.5 表示居中|积累开始前的反应倾向。$z>0.5$ 分别偏向 Match、Self、Good。模型内部的 `z_trans` 经逆 logit 转换后解释为 $z$。|
|$dc$|加在有符号漂移上的偏移量，两个刺激方向使用相同符号|积累过程中的方向性偏置。$dc>0$ 分别使漂移偏向 Match、Self、Good；它也可能吸收刺激难度不对称，不能单独解释为因果性偏好。|

这些心理学含义依赖模型假设。$z$ 改变积累的起点，$dc$ 改变积累过程中的平均方向，两者并非同一参数。Matching 的 Match 反应偏置不能直接解释为 Good-self 偏好。

### 候选模型与拟合

所有模型均分层估计 $a$、$t$ 和 $v$。$a$、$t$ 在任务内跨条件共享；$v$ 的条件依赖及五个候选模型如下。表中的“共享”指跨实验条件共享，不表示所有被试使用相同参数值。

|任务|$v$ 的条件依赖|共享的刺激方向|
|---|---|---|
|Matching|身份 × 道德属性|Match / Mismatch|
|Identity|道德属性|Self / Other|
|Valence|身份|Good / Bad|

|模型|$z$|$dc$|
|---|---|---|
|M0|固定为 0.5|固定为 0|
|M1|估计，跨条件共享|固定为 0|
|M2|固定为 0.5|估计，跨条件共享|
|M3|估计，跨条件共享|估计，跨条件共享|
|M4|估计，跨条件共享|与 $v$ 使用相同条件划分|

使用 HDDM 信息性先验（`informative=True`），离群混合概率固定为 0.05。每个模型初筛运行一条链，迭代 10000 次、舍弃前 5000 次，`thin=1`。综合 DIC、轨迹和预测检查，各任务暂选 M4，再追加三条配置相同的链。Identity、Valence 随后分别使用四条全新长链重拟合，每链迭代 60000 次、舍弃前 30000 次，不与原链合并。单链初筛不能证明收敛，M4 也不是已确定的唯一最佳模型。

## 分析结果

### Matching

M4 中 Good-self 的漂移率高于 Bad-self 和 Good-other，自我条件的 Good−Bad 差异也大于他人条件。

|漂移率对比|均值差|95%可信区间|
|---|---:|---|
|Good-self − Bad-self|0.666|[0.387, 0.942]|
|Good-self − Good-other|0.438|[0.160, 0.719]|
|Good-other − Bad-other|−0.009|[−0.278, 0.267]|
|身份 × 道德属性交互¹|0.675|[0.284, 1.059]|

¹ (Good-self−Bad-self)−(Good-other−Bad-other)。区间为 95% 等尾后验可信区间。

$z$ 为 0.568，[0.558, 0.578]，表明起点偏向 Match 反应；Good-self 与 Bad-self 的 $dc$ 差为 0.600，[0.427, 0.774]。这些偏置描述匹配判断中的反应倾向。

![Matching参数估计与条件对比](Output/hddm/matching/figures/matching_m4_parameter_contrasts.png)

**图1. Matching 的参数估计与条件对比。** 点为后验均值，线段为 95% 可信区间；$v$ 表示刺激相关漂移，$dc$ 表示漂移偏置，$z$ 表示相对起点。

### 收敛诊断

Matching M4 的四链诊断通过。Identity 和 Valence 的 M4 在独立长链重拟合后仍未收敛，因此不解释两项任务的参数估计。

|任务|最大 Rhat|最低 bulk ESS|最低 tail ESS|判断|
|---|---:|---:|---:|---|
|Matching|1.0022|917|472|通过|
|Identity（长链）|1.3076|10|17|未通过|
|Valence（长链）|1.5282|7|27|未通过|

采用秩归一化 Rhat<1.01、bulk/tail ESS 均≥400，并结合轨迹检查。Rhat 衡量链间一致性，ESS 衡量有效采样信息量。M0–M3 仅完成单链初筛，不能据此判断正式四链收敛。

![Matching四链轨迹与分布](Output/hddm/matching/figures/matching_m4_traces.png)

**图2. Matching 的四链轨迹与参数分布。** 不同颜色表示不同链，左列为采样轨迹，右列为参数分布。

## 模型检验与局限

**后验预测检验（PPC）。** 模型预测的慢反应时尾部比观测数据更长。实验反应窗口约为 0.9–1.2 秒，而预测未施加相同窗口；逐试次窗口尚未恢复，窗口对差异的贡献仍不确定。

![Matching后验预测检验](Output/hddm/matching/figures/matching_m4_ppc.png)

**图3. Matching 的观测数据与后验预测。** 黑色叉号为观测值；彩色点和区间表示预测均值及 95% 预测统计量区间。各面板显示选择比例与反应时分位数。两种补充离群生成约定用于检查预测敏感性，与拟合时的完整离群似然并不严格一致。

**敏感性分析。** 排除 RT≤0.2 秒后，主要群体估计变化很小，但一个个体参数的 tail ESS 为 342，未达到预设门槛。

**参数恢复。** 四个模拟情景均有收敛不足，尚不能确认恢复表现；每情景仅一个模拟数据集，也不足以评价重复抽样的区间覆盖率。

Matching 结果与 Good-self 加工优势的方向一致，但解释仍受拟合和恢复检验限制。本次 $v$ 在刺激方向间共享、$t$ 跨条件共享，并增加 $dc$、改变离群设定，与 2020 模型的估计量不完全相同，不能直接比较参数数值。Identity、Valence 未收敛，也不足以评价原文的跨任务结论。

## 文件结构与复现

`Data/` 保存三任务输入，`Scripts/hddm/` 包含共享计算代码与各任务入口，`Output/` 保存结果表和图。Matching 有 7 个任务级文件，Identity、Valence 各 5 个。[Notebook](results.ipynb) 只读取公开结果；[输出清单](Scripts/hddm/output_manifest.json) 列出 97 个公开结果文件及生成入口。

拟合使用 dockerHDDM（HDDM 1.0.1RC），固定镜像和采样设置见[任务配置](Scripts/hddm/matching/matching_config.json)；Identity、Valence 在各自目录中有对应配置。在该环境中，从本目录执行：

```sh
python Scripts/hddm/reproduce.py plan --output Output_new
python Scripts/hddm/reproduce.py run --output Output_new --workers 2
python Scripts/hddm/reproduce.py verify --output Output_new
```

`plan` 仅显示计划，`run` 执行完整分析。完整运行需要相邻的 `01_reproduction_2020/Data` 和尚不存在的输出目录。任务按顺序执行，`--workers` 接受 1–4，默认 1；配置中的期望值 5/9 不控制当前并行数。恢复分析固定最多两进程，preflight 顺序执行；内存不足时需减少并行任务。

各任务的 `*_run.py` 支持单链 `fit` 和批量阶段；`*_analyze.py` 的 `compare`、`diagnostics`、`ppc` 可独立调用，选项见 `--help`。检查已有公开文件可运行 `python Scripts/hddm/reproduce.py verify --output Output`，仅需 Python 标准库；它不重新计算统计量或判定收敛。Notebook 需要 pandas 和 IPython。

数据库与逐次后验样本不随仓库分发。后处理已使用保存的样本验证，当前完整拟合流程尚未从零重跑；新拟合应写入独立目录，不与设置或代码指纹不同的链混合。

## 参考资料

- [Hu 等（2020）：原始道德自我研究](https://doi.org/10.1525/collabra.301)
- [Pan 等（2025）：dockerHDDM](https://doi.org/10.1177/25152459241298700)
- [刺激编码方法代码](https://github.com/Chuan-Peng-Lab/ddm_stim_coding_posts)
