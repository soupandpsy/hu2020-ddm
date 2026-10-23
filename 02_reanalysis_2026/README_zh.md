# 第二部分：道德自我关联与知觉决策：HDDM 再分析

[English](README.md) · 中文 · [项目首页](../README_zh.md) · [分析结果 Notebook](results.ipynb)

## 研究目的

本分析使用 Hu 等（2020）的 Confirm 数据，采用刺激编码的层级漂移扩散模型（Hierarchical Drift Diffusion Model, HDDM），考察道德自我关联对知觉决策过程的影响。

Matching 任务主要检验 Good-self（好自我）与其他关联条件的证据积累差异；Identity 和 Valence 任务分别考察道德属性与身份信息对分类决策的影响。模型同时区分证据积累率、起点偏置和漂移偏置。

由于模型参数化与原研究不同，本项目属于探索性再分析，而非原模型的直接数值复现。

## 数据与实验任务

数据来自[原研究公开材料](https://github.com/hcp4715/moralSelf_ddm)。身份（Self、Other）与道德属性（Good、Bad）构成四种关联条件，三个任务分别拟合。

| 任务       | 判断内容               | 被试数 | 有效试次   |
| -------- | ------------------ | --- | ------ |
| Matching | 形状与标签是否匹配          | 42  | 24,527 |
| Identity | 刺激关联 Self 还是 Other | 41  | 14,493 |
| Valence  | 刺激关联 Good 还是 Bad   | 41  | 14,431 |

输入文件位于 `Data/`。主分析保留正反应时，单位为秒；Matching 排除了一个非法按键试次。

## 模型设定

### 刺激编码

模型使用 `HDDMStimCoding`，设置 `split_param="v"`。

`stim` 编码刺激对应的正确反应方向，`response` 编码被试实际选择的反应方向。两者相同时为正确反应，不同时为错误反应。与正确/错误反应编码不同，刺激编码保留具体反应类别对应的决策边界。

| 任务       | 上边界（response = 1） | 下边界（response = 0） |
| -------- | ----------------- | ----------------- |
| Matching | Match（匹配）         | Mismatch（不匹配）     |
| Identity | Self（自我）          | Other（他人）         |
| Valence  | Good（好）           | Bad（坏）            |

其中，`stim = 1` 表示正确反应对应上边界，`stim = 0` 表示正确反应对应下边界。

以朝向上边界为正方向，模型的有符号漂移为：

$$
\mu=(2\,\mathrm{stim}-1)v+dc
$$

两个刺激方向共享同一条件下的漂移幅度 $`v`$，但符号相反；漂移偏置 $`dc`$ 在两个方向上以相同符号加入。因此，$`v`$ 表示刺激相关的证据积累成分，而非对某一反应方向的总体倾向。

### 模型参数及心理学含义

| 参数                          | 计算含义                 | 在本研究中的解释                                                                                         |
| --------------------------- | -------------------- | ------------------------------------------------------------------------------------------------ |
| **$`v`$** Drift rate          | 单位时间内刺激相关证据的平均积累率    | Matching 中反映匹配判断的证据积累效率；Identity 中反映 Self/Other 分类；Valence 中反映 Good/Bad 分类。较大的正向 $`v`$ 表示更强的正确方向证据 |
| **$`a`$** Boundary separation | 两个决策边界之间的距离          | 反映完成判断所需的证据量。其他参数不变时，较大的 $`a`$ 通常对应更谨慎的决策                                                          |
| **$`t`$** Non-decision time   | 证据积累以外的反应时间成分，单位为秒   | 包括感觉编码、动作准备和反应执行等过程，模型不能进一步区分这些阶段，对应原研究的非决策时间 $`t_0`$                                              |
| **$`z`$** Starting point      | 证据积累开始时的相对位置，范围为 0–1 | 反映积累开始前的反应偏置；$`z=0.5`$ 表示起点居中                                                                      |
| **$`dc`$** Drift criterion    | 加在刺激相关漂移上的方向性偏移      | 反映证据积累过程中的反应偏置，也可能受到刺激难度不对称等因素影响                                                                 |

**漂移率 $`v`$** 是检验 Good-self 加工优势的主要参数。在 Matching 中，四种关联条件分别具有对应的漂移幅度参数。较大的 $`v`$ 表示更高的模型估计证据积累效率，但不能直接确定这一优势发生于哪个具体的知觉或神经加工阶段。

**起点 $`z`$** 表示证据积累开始前的反应倾向。在当前编码下，$`z>0.5`$ 分别表示起点更接近 Matching 的 Match、Identity 的 Self 或 Valence 的 Good 边界。Matching 中的 Match 反应偏置不能直接解释为 Good-self 偏好。

**漂移偏置 $`dc`$** 表示积累过程中的方向性偏移。正向 $`dc`$ 分别使漂移朝 Match、Self 或 Good 边界偏移，负向 $`dc`$ 则朝相反边界偏移。

$`z`$ 改变证据积累的初始位置，$`dc`$ 改变积累过程中的平均漂移方向。两种偏置能否从现有数据中可靠区分，仍需结合参数恢复结果判断。

相对起点为 $`z`$，绝对起点为 $`az`$；内部采样参数 `z_trans` 经逆 logit 转换后才作为 $`z`$ 解释。

![Matching M4 fitted-parameter simulations and effective drift](Output/hddm/matching/figures/matching_ddm_schematic.png)

**图 1. Matching M4 的模拟证据积累路径与有效漂移。**

A：根据 Good-self 条件的群体后验均值模拟 Match 与 Mismatch 刺激的证据积累路径。B：四种关联条件的有效漂移及其 95% 等尾可信区间，由配对的 v、dc 后验样本计算。模拟路径不代表实测证据轨迹，也不构成完整的后验预测检验。

### 候选模型与拟合

所有模型均分层估计 $`a`$、$`t`$ 和 $`v`$。其中，$`a`$ 和 $`t`$ 在任务内跨实验条件共享，$`v`$ 随相应任务的实验条件变化。跨条件共享不表示所有被试具有相同的参数值。

| 任务       | $`v`$ 的条件依赖          | 共享漂移幅度的刺激方向      |
| -------- | ------------------ | ---------------- |
| Matching | Identity × Valence | Match / Mismatch |
| Identity | Valence            | Self / Other     |
| Valence  | Identity           | Good / Bad       |

三个任务分别比较五个候选模型：

| 模型 | 起点 $`z`$     | 漂移偏置 $`dc`$  |
| -- | ---------- | ---------- |
| M0 | 固定为 0.5    | 固定为 0      |
| M1 | 自由估计，跨条件共享 | 固定为 0      |
| M2 | 固定为 0.5    | 自由估计，跨条件共享 |
| M3 | 自由估计，跨条件共享 | 自由估计，跨条件共享 |
| M4 | 自由估计，跨条件共享 | 随相应实验条件变化  |

模型使用 HDDM 信息性先验（`informative=True`），离群混合概率固定为 0.05。

每个候选模型首先运行一条 MCMC 链，迭代 10,000 次，舍弃前 5,000 次，`thin=1`。结合 DIC、采样轨迹和初步后验预测指标，各任务暂定选择 M4，并追加三条相同配置的链。

Identity 和 Valence 随后分别进行独立长链重拟合，每条链迭代 60,000 次、舍弃前 30,000 次。长链结果不与原有链合并。

由于初筛阶段每个模型仅有一条链，模型选择仍存在不确定性。

## 分析结果

### Matching

候选 M4 中，Good-self 的漂移幅度高于 Bad-self 和 Good-other，自我条件的 Good−Bad 差异也大于他人条件。

| 漂移幅度对比                  | 后验均值差  | 95% 可信区间        |
| ----------------------- | ------ | --------------- |
| Good-self − Bad-self    | 0.666  | [0.387, 0.942]  |
| Good-self − Good-other  | 0.438  | [0.160, 0.719]  |
| Good-other − Bad-other  | −0.009 | [−0.278, 0.267] |
| Identity × Valence 交互对比 | 0.675  | [0.284, 1.059]  |

交互对比定义为：

$$
(v_{GS}-v_{BS})-(v_{GO}-v_{BO})
$$

区间为 95% 等尾后验可信区间。

起点参数 $`z`$ 的后验均值为 0.568，95% 可信区间为 [0.558, 0.578]，对应向 Match 边界的起点偏置。

Good-self 与 Bad-self 的 $`dc`$ 差异为 0.600 [0.427, 0.774]，表明两种条件在朝向 Match 边界的漂移偏置上存在差异。

上述结果均为给定候选 M4 及其参数约束的后验估计。

![Matching M4](Output/hddm/matching/figures/matching_m4_parameter_contrasts.png)

**图 2. Matching M4 的参数估计与条件对比。**

点表示后验均值，线段表示 95% 可信区间。
### 收敛诊断

| 任务              | 最大 R-hat | 最低 bulk ESS | 最低 tail ESS | 结果  |
| --------------- | -------- | ----------- | ----------- | --- |
| Matching M4     | 1.0022   | 917         | 472         | 通过  |
| Identity M4（长链） | 1.3076   | 10          | 17          | 未通过 |
| Valence M4（长链）  | 1.5282   | 7           | 27          | 未通过 |

采用秩归一化 R-hat < 1.01、bulk ESS ≥ 400、tail ESS ≥ 400，并结合采样轨迹检查。

Matching M4 达到预设的四链收敛诊断标准。Identity 和 Valence 在独立长链重拟合后仍存在明显的收敛问题，因此不对其参数估计进行推断性解释。

M0—M3 仅有单链初筛结果，不能据此判断多链收敛。

![Matching M4](Output/hddm/matching/figures/matching_m4_traces.png)

**图 3. Matching M4 的四链采样轨迹与后验分布。**

## 模型检验与局限

### 后验预测检验（PPC）

Matching M4 预测的慢反应时尾部比观测数据更长。原实验的反应窗口约为 0.9–1.2 秒，而当前预测未施加相同的逐试次反应窗口，因此尚无法确定反应窗口对预测差异的影响。

此外，补充使用的两种离群值生成方式与拟合时的完整离群似然并不严格一致。

![Matching M4](Output/hddm/matching/figures/matching_m4_ppc.png)

**图 4. Matching M4 的观测与预测选择比例及反应时分位数。**

黑色叉号表示观测值，彩色点与区间表示后验预测统计量。
### 敏感性分析

排除 RT ≤ 0.2 秒后，Matching 的主要群体参数估计变化较小，但一个个体参数的 tail ESS 为 342，低于预设标准。

### 参数恢复

Matching M4 的四个模拟情景均存在采样收敛不足，尚不能确认参数恢复表现。每个情景仅使用一个模拟数据集，也不足以评估重复抽样下的区间覆盖率。

### 结果解释

Matching 的结果在方向上与 Good-self 证据积累优势一致，但仍受到模型选择、后验预测和参数恢复的限制。

与原研究相比，本次分析在刺激方向间共享 $`v`$，使 $`t`$ 跨条件共享，并引入 $`dc`$ 和不同的离群处理设置。因此，两项分析的模型估计量不完全相同，参数数值不宜直接比较。

Identity 和 Valence 尚未获得满足收敛诊断要求的结果，无法据此评价跨任务效应。

## 文件与复现方法

| 目录或文件                    | 内容                                |
| ------------------------ | --------------------------------- |
| `Data/`                  | Matching、Identity、Valence 的模型输入数据 |
| `Scripts/hddm/common.py` | 三任务共享的 HDDM 模型构建、拟合与分析函数          |
| `Scripts/hddm/matching/` | Matching 的模型拟合、分析、绘图、参数恢复与敏感性分析   |
| `Scripts/hddm/identity/` | Identity 的模型拟合、分析与绘图              |
| `Scripts/hddm/valence/`  | Valence 的模型拟合、分析与绘图               |
| `Output/`                | 模型比较、参数估计、诊断结果与图表                 |
| `results.ipynb`          | 已保存分析结果的汇总与展示                     |

公开结果文件及其来源记录在 [output_manifest.json](Scripts/hddm/output_manifest.json) 中。

模型拟合使用 dockerHDDM（HDDM 1.0.1RC），镜像与采样配置见任务配置文件。

在 `02_reanalysis_2026/` 目录中查看复现计划：

```bash
python Scripts/hddm/reproduce.py plan --output Output_new
```

在配置指定的 HDDM 环境中从头运行：

```bash
python Scripts/hddm/reproduce.py run --output Output_new --workers 2
```

完整运行需要相邻的 `01_reproduction_2020/Data` 目录和尚不存在的新输出目录。各任务顺序执行，`--workers` 支持 1–4 个并行任务，默认 1；内存不足时应减少并行数。

检查已公开结果文件：

```bash
python Scripts/hddm/reproduce.py verify --output Output
```

该命令仅检查文件完整性与基本结构，不重新计算模型参数，也不验证模型收敛。

当前完整拟合流程尚未在整理后的代码上进行端到端重跑。仓库提供模型输入、分析代码和汇总结果，但不包含完整 MCMC 数据库与逐次后验样本。

使用公开汇总结果生成图 1（需要 NumPy 和 Matplotlib，不需 HDDM 或本地后验数据库）：

```bash
python Scripts/hddm/matching/matching_plot.py --kind schematic --output Output --figures-output Figures_check
```

图像及模拟元数据写入 `Figures_check/`，避免覆盖已公开结果。

有效漂移表由 `python Scripts/hddm/matching/matching_analyze.py diagnostics --effective-drift --output Output` 生成，需要四条正式链的后验样本与 NumPy、pandas、ArviZ；绘图只读取公开汇总表。随机种子、数值步长、参数来源哈希和各路径的选择及反应时保存在 `Output/hddm/matching/figures/matching_ddm_simulation.json`。

### 图 1 的模拟设置

固定参数为 a = 1.224、z = 0.568、t = 0.464 秒、v = 2.214、dc = 0.321。每种刺激展示前两条路径，不按结局或外观筛选；路径从群体均值的乘积出发，在首次检测到边界接触时停止。图中时间从积累开始计，反应时另加 t。

路径使用单位扩散噪声的 Wiener 过程，以 Euler–Maruyama 步进（$`\Delta t=0.0001`$ 秒，种子 `2026102301`）近似模拟，并在最后一个时间步内插值确定边界接触。A 图不传播参数的后验不确定性或个体差异，仅展示 DDM 主体，不包含离群混合成分或实验反应窗口，因此不是完整的后验预测检验。$`\bar a\bar z`$ 是代入均值后的起点，不等于 $`E[az]`$。原有 PPC 与参数恢复限制仍适用。

## 参考资料

- [Hu et al. (2020), *Good me, bad me: Prioritization of the good-self during perceptual decision-making*](https://doi.org/10.1525/collabra.301)
- [Pan et al. (2025), dockerHDDM](https://doi.org/10.1177/25152459241298700)
- [HDDM 官方文档](https://hddm.readthedocs.io/en/latest/)
- [Stimulus-coding 方法代码](https://github.com/Chuan-Peng-Lab/ddm_stim_coding_posts)
- [Myers et al. (2022), A practical introduction to using the drift diffusion model](https://doi.org/10.3389/fpsyg.2022.1039172)
