# 道德自我关联与知觉决策：HDDM再分析

[English](README.md) · [完整结果与计算过程](results.ipynb)

## 研究问题与方法

本分析使用Hu等（2020）的验证性研究数据，考察与“好自我”关联的刺激是否比“坏自我”和“好他人”具有更快的证据积累速度，并区分证据积累与反应偏置。

实验条件由身份（Self、Other）和道德属性（Good、Bad）组合而成。Matching任务要求判断形状与标签是否匹配；Identity和Valence任务分别要求判断身份和道德属性。三个任务分别纳入42、41、41名被试，有效试次为24527、14493、14431。Matching排除了一个非法按键试次。

每个任务比较五个模型（M0–M4），随后对候选M4进行四链评估。M4估计各条件的漂移率v和漂移偏置dc，以及跨条件共享的起点z、边界间距a和非决策时间t。Matching中匹配与不匹配刺激共享v的幅度。初筛每个模型只有一条链，模型选择仍有不确定性。

## Matching结果

Good-self的漂移率高于Bad-self，也高于Good-other。自我条件中的Good−Bad差异大于他人条件，表明道德属性与身份的关系体现在证据积累速度上。

|漂移率对比|均值差|95%可信区间|
|---|---:|---|
|Good-self − Bad-self|0.666|[0.387, 0.942]|
|Good-self − Good-other|0.438|[0.160, 0.719]|
|Good-other − Bad-other|−0.009|[−0.278, 0.267]|
|身份×道德属性交互¹|0.675|[0.284, 1.059]|

¹ (Good-self−Bad-self)−(Good-other−Bad-other)。区间为95%等尾后验可信区间。

起点z为0.568，95%可信区间[0.558, 0.578]，表明存在朝向“匹配”反应的起点偏置。Good-self与Bad-self的dc差为0.600，[0.427, 0.774]。两类偏置描述的是匹配判断中的反应倾向，不能直接等同于对“好自我”的偏好。

![Matching参数估计与条件对比](Output/hddm/matching/figures/matching_m4_parameter_contrasts.png)

**图1. Matching的参数估计与条件对比。** 点为后验均值，线段为95%可信区间。v表示证据积累速度，dc表示漂移偏置，z表示相对起点。

## 估计稳定性

Matching的四条链得到接近的参数分布，收敛检查通过。Identity和Valence在增加迭代次数后仍有明显链间差异，因此不报告两项任务的参数解释。

|任务|最大Rhat|最低bulk ESS|最低tail ESS|判断|
|---|---:|---:|---:|---|
|Matching|1.0022|917|472|通过|
|Identity（长链）|1.3076|10|17|未通过|
|Valence（长链）|1.5282|7|27|未通过|

Rhat用于检查不同链的一致性，ESS表示采样包含的有效信息量。本分析采用Rhat<1.01、bulk/tail ESS均≥400，并结合轨迹检查。

![Matching四链轨迹与分布](Output/hddm/matching/figures/matching_m4_traces.png)

**图2. Matching的四链轨迹与参数分布。** 不同颜色表示不同链，左列为采样轨迹，右列为参数分布。

## 模型检验与局限

模型预测的慢反应时尾部比观测数据更长（图3）。实验有约0.9–1.2秒的反应窗口，而预测未施加相同窗口，两者的可比性受到限制。逐试次窗口尚未恢复，因此还不能判断窗口对差异的贡献。

![Matching后验预测检验](Output/hddm/matching/figures/matching_m4_ppc.png)

**图3. Matching的观测数据与后验预测。** 黑色叉号为观测值；彩色点和区间表示预测均值及95%预测统计量区间。各面板分别显示选择比例及反应时分位数。两种补充离群生成约定用于检查预测对设定的敏感性，与拟合时的完整离群似然并不严格一致。

排除RT≤0.2秒后，主要群体估计变化很小，但一个个体参数的tail ESS为342，未达到预设门槛。四个模拟情景的参数恢复分析也均有收敛不足，尚不能确认参数恢复表现。

因此，Matching结果提供了与Good-self加工优势方向一致的证据，但解释仍受模型拟合与恢复检验的限制。本次模型与2020模型的参数约束不同，数值不宜直接比较；Identity和Valence的结果也不足以评价原文的跨任务结论。

## 文件与运行

[Notebook](results.ipynb)包含模型比较、条件估计及检查结果；数据位于`Data/`，计算脚本位于`Scripts/hddm/`，结果位于`Output/`。Notebook读取已保存结果，不启动模型拟合。

在任务配置指定的dockerHDDM环境中，从本目录运行：

```sh
python Scripts/hddm/reproduce.py plan --output Output_new
python Scripts/hddm/reproduce.py run --output Output_new --workers 2
python Scripts/hddm/reproduce.py verify --output Output_new
```

完整运行需要相邻的`01_reproduction_2020/Data`和一个尚不存在的输出目录。[输出清单](Scripts/hddm/output_manifest.json)记录各文件的生成入口。后处理已使用历史样本验证；整理后的完整拟合流程尚未从零重跑。数据库与逐次样本不随仓库分发。

### 任务入口与历史来源

每个任务保留配置、模型定义、运行、分析和绘图入口；Matching 另有恢复和敏感性分析入口。`run fit` 拟合单链，其余 Run 阶段调度 smoke、screening、final 或独立长链。Analyze 的三个操作可单独执行：

```sh
python Scripts/hddm/matching/matching_run.py fit --model m4 --chain 1 --phase screening --output Output_new --plan
python Scripts/hddm/matching/matching_analyze.py compare --output Output_new --plan
python Scripts/hddm/matching/matching_analyze.py diagnostics --output Output_new --plan
python Scripts/hddm/matching/matching_analyze.py ppc --output Output_new --plan
```

将 `matching` 替换为 `identity` 或 `valence` 即可使用相应入口。去掉 `--plan` 才会执行，`--help` 显示各操作选项。正式拟合在指纹一致的初筛 r01 后追加 r02–r04；长链重拟合使用四条独立链。Compare 读取已有初筛输出，不重新运行 PPC。

比较表中的 `status` 和 `chains` 仅描述初筛；`selection_status` 记录暂定选择，`followup_status` 区分仅初筛（`screening_only`）、通过收敛诊断（`diagnostics_pass`）和未收敛（`nonconverged`）。通过收敛诊断不代表预测或参数恢复检验通过。

重构前公开版本为提交 `639aa66af934db39b9e2962a6c5154350266758e`，历史标签为 `reanalysis-before-entry-refactor-20261023`。它标记已公开的仓库快照，不代表全部历史拟合均由该版代码产生；实际拟合来源以原始 fingerprint 为准。清单中的 `producer` 指向当前入口，`historical_snapshot` 指向先前公开结果。链来源路径相对于 `Output/`，数据库不随仓库分发。新拟合会记录变化后的 `common.py` 哈希，不能追加到指纹不同的历史链；新拟合使用独立输出目录，历史续跑则需使用与原指纹一致的代码。

## 参考资料

[Hu等（2020）](https://doi.org/10.1525/collabra.301)；[Pan等（2025），dockerHDDM](https://doi.org/10.1177/25152459241298700)；[刺激编码方法代码](https://github.com/Chuan-Peng-Lab/ddm_stim_coding_posts)，commit `468965c8d4c991bb7633f3f38257a6bac978bce5`。
