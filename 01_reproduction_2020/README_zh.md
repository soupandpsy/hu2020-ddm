# 第一部分：原论文复现与历史 HDDM

[English](README.md) · [项目首页](../README_zh.md) · [本部分总 Notebook](01_reproduction_2020.ipynb)

## 研究范围

本部分复现《Good me, bad me: Prioritization of the good-self during perceptual decision-making》的主要计算分析，涵盖 Pilot 与 Confirm 行为统计、Pilot ex-Gaussian 和 Confirm HDDM，并单独保存敏感性分析。结果总体上复现了部分选定统计量和模型分析，但并非所有统计量、补充分析或模型设置均与论文完全一致。

原始材料见[作者仓库](https://github.com/hcp4715/moralSelf_ddm)。

## 复现结果

- Pilot 保留 29 名被试；Confirm Matching 保留 42 名，Categorization 保留 41 名。R 导出的 Confirm 输入与仓库内输入一致。15 项选定的 t 统计量在指定容差内一致（Pilot < .01；Confirm < .02），不代表严格符合论文打印值的舍入结果。
- Matching/Identity/Valence vtz p0 和 Identity vz p0 通过了配置中的四次运行诊断标准。Matching/Valence vz 只有一次运行，仅用于描述。**Identity vz p0.05 未通过诊断，不应合并其各次样本。**
- Figure 3C 将作者误用的 Identity 汇总改为 Valence 数据；误差线仍为均值 ± 1.96 SE。Pilot 个体线在 Subject × Identity 内连接 Good/Bad，Confirm 则连接每名被试的四个条件。
## 差异与局限

- Pilot Categorization RT 与作者 JASP 的条件列分配不同，对应七项 F/omega 不能视为等价比较。其余 33 项可比 ANOVA F 值在 1e-9 内一致。Omega 定义、包含贝叶斯因子及 Pilot Categorization Accuracy 的干扰因素设置也存在差异，因此不能声称所有统计量与论文相同。
- ex-Gaussian 使用 1,000 次 bootstrap，而论文报告 2,000 次。一个单元（7014，Val/Other/Bad，仅 3 个试次）未拟合；231 个完成拟合的单元筛选后得到 27 名符合条件的被试。与可取得的作者估计相比，mu/sigma/tau 最大差异约为 .0283/.0567/.0208 秒；不声称 ANOVA 等价。详见[各单元结果](Output/behavior/pilot/exgaussian/parameters_all.csv)及保留的警告记录。
- 补充 ex-Gaussian 贝叶斯与简单效应分析、Figure S3、完整 Table S2/S3 汇总和贝叶斯因子稳健性图尚未全部完成。论文未报告 Pilot DDM，本项目也未实现该部分。
- 本部分记录 2020 复现工作。上述模型通过的是当时配置中的诊断标准，不能自动等同于 2026 再分析采用的秩归一化 R-hat 和 ESS 标准。2026 再分析及其独立结论见[第二部分](../02_reanalysis_2026/README_zh.md)。


## 历史探索与阅读顺序

总 Notebook 按行为复现（I-1–I-5）→ 原论文 HDDM（I-6）→ M0–M4 扩展的路径示意、参数、轨迹与PPC（I-7–I-10）→ 限制展开。M0–M4 属于历史模型，代码和结果因路径依赖仍在 `../02_reanalysis_2026/`。Matching M4 通过诊断，但慢RT预测失配；Identity/Valence长链与四种恢复情景仍有收敛不足。行为复现不能直接证明唯一机制。

截至2026-10-28，第二部分H0–H3全数据拟合均通过当前采样诊断，但PPC存在RT、超时率及快速Good-self优势失配，五折比较尚未可靠完成。不能声称H2推翻v解释，或H3已是最佳/唯一机制。第三部分尚未正式开展。

## 必要复现方法

轻量查看：从仓库根目录安装 `requirements-notebooks.txt` 后执行：

```bash
python -m jupyter nbconvert --execute --to notebook --output-dir .preview 01_reproduction_2020/01_reproduction_2020.ipynb
```

数据口径见[Data说明](Data/README.md)。以下原论文复现命令从 `01_reproduction_2020/` 运行；R分析会产生新输出，绘图命令默认读取已保存结果。

```bash
Rscript --vanilla Scripts/pilot_behavior/pilot_analysis.R
Rscript --vanilla Scripts/confirm_behavior/confirm_analysis.R
Rscript --vanilla Scripts/pilot_behavior/pilot_plot.R --analysis-dir Output/behavior/pilot
Rscript --vanilla Scripts/confirm_behavior/confirm_plot.R --analysis-dir Output/behavior/confirm
python -B Scripts/hddm/figure04_plot.py
```

原论文计算环境为R 4.6.0、Python 3.7.12、HDDM 0.8.0；R包版本见结果中的 `dependencies.csv` 与 `sessionInfo.txt`。完整HDDM拟合使用 `Scripts/hddm/<task>/<model>/` 下的Run入口与JSON配置；例如Matching：

```bash
python -B Scripts/hddm/matching/vtz/matching_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
```

`CONFIRM_DIR`替换为Confirm分析输出；省略该参数则读取内置输入。`vtz`允许条件v/t并估计z，`vz`共享t；`p0/p0p05`为污染比例，`r01`为运行编号。完整数据库未发布，历史图4仅使用35个群体参数各9,000样本，不能代替多链诊断。

历史M0–M4使用其任务配置中的dockerHDDM环境，从 `02_reanalysis_2026/` 运行：

```bash
python Scripts/hddm/reproduce.py plan --output Output_new
python Scripts/hddm/reproduce.py run --output Output_new --workers 2
python Scripts/hddm/reproduce.py verify --output Output_new
```

完整采样昂贵且未在本轮重跑；展示Notebook只读取已存结果。后续问题与结果见[第二部分](../02_reanalysis_2026/README_zh.md)。
