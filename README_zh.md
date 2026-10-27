# Hu et al.（2020）：复现、机制再分析与研究进展

[English](README.md) · 状态：2026-10-28

本仓库围绕 Good-self 行为优势展开：先核对原论文行为结果及历史 HDDM 探索，再比较 v / dc 竞争解释。它是独立的探索性研究，不是原研究的官方更新。

| 部分 | 问题、阅读入口 | 当前状态 |
|---|---|---|
| 一：行为复现与历史 HDDM | [README](01_reproduction_2020/README_zh.md) → [复现 Notebook](01_reproduction_2020/reproduction_results.ipynb) → [历史 M0–M4](02_reanalysis_2026/results.ipynb) | 选定统计复现；模型设定差异、未收敛分支和补充分析缺口保留 |
| 二：H0–H3 v/dc 竞争机制 | [README](02_reanalysis_2026/README_zh.md) → [综合 Notebook](02_reanalysis_2026/matching_h0_h3_results.ipynb) | 全数据拟合与 PPC 完成；可靠五折预测比较尚未完成 |
| 三：后续研究 | 待明确方案和开展正式分析 | **尚未正式分析，无第三部分结果** |

阅读主线：研究问题 → 方法选择 → 诊断发现 → 方法调整 → 结果 → 限制 → 下一步。第一部分支持选定行为现象的复现，不等于原机制解释已获独立验证；第二部分四模型采样均通过但 PPC 均有失配，不能宣布 H2 推翻 v 解释或 H3 唯一胜出。五折为四模型 × 五折的 20 套训练拟合，修复仅一套完成，暂无可靠 ΔELPD 排名。

文件保持 `Data/`（输入）、`Scripts/`（代码）、`Output/`（原位结果）；历史 M0–M4 因路径兼容保留在 `02_reanalysis_2026`，阅读上归入第一部分背景。没有搬动或删除原始数据与历史结果。[目录与命名说明](02_reanalysis_2026/docs/ARCHIVE_AND_NAMING.md)区分 B0–B3 冒烟测试、M0–M4 探索、H0–H3 竞争模型，避免将真实标识误改。

[第二部分 README](02_reanalysis_2026/README_zh.md)提供方法、图表索引和轻量复现命令；[验证记录](02_reanalysis_2026/docs/VALIDATION.md)说明检查范围。已保存 Notebook 图表可直接阅读；完整后验和大型采样中间文件未发布，不能把展示复现当作采样全流程复现。

数据与原材料来自[作者仓库](https://github.com/hcp4715/moralSelf_ddm)；历史刺激编码实现参考[方法仓库](https://github.com/Chuan-Peng-Lab/ddm_stim_coding_posts)。
