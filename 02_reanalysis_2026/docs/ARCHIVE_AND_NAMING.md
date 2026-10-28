# 发布结构与历史归档

正式结果只保留两份：[第一部分](../../01_reproduction_2020/reproduction_results.ipynb)、[第二部分](../matching_h0_h3_results.ipynb)。第一部分含原论文行为复现和历史M0–M4；第二部分含行为诊断至H0–H3结果全链条。第三部分未正式分析。

所有其他本地Notebook已移到仓库 `_local_archive/notebooks/`，保持原相对路径。该目录由Git忽略；本地manifest逐文件记录SHA-256。旧文件从Git当前版本移除，但仍可通过提交历史恢复。原数据、CSV/PNG结果、后验与数据库不搬动、不覆盖。

第二部分不再读取旧Notebook，早期保存表图独立于 `docs/presentation_outputs.json`；采样器所需原代码单元逐字保存于 `Scripts/mechanism_source_cells.json`，编号保持兼容。阶段报告脚本停止向历史源Notebook追加内容，仍生成原Markdown/JSON报告。完整采样没有在本次重跑。

B0–B3是历史冒烟测试，M0–M4是历史候选，H0–H3是当前竞争约束，均保留真实标识。fold4表示第五折；四模型×五折是20套训练拟合，当前未完成可靠预测排名。

两份 `README_historical_m0_m4`（中英文历史说明）已移至本机 `_local_archive/documents/` 并记录哈希；有用的复现说明已并入第一部分 README，历史结果全部从第一部分总 Notebook 阅读。
