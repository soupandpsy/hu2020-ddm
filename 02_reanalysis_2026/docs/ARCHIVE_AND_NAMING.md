# 目录、归档与命名核查（2026-10-28）

[当前 README](../README_zh.md) · [综合 Notebook](../matching_h0_h3_results.ipynb)

## 唯一当前入口与历史顺序

- 第一部分：`01_reproduction_2020/README*`、`reproduction_results.ipynb`；历史模型方法背景为 `02_reanalysis_2026/README_historical_m0_m4*`、`results.ipynb`。原路径保留，避免断开脚本、图与引用。
- 第二部分当前入口：`02_reanalysis_2026/README*` 和 `matching_h0_h3_results.ipynb`。图表与表格直接读取既有 Output，不另建一份竞争性分析结果。
- 第二部分过程记录：behavior_review → mechanism_simulation → mechanism_audit → observation_recovery；`matching_model_design.md` 是历史设计及阶段进度。所有这些文件首部标明历史范围，代码、单元次序和保存输出不变。
- 本地 Output 内的 `continuation_20261026`（恢复修正）、`diagnostic_H2H3_20261026_1425`（早期真实拟合诊断）、`continuation_stage3_20261026`（H2 PPC/H3拟合）、`continuation_stage4_20261027`（四模型 PPC/首轮五折）、`continuation_stage5_20261027`（修复/敏感性）按时间归档，不重命名、不删除。
- 其中旧 REPORT、PROGRESS 和 `matching_reanalysis_current.ipynb` 是当时快照；“未推送”“未拟合”等句不代表本次整理后状态。部分历史报告仅本地保留，发布证据以 manifest 白名单为准。未来自动报告脚本不会改写本次新增综合 Notebook。
- 第三部分未正式分析，不创建带结果暗示的第三部分目录。

## 命名核查结果

| 标识 | 内容核查 | 处理 |
|---|---|---|
| B0/B1/B2/B3 | `Output/synthetic_smoke_20261020_143349/status.json` 明确记录 synthetic smoke、40 个保存样本；同组存在 B3.db、B3_summary.csv、B3_traces.csv | 真实历史测试标识，保留；不是 H3，不机械替换 |
| M0–M4 | 旧 HDDM 候选的 z/dc 结构筛选 | 保留为历史探索标识，不改称 H0–H3 |
| H0–H3 | Good-self−Bad-self 的 log(v)/dc 竞争约束 | 当前主线保持一致 |
| fold4 | 从0开始的第五折 | 保留机器标识，正文解释“第五折” |
| “五十个比较” | 实际为四模型×五折=20套训练拟合 | 当前正文统一写“五折预测比较”，不声称50次已完成 |
| results.ipynb / 旧 Part II 标题 | 内容是旧 M0–M4，并非 H0–H3 最新结果 | 文件保留兼容，Notebook 标题改为“历史 M0–M4”；README 分离当前/历史入口 |

未发现有内容证据支持把某个 b3 文件认定为 H3 口误，因此本次没有强行改名原始产物。十六进制哈希中出现 b3 也不是模型命名问题。图 I/II 编号用于当前导航，原论文 Figure 编号、模型标签及原图标题不改写。

## 整理要求与后续维护

按“问题→方法→诊断→调整→结果→限制→下一步”维护 README 与综合 Notebook。新结论必须指向保存的表/诊断；失败与未完成结果不可隐藏。每次采样另建输出目录，归档旧状态并更新唯一当前入口。不要重排 observation_recovery 单元：现有机制脚本依赖其编号。Git 只发布明确白名单；凭据、数据库、大型后验及逐事件密度保留本地。保留研究原始数据与有价值历史材料。

## 总 Notebook 合并更新

第二部分所有阶段已合入 `matching_h0_h3_results.ipynb`，包含目录与章内解释、早期保存表图以及正式结果。README默认只引导这一份，历史来源集中在总Notebook文末。原四份过程Notebook原位保留，避免破坏按固定路径/单元编号提取模型的代码；没有继续新增并列综合报告。
