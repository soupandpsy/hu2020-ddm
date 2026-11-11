# 分析输入

输入来自仓库已有的 Matching 数据与逐试次审计。`load_data()` 核验下列已发布来源的哈希及行对应关系，再应用 RT≥0.2 秒规则；不重新制定被试排除规则。

- [Part 1 原始记录](../../01_reproduction_2020/Data/raw/confirm/MS_rep_matchingTask_raw.csv)
- [Part 1 HDDM 输入](../../01_reproduction_2020/Data/matching_input.csv)
- [已核验的 Part 2 输入](../../02_reanalysis_2026/Data/matching_input.csv)
- [来源与按键映射审计](../../02_reanalysis_2026/Output/data_audit.json)

`processed/` 保存同一版本的 24,499 条有效响应、固定训练集 19,580 条、固定测试集 4,919 条和划分指纹。两个模型均只用训练集拟合、计算 DIC 和 PPC；测试集未用于本次比较。原划分保持不变，只移除了 RT<0.2 秒的记录。

42 名被试各有 Match/Mismatch × Good/Bad × Self/Other 八个条件。RT 单位为秒；实际范围 0.204–1.1685 秒。正确及错误反应均保留，无额外 RT 上限。`response=1/0` 表示 Match/Mismatch 选择，`stim=1/0` 表示 Match/Mismatch 刺激；正确率由二者是否一致计算，不能用上边界选择比例替代。

`source_row` 与 `event_id` 保留来源对应关系。Timeout 和非法按键不在有效响应输入中；当前模型未处理其观察机制，也未对 RT 下限筛选作条件似然校正。

划分清单中的状态字段记录 2026-11-07 制作数据时的情况；2026-11-11 的四链结果见各模型的 provenance 文件，不能把早期的 `new_mcmc_started: false` 理解为当前尚未拟合。
