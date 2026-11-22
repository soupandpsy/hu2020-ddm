# Hu et al. (2020) 正式实验行为分析

library(dplyr)
library(tidyr)
library(ggplot2)

raw_dir <- "data/raw"
derived_dir <- "data/derived"
results_dir <- "outputs/behavior"

# ---- 1. 读取原始试次 ----

matching_raw <- read.csv(file.path(raw_dir, "MS_rep_matchingTask_raw.csv"))
categ_raw <- read.csv(file.path(raw_dir, "MS_rep_categTask_raw.csv"))

# ---- 2. 筛选被试和试次 ----

# 返回试次数不符、某条件全错或总体正确率低于 50% 的被试编号。
subjects_to_exclude <- function(trials, condition, expected_trials) {
  by_condition <- trials |>
    group_by(Subject, across(all_of(c(condition, "Morality", "Identity")))) |>
    summarise(n_trials = n(), n_correct = sum(ACC == 1), .groups = "drop")

  invalid_conditions <- by_condition |>
    filter(n_trials != expected_trials | n_correct == 0) |>
    pull(Subject)

  low_accuracy <- by_condition |>
    group_by(Subject) |>
    summarise(accuracy = sum(n_correct) / sum(n_trials), .groups = "drop") |>
    filter(accuracy < 0.5) |>
    pull(Subject)

  union(invalid_conditions, low_accuracy)
}

# Matching 每条件应有 75 次，Categorization 每条件应有 90 次。
matching_excluded <- subjects_to_exclude(matching_raw, "Match", 75)
categ_excluded <- subjects_to_exclude(categ_raw, "Task", 90)

matching_retained <- filter(matching_raw, !Subject %in% matching_excluded)
categ_retained <- filter(categ_raw, !Subject %in% categ_excluded)

# 生成行为分析试次。无响应和无效键计为错误，剔除 RT ≤ 200 ms 的正确试次。
prepare_behavior_trials <- function(trials) {
  trials |>
    mutate(
      correct = ACC == 1,
      Valence = factor(Morality, levels = c("Good", "Bad")),
      Identity = factor(Identity, levels = c("Self", "Other"))
    ) |>
    filter(!(correct & RT <= 0.2))
}

matching_behavior <- prepare_behavior_trials(matching_retained)
categ_behavior <- prepare_behavior_trials(categ_retained)

# ---- 3. 计算每位被试的条件指标 ----

# 平均 RT 取正确匹配试次，d′ 取命中率和虚报率的正态分位数之差。
# 不匹配试次中，无响应、无效键和答错均计入虚报数。
matching_cells <- matching_behavior |>
  group_by(Subject, Valence, Identity) |>
  summarise(
    rt_ms = 1000 * mean(RT[Match == "match" & correct]),
    n_match = sum(Match == "match"),
    n_mismatch = sum(Match == "mismatch"),
    hit_rate = sum(Match == "match" & correct) / n_match,
    false_alarm_rate = sum(Match == "mismatch" & !correct) / n_mismatch,
    .groups = "drop"
  ) |>
  # 用半个试次修正命中率为 1 或虚报率为 0 的条件。
  mutate(
    hit_rate = if_else(hit_rate == 1, 1 - 0.5 / n_match, hit_rate),
    false_alarm_rate = if_else(false_alarm_rate == 0, 0.5 / n_mismatch, false_alarm_rate),
    dprime = qnorm(hit_rate) - qnorm(false_alarm_rate)
  ) |>
  select(Subject, Valence, Identity, rt_ms, dprime)

# 分 Id 和 Val 汇总正确反应的平均 RT 与全部保留试次的正确率。
categ_cells <- categ_behavior |>
  group_by(Subject, Task, Valence, Identity) |>
  summarise(
    rt_ms = 1000 * mean(RT[correct]),
    accuracy = mean(correct),
    .groups = "drop"
  )

# 合并 Id 和 Val 的试次，计算配对对比使用的 RT 与正确率。
categ_collapsed <- categ_behavior |>
  group_by(Subject, Valence, Identity) |>
  summarise(
    rt_ms = 1000 * mean(RT[correct]),
    accuracy = mean(correct),
    .groups = "drop"
  )

# ---- 4. 生成 ANOVA 表与贝叶斯模型比较表 ----

# 对一个行为指标拟合重复测量 ANOVA 和贝叶斯 ANOVA，返回两张表。
analyze_anova <- function(cells, measure, within) {
  data <- cells |>
    select(Subject, all_of(within), value = all_of(measure)) |>
    mutate(across(all_of(c("Subject", within)), factor)) |>
    as.data.frame()

  fit <- afex::aov_ez(
    id = "Subject", dv = "value", data = data, within = within, type = 3,
    anova_table = list(correction = "none", es = "none")
  )

  tab <- fit$anova_table
  anova <- tibble::tibble(
    effect = rownames(tab),
    df_effect = tab[["num Df"]],
    df_error = tab[["den Df"]],
    F = tab[["F"]],
    p = tab[["Pr(>F)"]],
    error_ms = tab[["MSE"]]
  )

  # 用效应、残差与被试平方和计算 ω²，负估计截为 0。
  subject_means <- tapply(data$value, data$Subject, mean)
  n_subjects <- length(subject_means)
  subject_ms <- (nrow(data) / n_subjects) * var(subject_means)
  subject_ss <- (n_subjects - 1) * subject_ms
  effect_ss <- anova$F * anova$df_effect * anova$error_ms
  error_ss <- anova$df_error * anova$error_ms

  anova$omega2 <- pmax(
    0,
    (effect_ss - anova$df_effect * anova$error_ms) /
      (effect_ss + error_ss + subject_ss + subject_ms)
  )

  # 各候选模型均含被试随机效应。固定和随机效应的先验尺度为 0.5 和 1。
  formula <- reformulate(c(paste(within, collapse = " * "), "Subject"), response = "value")
  bf <- BayesFactor::anovaBF(
    formula, data = data, whichRandom = "Subject", whichModels = "withmain",
    rscaleFixed = 0.5, rscaleRandom = 1, iterations = 10000, progress = FALSE
  )
  # BF10 以仅含被试随机效应的模型为参照。
  bf_table <- BayesFactor::extractBF(bf)
  models <- tibble::tibble(
    model = c("Subject", rownames(bf_table)),
    BF10 = c(1, bf_table$bf)
  )

  # 等模型先验下，纳入 BF 是含该效应与不含该效应的平均 BF10 之比。
  terms <- strsplit(models$model, " + ", fixed = TRUE)
  anova$BF_inclusion <- vapply(anova$effect, function(effect) {
    included <- vapply(terms, function(model_terms) effect %in% model_terms, logical(1))
    mean(models$BF10[included]) / mean(models$BF10[!included])
  }, numeric(1))

  list(anova = anova, models = models)
}

set.seed(2020)

matching_design <- c("Valence", "Identity")
categ_design <- c("Task", "Valence", "Identity")

analyses <- list(
  matching_rt     = analyze_anova(matching_cells, "rt_ms", matching_design),
  matching_dprime = analyze_anova(matching_cells, "dprime", matching_design),
  categ_rt        = analyze_anova(categ_cells, "rt_ms", categ_design),
  categ_accuracy  = analyze_anova(categ_cells, "accuracy", categ_design)
)

anova_results <- bind_rows(lapply(analyses, function(x) x$anova), .id = "outcome")
bayes_models <- bind_rows(lapply(analyses, function(x) x$models), .id = "outcome")

# 将论文报告的 F、分母自由度和纳入 BF 合并为对照列。
paper_anova <- tibble::tribble(
  ~outcome,          ~effect,            ~paper_F, ~paper_df_error, ~paper_BF_inclusion,
  "matching_rt",     "Valence",             57.88,              40,              2.3e8,
  "matching_rt",     "Valence:Identity",    14.65,              41,              133.9,
  "matching_dprime", "Valence",              5.71,              41,                8.5,
  "matching_dprime", "Valence:Identity",    12.03,              41,               24.7,
  "categ_rt",        "Valence",             15.30,              40,              407,
  "categ_rt",        "Identity",            41.60,              40,              879324,
  "categ_rt",        "Valence:Identity",    41.60,              40,               54.5,
  "categ_accuracy",  "Valence:Identity",     8.76,              40,               14.7
)

anova_results <- left_join(anova_results, paper_anova, by = c("outcome", "effect"))

# ---- 5. 生成配对对比表 ----

# 按被试配对两个条件，计算均值、单侧与双侧 t 检验、d_z 和方向 Bayes factor。
paired_comparison <- function(cells, measure, first, second, direction) {
  wide <- cells |>
    select(Subject, Valence, Identity, value = all_of(measure)) |>
    unite("condition", Valence, Identity) |>
    pivot_wider(names_from = condition, values_from = value)

  x <- wide[[first]]
  y <- wide[[second]]
  difference <- x - y
  test <- t.test(x, y, paired = TRUE)
  directional_test <- t.test(x, y, paired = TRUE, alternative = direction)

  # 方向 BF 对 RT 检验负差值，对 d′ 与正确率检验正差值。
  alternative <- if (direction == "less") c(-Inf, 0) else c(0, Inf)
  bf <- BayesFactor::ttestBF(
    x, y, paired = TRUE, rscale = 0.707, nullInterval = alternative
  )

  tibble::tibble(
    mean_first = mean(x),
    mean_second = mean(y),
    mean_difference = mean(difference),
    t = unname(test$statistic),
    df = unname(test$parameter),
    p_two_sided = test$p.value,
    p_directional = directional_test$p.value,
    dz = mean(difference) / sd(difference),
    BF_directional = BayesFactor::extractBF(bf, onlybf = TRUE)[1]
  )
}

# 指定配对条件、检验方向和论文对照值。
contrast_plan <- tibble::tribble(
  ~outcome,          ~measure,   ~first,       ~second,      ~direction, ~paper_t, ~paper_BF,
  "matching_rt",     "rt_ms",    "Good_Self",  "Bad_Self",   "less",         -8.42,    1.19e8,
  "matching_rt",     "rt_ms",    "Good_Other", "Bad_Other",  "less",         -2.37,       4.01,
  "matching_rt",     "rt_ms",    "Good_Self",  "Good_Other", "less",         -3.34,      35.50,
  "matching_dprime", "dprime",   "Good_Self",  "Bad_Self",   "greater",       4.30,     472.80,
  "matching_dprime", "dprime",   "Good_Self",  "Good_Other", "greater",       2.67,       7.38,
  "categ_rt",        "rt_ms",    "Good_Self",  "Bad_Self",   "less",         -3.18,      23.90,
  "categ_rt",        "rt_ms",    "Good_Self",  "Good_Other", "less",         -4.37,      57.00,
  "categ_accuracy",  "accuracy", "Good_Self",  "Bad_Self",   "greater",       3.06,      17.90,
  "categ_accuracy",  "accuracy", "Good_Self",  "Good_Other", "greater",       3.77,     106.70
)

contrast_data <- list(
  matching_rt = matching_cells,
  matching_dprime = matching_cells,
  categ_rt = categ_collapsed,
  categ_accuracy = categ_collapsed
)

comparisons <- vector("list", nrow(contrast_plan))
for (i in seq_len(nrow(contrast_plan))) {
  contrast <- contrast_plan[i, ]
  comparisons[[i]] <- paired_comparison(
    contrast_data[[contrast$outcome]],
    contrast$measure, contrast$first, contrast$second, contrast$direction
  )
}

contrast_results <- bind_cols(contrast_plan, bind_rows(comparisons))

# ---- 6. 绘制 Figure 3A–C ----

# 将同一任务的 RT 和另一项指标画在两个并排子图中。
behavior_panel <- function(cells, measure, measure_label, title) {
  points <- cells |>
    transmute(Subject, Valence, Identity, rt_ms, performance = .data[[measure]]) |>
    pivot_longer(c(rt_ms, performance), names_to = "metric", values_to = "value") |>
    mutate(
      metric = factor(
        metric, levels = c("rt_ms", "performance"),
        labels = c("RT (ms)", measure_label)
      ),
      x = as.integer(Valence) + if_else(Identity == "Self", -0.12, 0.12)
    )

  # 以被试为单位计算每个条件的均值和标准误。
  means <- points |>
    group_by(metric, Identity, x) |>
    summarise(
      mean_value = mean(value),
      se = sd(value) / sqrt(n()),
      .groups = "drop"
    )

  # 灰线连接同一被试的条件点。彩线表示组均值，误差条为均值 ± 1.96 SE。
  ggplot(points, aes(x, value)) +
    geom_line(aes(group = Subject), colour = "grey50", alpha = 0.12) +
    geom_point(colour = "grey50", alpha = 0.16) +
    geom_line(data = means, aes(y = mean_value, colour = Identity, group = Identity)) +
    geom_point(data = means, aes(y = mean_value, colour = Identity), size = 2.8) +
    geom_errorbar(
      data = means,
      aes(y = mean_value, ymin = mean_value - 1.96 * se,
          ymax = mean_value + 1.96 * se, colour = Identity),
      width = 0.06
    ) +
    scale_x_continuous(breaks = 1:2, labels = c("Good", "Bad")) +
    scale_colour_brewer(palette = "Dark2") +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    labs(title = title, x = NULL, y = NULL, colour = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
}

figures <- list(
  Figure_3A = behavior_panel(matching_cells, "dprime", "d-prime", "A. Matching"),
  Figure_3B = behavior_panel(filter(categ_cells, Task == "Id"), "accuracy",
                             "Accuracy", "B. Identity categorization"),
  Figure_3C = behavior_panel(filter(categ_cells, Task == "Val"), "accuracy",
                             "Accuracy", "C. Valence categorization")
)

# ---- 7. 生成三份 HDDM 试次输入表 ----

# 从保留被试的原始试次中选取 ACC 为 0 或 1 的记录，保留快速反应。RT 以秒计。
# stim 记录正确边界，response 记录实际选择。答错时 response 取另一边。
# Matching 的上边界为 Match（stim = 1）。
matching_hddm <- matching_retained |>
  filter(ACC %in% c(0, 1)) |>
  transmute(
    subj_idx = Subject, match = Match, val = Morality, id = Identity,
    stim = as.integer(Match == "match"),
    response = if_else(ACC == 1, stim, 1L - stim),
    rt = RT
  )

# 将指定形状映射到上边界，生成一个分类任务的试次输入表。
encode_categorization <- function(trials, task_code, upper_shapes) {
  trials |>
    filter(Task == task_code, ACC %in% c(0, 1)) |>
    transmute(
      subj_idx = Subject, val = Morality, id = Identity,
      stim = as.integer(Shape %in% upper_shapes),
      response = if_else(ACC == 1, stim, 1L - stim),
      rt = RT
    )
}

# Identity 的上边界为 Self，Valence 的上边界为 Good。
identity_hddm <- encode_categorization(categ_retained, "Id", c("moralSelf", "immoralSelf"))
valence_hddm <- encode_categorization(categ_retained, "Val", c("moralSelf", "moralOther"))

# ---- 8. 导出表和图 ----

# 三张统计表，CSV。
write.csv(anova_results, file.path(results_dir, "anova.csv"), row.names = FALSE)
write.csv(bayes_models, file.path(results_dir, "bayes_models.csv"), row.names = FALSE)
write.csv(contrast_results, file.path(results_dir, "planned_contrasts.csv"), row.names = FALSE)

# Figure 3A、3B、3C 各一个 PDF。
for (name in names(figures)) {
  ggsave(
    filename = file.path(results_dir, paste0(name, ".pdf")),
    plot = figures[[name]], width = 8, height = 4, units = "in"
  )
}

# 三份 HDDM 试次输入表，CSV。
write.csv(matching_hddm, file.path(derived_dir, "matching_hddm.csv"), row.names = FALSE)
write.csv(identity_hddm, file.path(derived_dir, "identity_hddm.csv"), row.names = FALSE)
write.csv(valence_hddm, file.path(derived_dir, "valence_hddm.csv"), row.names = FALSE)
