# Hu et al. (2020): pilot behavioral analysis.
# Figures are generated separately by pilot_plot.R.
# Required packages must already be installed; .rlib is optional and local.
SCRIPT <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), mustWork=TRUE)
ROOT <- dirname(dirname(dirname(SCRIPT)))
args <- commandArgs(TRUE)
arg <- function(name, default=NULL) {
  i <- match(name,args)
  if (is.na(i)) return(default)
  if (i==length(args) || startsWith(args[i+1],"--")) stop("Missing value for ", name)
  args[i+1]
}
allowed <- c("--data-dir","--output-dir","--analysis-dir")
if (length(args) && (length(args) %% 2 != 0 || any(!args[seq(1,length(args),by=2)] %in% allowed))) {
  if(length(args)) stop("Use --data-dir, --output-dir or --analysis-dir followed by a path")
}
LIB <- file.path(ROOT,".rlib")
.libPaths(c(LIB,.libPaths()))
new_output <- function(path) {
  if(dir.exists(path) && length(list.files(path,all.files=TRUE,no..=TRUE))) stop("Output directory is not empty: ",path)
  if(file.exists(path) && !dir.exists(path)) stop("Output is a file: ",path)
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  normalizePath(path,mustWork=TRUE)
}
stamp <- function() format(Sys.time(),"%Y%m%d_%H%M%OS6")
set.seed(20261006)

STUDY <- "pilot"
DATA <- normalizePath(arg("--data-dir",file.path(ROOT,"Data","raw",STUDY)),mustWork=TRUE)
PKGS <- c("Rcpp","dplyr","tidyr","ggplot2","afex","BayesFactor","effectsize","retimes")
missing <- PKGS[!vapply(PKGS,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing R packages: ",paste(missing,collapse=", "),". Install in .rlib or your configured library.")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(afex);library(BayesFactor)})
options(contrasts=c("contr.sum","contr.poly"),BFprogress=FALSE)
OUT <- new_output(arg("--output-dir",file.path(ROOT,"Output","behavior",STUDY,paste0("analysis_",stamp()))))
globalCallingHandlers(warning=function(w) cat(conditionMessage(w), "\n", file=file.path(OUT,"warnings.txt"), append=TRUE))
for(n in c("hddm_input","exgaussian")) dir.create(file.path(OUT,n),showWarnings=FALSE)
find_input <- function(name) { p<-file.path(DATA,name); if(!file.exists(p)) stop("Missing raw input: ",p); normalizePath(p,mustWork=TRUE) }
save_csv <- function(x,name,dir=OUT) write.csv(x,file.path(dir,name),row.names=FALSE)
warn_expected <- function(label,actual,expected) {if(!identical(as.integer(actual),as.integer(expected))) stop(label,": ",actual," != ",expected)}

# 1. Statistical helpers
# ==============================================================================

dprime_value <- function(hit, miss, fa, cr) {
  hr <- hit / (hit + miss)
  fr <- fa  / (fa + cr)
  # 忠实于作者代码：只修正 hit rate == 1 和 false-alarm rate == 0
  if (hr == 1) hr <- 1 - 1 / (2 * (hit + miss))
  if (fr == 0) fr <- 1 / (2 * (fa + cr))
  qnorm(hr) - qnorm(fr)
}

dz_ci <- function(t_value, df, n) {
  out <- tryCatch(
    sort(vapply(c(.975, .025), function(q) {
      uniroot(function(ncp) pt(t_value, df, ncp = ncp) - q,
              interval = c(-100, 100))$root / sqrt(n)
    }, numeric(1))),
    error = function(e) c(NA_real_, NA_real_)
  )
  names(out) <- c("low", "high")
  out
}

planned_tests <- function(long4, outcome) {
  wide <- long4 |>
    select(Subject, Val, Id, y) |>
    tidyr::unite(cond, Val, Id, sep = "_") |>
    tidyr::pivot_wider(names_from = cond, values_from = y)

  pairs <- list(
    GS_BS = c("Good_Self",  "Bad_Self"),
    GO_BO = c("Good_Other", "Bad_Other"),
    GS_GO = c("Good_Self",  "Good_Other"),
    BS_BO = c("Bad_Self",   "Bad_Other"),
    GS_BO = c("Good_Self",  "Bad_Other")  # supplementary/congruency check
  )

  direction <- if (grepl("RT", outcome)) "less" else "greater"

  bind_rows(lapply(names(pairs), function(nm) {
    a <- wide[[pairs[[nm]][1]]]
    b <- wide[[pairs[[nm]][2]]]
    ok <- complete.cases(a, b); a <- a[ok]; b <- b[ok]
    z <- a - b; n <- length(z)
    tt2 <- t.test(a, b, paired = TRUE)
    tt1 <- t.test(a, b, paired = TRUE, alternative = direction)
    ci  <- dz_ci(unname(tt2$statistic), n - 1, n)

    bf2 <- as.numeric(BayesFactor::extractBF(
      BayesFactor::ttestBF(a, b, paired = TRUE, rscale = .707),
      onlybf = TRUE
    ))
    bfd <- as.numeric(BayesFactor::extractBF(
      BayesFactor::ttestBF(
        a, b, paired = TRUE, rscale = .707,
        nullInterval = if (direction == "less") c(-Inf, 0) else c(0, Inf)
      ),
      onlybf = TRUE
    ))[1]

    data.frame(
      outcome = outcome, contrast = nm, n = n,
      mean_1 = mean(a), mean_2 = mean(b),
      mean_difference = mean(z),
      t = unname(tt2$statistic), df = n - 1,
      p_two_sided = tt2$p.value,
      p_one_sided = tt1$p.value,
      one_sided_direction = direction,
      dz = mean(z) / sd(z),
      dz_low = ci["low"], dz_high = ci["high"],
      BF10_two_sided = bf2,
      BF_directional = bfd
    )
  }))
}

bayes_inclusion <- function(bf, within, outcome) {
  tab <- BayesFactor::extractBF(bf)
  tab$model <- rownames(tab)
  logs <- c(0, bf@bayesFactor$bf)
  names(logs) <- c("Subject", tab$model)
  w <- exp(logs - max(logs)); w <- w / sum(w)

  canon <- function(x) paste(sort(strsplit(x, ":", fixed = TRUE)[[1]]), collapse = ":")
  terms_by_model <- lapply(strsplit(gsub(" ", "", names(logs)), "\\+"),
                           function(z) unname(vapply(z, canon, "")))
  effects <- attr(terms(as.formula(paste("~", paste(within, collapse = "*")))),
                  "term.labels")

  bind_rows(lapply(effects, function(e) {
    ec <- canon(e)
    inc <- vapply(terms_by_model, function(z) ec %in% z, logical(1))
    prior <- mean(inc); post <- sum(w[inc])
    allbf <- (post / (1 - post)) / (prior / (1 - prior))

    matched <- lapply(which(inc), function(i) {
      target <- setdiff(terms_by_model[[i]], ec)
      j <- which(vapply(terms_by_model, function(z) setequal(z, target), logical(1)))
      if (length(j)) c(i, j[1]) else NULL
    })
    matched <- Filter(Negate(is.null), matched)
    mbf <- if (length(matched)) {
      ii <- unique(vapply(matched, `[`, integer(1), 1))
      jj <- unique(vapply(matched, `[`, integer(1), 2))
      (sum(w[ii]) / sum(w[jj])) / (length(ii) / length(jj))
    } else NA_real_

    data.frame(outcome = outcome, effect = e,
               BFincl_all = allbf, BFincl_matched = mbf)
  }))
}

run_rm_anova <- function(long, within, outcome) {
  dat <- long
  for (v in c("Subject", within)) dat[[v]] <- factor(dat[[v]])

  fit <- afex::aov_ez(
    id = "Subject", dv = "y", data = dat, within = within,
    type = 3, return = "afex_aov",
    anova_table = list(correction = "none", es = "ges")
  )
  tab <- as.data.frame(fit$anova_table)
  tab$effect <- rownames(tab); rownames(tab) <- NULL
  tab$outcome <- outcome

  omega <- tryCatch(
    as.data.frame(effectsize::omega_squared(
      fit, partial = FALSE, ci = .95, alternative = "two.sided"
    )),
    error = function(e) data.frame(error = conditionMessage(e))
  )
  omega$outcome <- outcome

  form <- as.formula(paste("y ~", paste(within, collapse = "*"), "+ Subject"))
  bf <- BayesFactor::anovaBF(
    form, data = dat, whichRandom = "Subject", whichModels = "withmain",
    rscaleFixed = .5, rscaleRandom = 1,
    iterations = 10000, progress = FALSE
  )
  models <- BayesFactor::extractBF(bf)
  models$model <- rownames(models); rownames(models) <- NULL
  models$outcome <- outcome

  list(frequentist = tab, omega = omega, bayes_models = models,
       bayes_inclusion = bayes_inclusion(bf, within, outcome))
}

cross_correlations <- function(effects) {
  effect_names <- c("goodslf_badslf", "goodoth_badoth",
                    "goodslf_goodoth", "badslf_badoth")
  bind_rows(lapply(effect_names, function(eff) {
    bind_rows(lapply(c("Val", "Id"), function(task) {
      bind_rows(lapply(c("RT", "accuracy"), function(metric) {
        xname <- paste0(if (metric == "RT") "RT_" else "d_", eff)
        yname <- paste0(task, if (metric == "RT") "_RT_" else "_ACC_", eff)
        ok <- complete.cases(effects[, c(xname, yname)])
        x <- effects[[xname]][ok]; y <- effects[[yname]][ok]
        ct <- cor.test(x, y)
        data.frame(
          effect = eff, task = task, metric = metric, n = sum(ok),
          r = unname(ct$estimate),
          low = ct$conf.int[1], high = ct$conf.int[2],
          p = ct$p.value,
          BF10 = as.numeric(BayesFactor::extractBF(
            BayesFactor::correlationBF(x, y, rscale = 1), onlybf = TRUE
          ))
        )
      }))
    }))
  }))
}

# 2. Load raw data
# ==============================================================================

MATCH_FILE <- find_input("MS_matchTask_raw.csv")
CATEG_FILE <- find_input("MS_categTask_raw.csv")

m_raw <- read.csv(MATCH_FILE, stringsAsFactors = FALSE)
c_raw <- read.csv(CATEG_FILE, stringsAsFactors = FALSE)

m_raw <- m_raw |>
  rename(Subject = SubjectID, ACC = Accuracy) |>
  mutate(
    Morality = ifelse(Morality == "moral", "Good", "Bad"),
    Identity = ifelse(Identity == "self", "Self", "Other"),
    Match = factor(Match, levels = c("match", "nonmatch")),
    Morality = factor(Morality, levels = c("Good", "Bad")),
    Identity = factor(Identity, levels = c("Self", "Other"))
  )

# 作者原代码明确排除 importance-based blocks
c_raw <- c_raw |>
  filter(Task %in% c("self", "moral")) |>
  rename(Subject = SubjectID, ACC = Accuracy) |>
  mutate(
    Morality = ifelse(Morality == "moral", "Good", "Bad"),
    Identity = ifelse(Identity == "self", "Self", "Other"),
    Task = ifelse(Task == "self", "Id", "Val"),
    Task = factor(Task, levels = c("Id", "Val")),
    Morality = factor(Morality, levels = c("Good", "Bad")),
    Identity = factor(Identity, levels = c("Self", "Other"))
  )


# ==============================================================================
# 3. Pilot preprocessing
# ==============================================================================

# Rule 1: procedural failures (作者硬编码)
excl_procedure <- c(7, 8, 2027, 7035)
m0 <- m_raw |> filter(!Subject %in% excl_procedure)
c0 <- c_raw |> filter(!Subject %in% excl_procedure)

# Matching: 每个被试前 48 trials 是 practice，按原始行顺序删除
m0 <- m0 |>
  group_by(Subject) |>
  mutate(.trial_order = row_number()) |>
  filter(.trial_order > 48) |>
  select(-.trial_order) |>
  ungroup()

# Rule 2: overall accuracy < .50；no response (-1) 计为错误
subject_acc <- function(dat) {
  dat |>
    mutate(.acc = ifelse(ACC == 1, 1, 0)) |>
    group_by(Subject) |>
    summarise(ACC = mean(.acc), .groups = "drop")
}
excl_low_m <- subject_acc(m0) |> filter(ACC < .5) |> pull(Subject)
excl_low_c <- subject_acc(c0) |> filter(ACC < .5) |> pull(Subject)

m_valid <- m0 |> filter(!Subject %in% excl_low_m)
# 与作者代码一致：matching 低准确率被试也从 categorization 排除
c_valid <- c0 |> filter(!Subject %in% c(excl_low_m, excl_low_c))

warn_expected("Pilot matching valid N", n_distinct(m_valid$Subject), 29)
warn_expected("Pilot categorization valid N", n_distinct(c_valid$Subject), 29)

exclusions <- data.frame(
  rule = c(rep("procedural_failure", length(excl_procedure)),
           rep("matching_overall_ACC_lt_.5", length(excl_low_m)),
           rep("categorization_overall_ACC_lt_.5", length(excl_low_c))),
  Subject = c(excl_procedure, excl_low_m, excl_low_c)
)
save_csv(exclusions, "pilot_exclusions.csv")

sample_flow <- data.frame(
  stage = c("raw_matching", "raw_categorization_after_importance_removed",
            "after_procedural_exclusion", "final_matching", "final_categorization"),
  participants = c(n_distinct(m_raw$Subject), n_distinct(c_raw$Subject),
                   n_distinct(m0$Subject), n_distinct(m_valid$Subject),
                   n_distinct(c_valid$Subject))
)
save_csv(sample_flow, "pilot_sample_flow.csv")

# Traditional analysis:
# 1) no response -> incorrect
# 2) seconds -> milliseconds
# 3) pilot 原代码删除所有 RT <= 200 ms trials
m_trad <- m_valid |>
  mutate(ACC = ifelse(ACC == -1, 0, ACC), RT = RT * 1000) |>
  filter(RT > 200)

c_trad <- c_valid |>
  mutate(ACC = ifelse(ACC == -1, 0, ACC), RT = RT * 1000) |>
  filter(RT > 200)

save_csv(m_trad, "pilot_matching_clean_traditional.csv")
save_csv(c_trad, "pilot_categorization_clean_traditional.csv")

# HDDM inputs：从 valid raw data 出发，保留 seconds，不删除 fast correct trials
m_hddm_match <- m_valid |>
  filter(ACC %in% c(0, 1), Match == "match") |>
  transmute(subj_idx = Subject, val = Morality, id = Identity,
            rt = RT, response = ACC)
m_hddm_nonmatch <- m_valid |>
  filter(ACC %in% c(0, 1), Match == "nonmatch") |>
  transmute(subj_idx = Subject, val = Morality, id = Identity,
            rt = RT, response = ACC)

c_hddm_val <- c_valid |>
  filter(ACC %in% c(0, 1), Task == "Val") |>
  mutate(
    response = ifelse((trialType == "moral" & ACC == 1) |
                      (trialType == "immoral" & ACC == 0), 1, 0),
    stim = ifelse(trialType == "moral", 1, 0)
  ) |>
  transmute(subj_idx = Subject, val = Morality, id = Identity,
            stim, response, rt = RT)

c_hddm_id <- c_valid |>
  filter(ACC %in% c(0, 1), Task == "Id") |>
  mutate(
    response = ifelse((trialType == "self" & ACC == 1) |
                      (trialType == "other" & ACC == 0), 1, 0),
    stim = ifelse(trialType == "self", 1, 0)
  ) |>
  transmute(subj_idx = Subject, val = Morality, id = Identity,
            stim, response, rt = RT)

save_csv(m_hddm_match,    "MS_match_hddm.csv", file.path(OUT, "hddm_input"))
save_csv(m_hddm_nonmatch, "MS_mismatch_hddm.csv", file.path(OUT, "hddm_input"))
save_csv(c_hddm_val,      "MS_categ_val_hddm_stim.csv", file.path(OUT, "hddm_input"))
save_csv(c_hddm_id,       "MS_categ_id_hddm_stim.csv", file.path(OUT, "hddm_input"))

# ex-Gaussian input：作者代码使用 valid raw data，correct & RT > .2 s
c_exg <- c_valid |> filter(ACC == 1, RT > .2)
save_csv(c_exg, "MS_categ_exG.csv", file.path(OUT, "exgaussian"))


# ==============================================================================
# 4. Behavioral variables
# ==============================================================================

# Matching ACC
m_acc <- m_trad |>
  group_by(Subject, Match, Morality, Identity) |>
  summarise(N = n(), countN = sum(ACC), ACC = mean(ACC), .groups = "drop")

m_acc_w <- m_acc |>
  mutate(key = paste("ACC", Match, Morality, Identity, sep = "_")) |>
  select(Subject, key, ACC) |>
  pivot_wider(names_from = key, values_from = ACC)

# Matching d'
m_dp <- m_trad |>
  group_by(Subject, Age, Gender, Morality, Identity) |>
  summarise(
    hit  = sum(Match == "match"    & ACC == 1),
    miss = sum(Match == "match"    & ACC == 0),
    CR   = sum(Match == "nonmatch" & ACC == 1),
    FA   = sum(Match == "nonmatch" & ACC == 0),
    .groups = "drop"
  ) |>
  rowwise() |>
  mutate(dprime = dprime_value(hit, miss, FA, CR)) |>
  ungroup()

m_dp_w <- m_dp |>
  mutate(key = paste("d", Morality, Identity, sep = "_")) |>
  select(Subject, Age, Gender, key, dprime) |>
  pivot_wider(names_from = key, values_from = dprime)

# Matching RT：correct matching/mismatching 都汇总；主分析只用 match
m_rt <- m_trad |>
  filter(ACC == 1) |>
  group_by(Subject, Match, Morality, Identity) |>
  summarise(RT = mean(RT), .groups = "drop")

m_rt_w <- m_rt |>
  mutate(key = paste("RT", Match, Morality, Identity, sep = "_")) |>
  select(Subject, key, RT) |>
  pivot_wider(names_from = key, values_from = RT)

match_wide <- m_acc_w |>
  inner_join(m_dp_w, by = "Subject") |>
  inner_join(m_rt_w, by = "Subject")
save_csv(match_wide, "MS_match_behav_wide.csv")

# Categorization task-specific ACC / RT
c_acc <- c_trad |>
  group_by(Subject, Task, Morality, Identity) |>
  summarise(N = n(), countN = sum(ACC), ACC = mean(ACC), .groups = "drop")

c_rt <- c_trad |>
  filter(ACC == 1) |>
  group_by(Subject, Task, Morality, Identity) |>
  summarise(RT = mean(RT), .groups = "drop")

c_acc_w <- c_acc |>
  mutate(key = paste("ACC", Task, Morality, Identity, sep = "_")) |>
  select(Subject, key, ACC) |>
  pivot_wider(names_from = key, values_from = ACC)

c_rt_w <- c_rt |>
  mutate(key = paste("RT", Task, Morality, Identity, sep = "_")) |>
  select(Subject, key, RT) |>
  pivot_wider(names_from = key, values_from = RT)

categ_wide <- inner_join(c_acc_w, c_rt_w, by = "Subject")
save_csv(categ_wide, "MS_categ_behav_wide.csv")

# Collapsed across Task：按作者代码直接从 trials pool 计算，不是简单平均 task means
c_acc_no <- c_trad |>
  group_by(Subject, Morality, Identity) |>
  summarise(ACC = mean(ACC), .groups = "drop")
c_rt_no <- c_trad |>
  filter(ACC == 1) |>
  group_by(Subject, Morality, Identity) |>
  summarise(RT = mean(RT), .groups = "drop")

c_acc_no_w <- c_acc_no |>
  mutate(key = paste("ACC", Morality, Identity, sep = "_")) |>
  select(Subject, key, ACC) |>
  pivot_wider(names_from = key, values_from = ACC)
c_rt_no_w <- c_rt_no |>
  mutate(key = paste("RT", Morality, Identity, sep = "_")) |>
  select(Subject, key, RT) |>
  pivot_wider(names_from = key, values_from = RT)

categ_no_wide <- inner_join(c_acc_no_w, c_rt_no_w, by = "Subject")
save_csv(categ_no_wide, "MS_categ_behav_noTask_wide.csv")

# Cross-task effect scores；RT 方向按作者定义为 comparator - focal，
# 使 positive 值表示 good/self processing advantage。
match_effects <- match_wide |>
  transmute(
    Subject,
    d_goodslf_goodoth = d_Good_Self - d_Good_Other,
    d_goodslf_badslf  = d_Good_Self - d_Bad_Self,
    d_goodoth_badoth  = d_Good_Other - d_Bad_Other,
    d_badslf_badoth   = d_Bad_Self - d_Bad_Other,
    RT_goodslf_goodoth = RT_match_Good_Other - RT_match_Good_Self,
    RT_goodslf_badslf  = RT_match_Bad_Self - RT_match_Good_Self,
    RT_goodoth_badoth  = RT_match_Bad_Other - RT_match_Good_Other,
    RT_badslf_badoth   = RT_match_Bad_Self - RT_match_Bad_Other
  )

categ_effects <- categ_wide |>
  transmute(
    Subject,
    Val_RT_goodslf_goodoth = RT_Val_Good_Other - RT_Val_Good_Self,
    Val_RT_goodslf_badslf  = RT_Val_Bad_Self - RT_Val_Good_Self,
    Val_RT_goodoth_badoth  = RT_Val_Bad_Other - RT_Val_Good_Other,
    Val_RT_badslf_badoth   = RT_Val_Bad_Self - RT_Val_Bad_Other,
    Id_RT_goodslf_goodoth  = RT_Id_Good_Other - RT_Id_Good_Self,
    Id_RT_goodslf_badslf   = RT_Id_Bad_Self - RT_Id_Good_Self,
    Id_RT_goodoth_badoth   = RT_Id_Bad_Other - RT_Id_Good_Other,
    Id_RT_badslf_badoth    = RT_Id_Bad_Self - RT_Id_Bad_Other,
    Val_ACC_goodslf_goodoth = ACC_Val_Good_Self - ACC_Val_Good_Other,
    Val_ACC_goodslf_badslf  = ACC_Val_Good_Self - ACC_Val_Bad_Self,
    Val_ACC_goodoth_badoth  = ACC_Val_Good_Other - ACC_Val_Bad_Other,
    Val_ACC_badslf_badoth   = ACC_Val_Bad_Self - ACC_Val_Bad_Other,
    Id_ACC_goodslf_goodoth  = ACC_Id_Good_Self - ACC_Id_Good_Other,
    Id_ACC_goodslf_badslf   = ACC_Id_Good_Self - ACC_Id_Bad_Self,
    Id_ACC_goodoth_badoth   = ACC_Id_Good_Other - ACC_Id_Bad_Other,
    Id_ACC_badslf_badoth    = ACC_Id_Bad_Self - ACC_Id_Bad_Other
  )

cross_effects <- inner_join(match_effects, categ_effects, by = "Subject")
save_csv(cross_effects, "MS_cross_taskeffect_wide.csv")


# ==============================================================================
# 5. Behavioral statistical analyses
# ==============================================================================

match_rt_long <- m_rt |>
  filter(Match == "match") |>
  transmute(Subject, Val = Morality, Id = Identity, y = RT)
match_d_long <- m_dp |>
  transmute(Subject, Val = Morality, Id = Identity, y = dprime)
categ_rt_long <- c_rt |>
  transmute(Subject, Task, Val = Morality, Id = Identity, y = RT)
categ_acc_long <- c_acc |>
  transmute(Subject, Task, Val = Morality, Id = Identity, y = ACC)
categ_rt_no_long <- c_rt_no |>
  transmute(Subject, Val = Morality, Id = Identity, y = RT)
categ_acc_no_long <- c_acc_no |>
  transmute(Subject, Val = Morality, Id = Identity, y = ACC)

save_csv(match_rt_long, "pilot_match_RT_long.csv")
save_csv(match_d_long, "pilot_match_d_long.csv")
save_csv(categ_rt_long, "pilot_categ_RT_long.csv")
save_csv(categ_acc_long, "pilot_categ_ACC_long.csv")

save_csv(categ_rt_no_long, "pilot_categ_RT_noTask_long.csv")
save_csv(categ_acc_no_long, "pilot_categ_ACC_noTask_long.csv")

tests <- bind_rows(
  planned_tests(match_rt_long,      "match_RT"),
  planned_tests(match_d_long,       "match_d"),
  planned_tests(categ_rt_no_long,   "categ_RT"),
  planned_tests(categ_acc_no_long,  "categ_ACC")
)
save_csv(tests, "pilot_planned_tests.csv")

anova_results <- list(
  match_RT  = run_rm_anova(match_rt_long,  c("Val", "Id"), "match_RT"),
  match_d   = run_rm_anova(match_d_long,   c("Val", "Id"), "match_d"),
  categ_RT  = run_rm_anova(categ_rt_long,  c("Task", "Val", "Id"), "categ_RT"),
  categ_ACC = run_rm_anova(categ_acc_long, c("Task", "Val", "Id"), "categ_ACC")
)

save_csv(bind_rows(lapply(anova_results, `[[`, "frequentist")),
         "pilot_anova_frequentist.csv")
save_csv(bind_rows(lapply(anova_results, `[[`, "omega")),
         "pilot_anova_omega_current_R.csv")
save_csv(bind_rows(lapply(anova_results, `[[`, "bayes_models")),
         "pilot_anova_bayes_models.csv")
save_csv(bind_rows(lapply(anova_results, `[[`, "bayes_inclusion")),
         "pilot_anova_bayes_inclusion.csv")

cors <- cross_correlations(cross_effects)
save_csv(cors, "pilot_cross_task_correlations.csv")


# ==============================================================================
# 6. Pilot ex-Gaussian refit
# ==============================================================================

# 原 preprocessing 只准备 MS_categ_exG.csv；这里使用当前复现流程重新 fit。
# 每个 Subject × Task × Identity × Morality cell 使用 retimes::timefit(iter=1000)。
set.seed(20261006)
ex_groups <- split(
  c_exg,
  interaction(c_exg$Subject, c_exg$Task, c_exg$Identity, c_exg$Morality,
              drop = TRUE)
)

ex_rows <- lapply(ex_groups, function(d) {
  fit_error <- ""
  fit <- if (nrow(d) >= 10) {
    tryCatch(retimes::timefit(d$RT, iter = 1000), error = function(e) {fit_error <<- conditionMessage(e); NULL})
  } else {fit_error <- "Fewer than 10 trials"; NULL}
  par <- if (is.null(fit)) rep(NA_real_, 3) else fit@par
  data.frame(
    Subject = d$Subject[1], Task = d$Task[1],
    Identity = d$Identity[1], Morality = d$Morality[1],
    n = nrow(d), fit_status = if (is.null(fit) || any(!is.finite(c(par,fit@AIC,fit@BIC))) || any(par[2:3] <= 0)) "failed" else "completed", fit_error = fit_error, mu = par[1], sigma = par[2], tau = par[3],
    AIC = if (is.null(fit)) NA_real_ else fit@AIC,
    BIC = if (is.null(fit)) NA_real_ else fit@BIC
  )
})
ex_par <- bind_rows(ex_rows)
save_csv(ex_par, "parameters_all.csv", file.path(OUT, "exgaussian"))

# 作者后续旧流程使用 n > 36 且 8 cells 全部有效的被试
ex_valid <- ex_par |>
  filter(n > 36, fit_status == "completed", complete.cases(mu, sigma, tau)) |>
  group_by(Subject) |>
  filter(n() == 8) |>
  ungroup()
save_csv(ex_valid, "parameters_valid.csv", file.path(OUT, "exgaussian"))

if (nrow(ex_valid)) {
  ex_anova <- bind_rows(lapply(c("mu", "sigma", "tau"), function(param) {
    d <- ex_valid
    for (v in c("Subject", "Task", "Identity", "Morality")) d[[v]] <- factor(d[[v]])
    fit <- afex::aov_ez(
      id = "Subject", dv = param, data = d,
      within = c("Task", "Identity", "Morality"), type = 3,
      anova_table = list(correction = "none", es = "ges")
    )
    z <- as.data.frame(fit$anova_table)
    z$effect <- rownames(z); rownames(z) <- NULL
    z$parameter <- param
    z
  }))
  save_csv(ex_anova, "exgaussian_anova.csv", file.path(OUT, "exgaussian"))
}


# ==============================================================================
# Figures are generated independently by pilot_plot.R

# 8. Reproduction checks + environment record
# ==============================================================================

# 关键 t 值来自论文正文；用于快速判断 preprocessing/statistics 是否走偏。
paper_t <- data.frame(
  outcome = c("match_RT","match_RT","match_d","categ_RT","categ_ACC","categ_ACC"),
  contrast = c("GS_BS","GO_BO","GS_BS","GS_BS","GS_BS","GS_GO"),
  paper_t = c(-5.669, -3.164, 3.260, -3.000, 3.591, 3.430)
)
checks <- inner_join(tests, paper_t, by = c("outcome","contrast")) |>
  mutate(abs_difference = abs(t - paper_t),
         within_target_tolerance = abs_difference < .01) |>
  select(outcome, contrast, paper_t, reproduced_t = t,
         abs_difference, within_target_tolerance)
save_csv(checks, "pilot_key_reproduction_checks.csv")

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
deps <- data.frame(
  Package = PKGS,
  Version = vapply(PKGS, function(p) as.character(packageVersion(p)), character(1))
)
save_csv(deps, "dependencies.csv")

message("Pilot reproduction complete. Results: ", OUT)

save_csv(data.frame(study=STUDY,input=basename(c(MATCH_FILE,CATEG_FILE))),"input_provenance.csv")
writeLines("Analysis finished; inspect statistical checks and any ex-Gaussian failures before interpretation.",file.path(OUT,"completed.txt"))
