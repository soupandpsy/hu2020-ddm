# Hu et al. (2020): confirm behavioral analysis.
# Figures are generated separately by confirm_plot.R.
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

STUDY <- "confirm"
DATA <- normalizePath(arg("--data-dir",file.path(ROOT,"Data","raw",STUDY)),mustWork=TRUE)
PKGS <- c("Rcpp","dplyr","tidyr","ggplot2","afex","BayesFactor","effectsize")
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

MATCH_FILE <- find_input("MS_rep_matchingTask_raw.csv")
CATEG_FILE <- find_input("MS_rep_categTask_raw.csv")

m_raw <- read.csv(MATCH_FILE, stringsAsFactors = FALSE) |>
  mutate(
    Match = factor(Match, levels = c("match", "mismatch")),
    Morality = factor(Morality, levels = c("Good", "Bad")),
    Identity = factor(Identity, levels = c("Self", "Other"))
  )

c_raw <- read.csv(CATEG_FILE, stringsAsFactors = FALSE) |>
  mutate(
    Task = factor(Task, levels = c("Id", "Val")),
    Morality = factor(Morality, levels = c("Good", "Bad")),
    Identity = factor(Identity, levels = c("Self", "Other"))
  )


# ==============================================================================
# 3. Confirmatory preprocessing
# ==============================================================================

# Rule 1: 每个 condition 的 trial number 必须完整
excl1_m <- m_raw |>
  mutate(.acc = ifelse(ACC == 1, 1, 0)) |>
  count(Subject, Match, Identity, Morality, name = "N") |>
  filter(N != 75) |>
  distinct(Subject) |>
  pull(Subject)

excl1_c <- c_raw |>
  mutate(.acc = ifelse(ACC == 1, 1, 0)) |>
  count(Subject, Task, Identity, Morality, name = "N") |>
  filter(N != 90) |>
  distinct(Subject) |>
  pull(Subject)

# Rule 2: overall ACC < .50；作者代码把所有非 1 response 计为错误
subject_acc <- function(dat) {
  dat |>
    mutate(.acc = ifelse(ACC == 1, 1, 0)) |>
    group_by(Subject) |>
    summarise(ACC = mean(.acc), .groups = "drop")
}
excl2_m <- subject_acc(m_raw) |> filter(ACC < .5) |> pull(Subject)
excl2_c <- subject_acc(c_raw) |> filter(ACC < .5) |> pull(Subject)

# Rule 3: 任一 experimental condition ACC == 0
excl3_m <- m_raw |>
  mutate(.acc = ifelse(ACC == 1, 1, 0)) |>
  group_by(Subject, Match, Identity, Morality) |>
  summarise(ACC = mean(.acc), .groups = "drop") |>
  filter(ACC == 0) |>
  distinct(Subject) |>
  pull(Subject)

excl3_c <- c_raw |>
  mutate(.acc = ifelse(ACC == 1, 1, 0)) |>
  group_by(Subject, Task, Identity, Morality) |>
  summarise(ACC = mean(.acc), .groups = "drop") |>
  filter(ACC == 0) |>
  distinct(Subject) |>
  pull(Subject)

excl_m <- unique(c(excl1_m, excl2_m, excl3_m))
excl_c <- unique(c(excl1_c, excl2_c, excl3_c))

warn_expected("Confirm matching excluded N", length(excl_m), 2)
warn_expected("Confirm categorization excluded N", length(excl_c), 3)

# Traditional analysis：
# 1) 所有非 ACC==1 统一作为 incorrect
# 2) 只删除 correct & RT <= .2 s
# 3) 行为分析统一把 RT 转换到 ms
m_trad <- m_raw |>
  filter(!Subject %in% excl_m) |>
  mutate(ACC = ifelse(ACC == 1, 1, 0)) |>
  filter(!(RT <= .2 & ACC == 1)) |>
  mutate(RT = RT * 1000)

c_trad <- c_raw |>
  filter(!Subject %in% excl_c) |>
  mutate(ACC = ifelse(ACC == 1, 1, 0)) |>
  filter(!(RT <= .2 & ACC == 1)) |>
  mutate(RT = RT * 1000)

warn_expected("Confirm matching valid N", n_distinct(m_trad$Subject), 42)
warn_expected("Confirm categorization valid N", n_distinct(c_trad$Subject), 41)

make_exclusion_rows <- function(task, rule, ids) {
  data.frame(
    task = rep(task, length(ids)),
    rule = rep(rule, length(ids)),
    Subject = ids
  )
}

exclusions <- bind_rows(
  make_exclusion_rows("matching", "wrong_trial_number", excl1_m),
  make_exclusion_rows("matching", "overall_ACC_lt_.5", excl2_m),
  make_exclusion_rows("matching", "zero_ACC_condition", excl3_m),
  make_exclusion_rows("categorization", "wrong_trial_number", excl1_c),
  make_exclusion_rows("categorization", "overall_ACC_lt_.5", excl2_c),
  make_exclusion_rows("categorization", "zero_ACC_condition", excl3_c)
)
save_csv(exclusions, "confirm_exclusions.csv")

sample_flow <- data.frame(
  task = c("matching", "categorization"),
  raw_N = c(n_distinct(m_raw$Subject), n_distinct(c_raw$Subject)),
  excluded_N = c(length(excl_m), length(excl_c)),
  final_N = c(n_distinct(m_trad$Subject), n_distinct(c_trad$Subject))
)
save_csv(sample_flow, "confirm_sample_flow.csv")
save_csv(m_trad, "confirm_matching_clean_traditional.csv")
save_csv(c_trad, "confirm_categorization_clean_traditional.csv")


# ==============================================================================
# 4. HDDM inputs (seconds; no <200 ms exclusion)
# ==============================================================================

m_hddm_match <- m_raw |>
  filter(!Subject %in% excl_m, ACC %in% c(0, 1), Match == "match") |>
  transmute(subj_idx = Subject, match = Match, val = Morality,
            id = Identity, rt = RT, response = ACC)

m_hddm_nonmatch <- m_raw |>
  filter(!Subject %in% excl_m, ACC %in% c(0, 1), Match == "mismatch") |>
  transmute(subj_idx = Subject, val = Morality, id = Identity,
            rt = RT, response = ACC)

m_hddm_stim <- m_raw |>
  filter(!Subject %in% excl_m, ACC %in% c(0, 1)) |>
  mutate(
    response = ifelse((Match == "match" & ACC == 1) |
                      (Match == "mismatch" & ACC == 0), 1, 0),
    stim = ifelse(Match == "match", 1, 0)
  ) |>
  transmute(subj_idx = Subject, match = Match, val = Morality,
            id = Identity, stim, response, rt = RT)

c_hddm_val <- c_raw |>
  filter(!Subject %in% excl_c, ACC %in% c(0, 1), Task == "Val") |>
  mutate(
    response = ifelse((Morality == "Good" & ACC == 1) |
                      (Morality == "Bad" & ACC == 0), 1, 0),
    stim = ifelse(Shape %in% c("moralSelf", "moralOther"), 1, 0)
  ) |>
  transmute(subj_idx = Subject, val = Morality, id = Identity,
            stim, response, rt = RT)

c_hddm_id <- c_raw |>
  filter(!Subject %in% excl_c, ACC %in% c(0, 1), Task == "Id") |>
  mutate(
    response = ifelse((Identity == "Self" & ACC == 1) |
                      (Identity == "Other" & ACC == 0), 1, 0),
    stim = ifelse(Shape %in% c("moralSelf", "immoralSelf"), 1, 0)
  ) |>
  transmute(subj_idx = Subject, val = Morality, id = Identity,
            stim, response, rt = RT)

save_csv(m_hddm_match,    "MS_match_hddm.csv", file.path(OUT, "hddm_input"))
save_csv(m_hddm_nonmatch, "MS_mismatch_hddm.csv", file.path(OUT, "hddm_input"))
save_csv(m_hddm_stim,     "MS_match_hddm_stim.csv", file.path(OUT, "hddm_input"))
save_csv(c_hddm_val,      "MS_categ_val_hddm_stim.csv", file.path(OUT, "hddm_input"))
save_csv(c_hddm_id,       "MS_categ_id_hddm_stim.csv", file.path(OUT, "hddm_input"))


# ==============================================================================
# 5. Behavioral variables
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
  group_by(Subject, Age, Sex, Morality, Identity) |>
  summarise(
    hit  = sum(Match == "match"    & ACC == 1),
    miss = sum(Match == "match"    & ACC == 0),
    CR   = sum(Match == "mismatch" & ACC == 1),
    FA   = sum(Match == "mismatch" & ACC == 0),
    .groups = "drop"
  ) |>
  rowwise() |>
  mutate(dprime = dprime_value(hit, miss, FA, CR)) |>
  ungroup()

m_dp_w <- m_dp |>
  mutate(key = paste("d", Morality, Identity, sep = "_")) |>
  select(Subject, Age, Sex, key, dprime) |>
  pivot_wider(names_from = key, values_from = dprime)

# Matching RT
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

# Collapsed across Task：直接 pool trials，与作者代码一致
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

# Cross-task effects
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
# 6. Behavioral statistical analyses
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

save_csv(match_rt_long, "confirm_match_RT_long.csv")
save_csv(match_d_long, "confirm_match_d_long.csv")
save_csv(categ_rt_long, "confirm_categ_RT_long.csv")
save_csv(categ_acc_long, "confirm_categ_ACC_long.csv")

save_csv(categ_rt_no_long, "confirm_categ_RT_noTask_long.csv")
save_csv(categ_acc_no_long, "confirm_categ_ACC_noTask_long.csv")

tests <- bind_rows(
  planned_tests(match_rt_long,      "match_RT"),
  planned_tests(match_d_long,       "match_d"),
  planned_tests(categ_rt_no_long,   "categ_RT"),
  planned_tests(categ_acc_no_long,  "categ_ACC")
)
save_csv(tests, "confirm_planned_tests.csv")

anova_results <- list(
  match_RT  = run_rm_anova(match_rt_long,  c("Val", "Id"), "match_RT"),
  match_d   = run_rm_anova(match_d_long,   c("Val", "Id"), "match_d"),
  categ_RT  = run_rm_anova(categ_rt_long,  c("Task", "Val", "Id"), "categ_RT"),
  categ_ACC = run_rm_anova(categ_acc_long, c("Task", "Val", "Id"), "categ_ACC")
)

save_csv(bind_rows(lapply(anova_results, `[[`, "frequentist")),
         "confirm_anova_frequentist.csv")
save_csv(bind_rows(lapply(anova_results, `[[`, "omega")),
         "confirm_anova_omega_current_R.csv")
save_csv(bind_rows(lapply(anova_results, `[[`, "bayes_models")),
         "confirm_anova_bayes_models.csv")
save_csv(bind_rows(lapply(anova_results, `[[`, "bayes_inclusion")),
         "confirm_anova_bayes_inclusion.csv")

cors <- cross_correlations(cross_effects)
save_csv(cors, "confirm_cross_task_correlations.csv")


# ==============================================================================
# Figures are generated independently by confirm_plot.R

# 8. Reproduction checks + environment record
# ==============================================================================

# 关键 t 值来自论文正文。论文中 p 值存在单/双侧报告混用，因此这里主要核 t。
paper_t <- data.frame(
  outcome = c("match_RT","match_RT","match_RT",
              "match_d","match_d",
              "categ_RT","categ_RT",
              "categ_ACC","categ_ACC"),
  contrast = c("GS_BS","GO_BO","GS_GO",
               "GS_BS","GS_GO",
               "GS_BS","GS_GO",
               "GS_BS","GS_GO"),
  paper_t = c(-8.42, -2.37, -3.34,
               4.30, 2.67,
              -3.18, -4.37,
               3.06, 3.77)
)
checks <- inner_join(tests, paper_t, by = c("outcome","contrast")) |>
  mutate(abs_difference = abs(t - paper_t),
         within_target_tolerance = abs_difference < .02) |>
  select(outcome, contrast, paper_t, reproduced_t = t,
         abs_difference, within_target_tolerance)
save_csv(checks, "confirm_key_reproduction_checks.csv")

# 记录论文正文一个已知可疑点：
# matching RT Valence effect 正文写 F(1,40)=57.88，
# 但 N=42 的完整 within-subject design 应为 denominator df=41，
# 原始数据重算也应得到约 F=57.88, df=41。
notes <- c(
  "Behavioral RT outputs in this script are standardized to milliseconds.",
  "Confirm matching final N should be 42; categorization final N should be 41.",
  "Paper reports matching RT Valence F(1,40)=57.88, but N=42 implies df=41; treat the printed df=40 as a likely typo if reproduced F is ~57.88.",
  "Paper confirm categorization Good-self vs Bad-self accuracy prints negative d despite positive t and positive CI; preserve computed dz sign rather than copying the likely typo.",
  "HDDM input files are exported, but HDDM fitting is not performed in this R script."
)
writeLines(notes, file.path(OUT, "reproduction_notes.txt"))

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
deps <- data.frame(
  Package = PKGS,
  Version = vapply(PKGS, function(p) as.character(packageVersion(p)), character(1))
)
save_csv(deps, "dependencies.csv")

message("Confirmatory reproduction complete. Results: ", OUT)

save_csv(data.frame(study=STUDY,input=basename(c(MATCH_FILE,CATEG_FILE))),"input_provenance.csv")
writeLines("Analysis finished; inspect statistical checks and any ex-Gaussian failures before interpretation.",file.path(OUT,"completed.txt"))
