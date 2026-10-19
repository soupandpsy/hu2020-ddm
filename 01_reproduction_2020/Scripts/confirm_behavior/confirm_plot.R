
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
INPUT <- arg("--analysis-dir")
if(is.null(INPUT)) stop("--analysis-dir is required")
INPUT <- normalizePath(INPUT,mustWork=TRUE)
if(!file.exists(file.path(INPUT,"completed.txt"))) stop("Analysis is incomplete")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2)})
OUT <- new_output(arg("--output-dir",file.path(INPUT,paste0("plots_",stamp()))))
dir.create(file.path(OUT,"figures"))
match_rt_long <- read.csv(file.path(INPUT,"confirm_match_RT_long.csv"))
match_d_long <- read.csv(file.path(INPUT,"confirm_match_d_long.csv"))
categ_rt_long <- read.csv(file.path(INPUT,"confirm_categ_RT_long.csv"))
categ_acc_long <- read.csv(file.path(INPUT,"confirm_categ_ACC_long.csv"))
categ_rt_no_long <- read.csv(file.path(INPUT,"confirm_categ_RT_noTask_long.csv"))
categ_acc_no_long <- read.csv(file.path(INPUT,"confirm_categ_ACC_noTask_long.csv"))
plot_two_dv <- function(rt_long, other_long, other_label, title, basename) {
  pdat <- bind_rows(
    rt_long |> transmute(Subject, Val, Id, DV = "RT (ms)", value = y),
    other_long |> transmute(Subject, Val, Id, DV = other_label, value = y)
  ) |>
    mutate(
      # 固定 facet 顺序：RT 始终在左，d-prime / Accuracy 在右
      DV = factor(DV, levels = c("RT (ms)", other_label)),
      Val = factor(Val, levels = c("Good", "Bad")),
      Id = factor(Id, levels = c("Self", "Other")),
      x = as.numeric(Val) + ifelse(Id == "Self", -0.12, 0.12)
    )

  sm <- pdat |>
    group_by(DV, Val, Id) |>
    summarise(mean = mean(value), sd = sd(value), n = n(),
              se = sd / sqrt(n), .groups = "drop") |>
    mutate(x = as.numeric(Val) + ifelse(Id == "Self", -0.12, 0.12))

  p <- ggplot() +
    geom_line(data = pdat,
              aes(x = x, y = value, group = Subject),
              linewidth = .35, alpha = .08) +
    geom_point(data = pdat, aes(x = x, y = value),
               size = 1.2, alpha = .10) +
    geom_line(data = sm, aes(x = x, y = mean, group = Id, linetype = Id),
              linewidth = .9) +
    geom_point(data = sm, aes(x = x, y = mean, shape = Id), size = 3) +
    geom_errorbar(data = sm,
                  aes(x = x, ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                  width = .04, linewidth = .7) +
    scale_x_continuous(breaks = c(1, 2), labels = c("Good", "Bad")) +
    facet_wrap(~DV, scales = "free_y", nrow = 1) +
    labs(title = title, x = NULL, y = NULL, shape = NULL, linetype = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(OUT, "figures", paste0(basename, ".pdf")), p,
         width = 9, height = 5.5)
  ggsave(file.path(OUT, "figures", paste0(basename, ".png")), p,
         width = 9, height = 5.5, dpi = 180)
  invisible(p)
}
plot_two_dv(
  match_rt_long, match_d_long, "d-prime",
  "A. Matching task", "confirm_matching_behavior"
)

# Figure 3B：Identity-based categorization
plot_two_dv(
  categ_rt_long |> filter(Task == "Id"),
  categ_acc_long |> filter(Task == "Id"),
  "Accuracy",
  "B. Identity-based categorization task",
  "confirm_identity_behavior"
)

# Figure 3C：Valence-based categorization
# 原作者 MS_rep_preproc.r 在 group mean/error bars 中误用了 df.C1_Id.sum_p；
# 这里使用真正的 Valence-task participant means。
plot_two_dv(
  categ_rt_long |> filter(Task == "Val"),
  categ_acc_long |> filter(Task == "Val"),
  "Accuracy",
  "C. Valence-based categorization task",
  "confirm_valence_behavior"
)
writeLines(INPUT,file.path(OUT,"analysis_source.txt"))
message("Figures: ",OUT)
