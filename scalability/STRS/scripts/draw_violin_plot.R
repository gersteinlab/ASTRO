#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(fs)
})


base_dir <- "./"
newpipe  <- file.path(base_dir, "result/")
outdir   <- file.path(base_dir, "result/")
dir_create(outdir)


files <- c(
  astro = file.path(newpipe, "astrodownsample.txt"),
  strs  = file.path(newpipe, "theirdownsample.txt")
)

down_value <- 0.5
read_one <- function(path, data_type){
  stopifnot(file_exists(path))
  df <- read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "\"", comment.char = "")
  need <- c("rep","resolution","silhouette","CH","DB")
  miss <- setdiff(need, names(df))
  if(length(miss) > 0){
    stop(sprintf("缺少列：%s（文件：%s）", paste(miss, collapse = ", "), path))
  }
  df$data_type <- data_type
  df
}

d_ast <- read_one(files[["astro"]], "astro")
d_str <- read_one(files[["strs"]],  "strs")
d_all <- bind_rows(d_ast, d_str)

make_metric <- function(df, metric_col, out_tsv){
  out <- df %>%
    transmute(
      down = down_value,
      data_type = .data[["data_type"]],
      value = as.numeric(.data[[metric_col]])
    )
  write_tsv(out, out_tsv)
  message("Wrote: ", out_tsv, "  (n=", nrow(out), ")")
}

make_metric(d_all, "silhouette",
            file.path(outdir, "df_silhouette1014.tsv"))
make_metric(d_all, "CH",
            file.path(outdir, "df_calinski1014.tsv"))
make_metric(d_all, "DB",
            file.path(outdir, "df_davies1014.tsv"))


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)   # stat_compare_means, theme_pubclean, geom_bracket
})

args <- commandArgs(trailingOnly = TRUE)
outname <- ifelse(length(args) >= 1, args[1], "/vast/palmer/scratch/jun_lu/dz287/colab/trySTRS/benchmarking/violin_scores.png")

read_metric <- function(path, metric_name){
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    rename(down = down, data_type = data_type, value = value) %>%
    mutate(Metric = metric_name)
}

df_sil <- read_metric("result/df_silhouette1014.tsv", "Silhouette")
df_cal <- read_metric("result/df_calinski1014.tsv", "Calinski-Harabasz")
df_dav <- read_metric("result/df_davies1014.tsv", "Davies-Bouldin")

test.data <- bind_rows(df_sil, df_dav, df_cal) %>%
  filter(as.numeric(down) == 0.5 | down == "0.5") %>%
  mutate(
    Score = value,
    Sample = ifelse(data_type == "astro", "ASTRO", "STRS"),
    Sample = factor(Sample, levels = c("ASTRO", "STRS")),
    group = paste0("0.5_", data_type),
    sample = as.character(as.numeric(factor(group, levels = c("0.5_astro","0.5_STRS")))),
    Metric = factor(Metric, levels = c("Silhouette", "Calinski-Harabasz", "Davies-Bouldin"))
  ) %>%
  group_by(Metric, data_type) %>%
  slice_head(n = 40) %>%
  ungroup()

ypos_df <- test.data %>%
  group_by(Metric) %>%
  summarise(ypos = min(Score, na.rm = TRUE) - 0.05, .groups = "drop")

#dp <- ggplot(test.data, aes(x = sample, y = Score, fill = Sample)) +
dp <- ggplot(test.data, aes(x = Sample, y = Score, fill = Sample)) +
  geom_violin(trim = FALSE,scale = "width", width = 0.8) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(x = NULL, y = "Score") +
  scale_fill_brewer(palette = "Set3", na.value = "#222A35") +
  theme_pubclean() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ Metric, nrow = 1, scales = "free_y") +
  stat_compare_means(
    #comparisons = list(c("1", "2")),
    comparisons = list(c("ASTRO", "STRS")),
    method = "wilcox.test",
    method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    label = "p.format"
  )

ggsave(outname, dp, units = "px", width = 2400, height = 1200)

file.remove("result/df_silhouette1014.tsv", "result/df_calinski1014.tsv", "result/df_davies1014.tsv")

