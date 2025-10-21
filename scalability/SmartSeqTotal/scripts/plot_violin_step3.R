#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)   # stat_compare_means, theme_pubclean, geom_bracket
})

args <- commandArgs(trailingOnly = TRUE)
outname <- ifelse(length(args) >= 1, args[1], "result/violin_scores1012.png")

read_metric <- function(path, metric_name){
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    rename(down = down, data_type = data_type, value = value) %>%
    mutate(Metric = metric_name)
}

df_sil <- read_metric("cache/df_silhouette1012.tsv", "Silhouette")
df_cal <- read_metric("cache/df_calinski1012.tsv", "Calinski-Harabasz")
df_dav <- read_metric("cache/df_davies1012.tsv", "Davies-Bouldin")

test.data <- bind_rows(df_sil, df_dav, df_cal) %>%
  filter(as.numeric(down) == 0.5 | down == "0.5") %>%
  mutate(
    Score = value,
    Sample = ifelse(data_type == "astro", "ASTRO", "spaceranger"),
    Sample = factor(Sample, levels = c("ASTRO", "spaceranger")),
    group = paste0("0.5_", data_type),
    sample = as.character(as.numeric(factor(group, levels = c("0.5_astro","0.5_spaceranger")))),
    Metric = factor(Metric, levels = c("Silhouette", "Calinski-Harabasz", "Davies-Bouldin"))
  ) %>%
  group_by(Metric, data_type) %>%
  slice_head(n = 20) %>%
  ungroup()

ypos_df <- test.data %>%
  group_by(Metric) %>%
  summarise(ypos = min(Score, na.rm = TRUE) - 0.05, .groups = "drop")

dp <- ggplot(test.data, aes(x = sample, y = Score, fill = Sample)) +
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
    comparisons = list(c("1", "2")),
    method = "wilcox.test",
    method.args = list(exact = FALSE, correct = TRUE),
    size = 3,
    label = "p.format"
  )

ggsave(outname, dp, units = "px", width = 2400, height = 1200)
