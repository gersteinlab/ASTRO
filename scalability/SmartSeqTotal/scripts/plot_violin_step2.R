# convert_to_violin_format.R
library(dplyr)

# 读取两个 CSV
#astro_data <- read.csv("/vast/palmer/scratch/jun_lu/dz287/colab/yh742/smart/benchmarking/astro_metrics_summary.csv")
#slide_data <- read.csv("/vast/palmer/scratch/jun_lu/dz287/colab/yh742/smart/benchmarking/smart_metrics_summary.csv")  # smart 数据

#astro_data <- read.csv("/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/benchmarking/boxpic927-1/metrics_summary_astro.csv")
#slide_data <- read.csv("/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/benchmarking/boxpic927-2/metrics_summary.csv")  # smart 数据

astro_data <- read.csv("/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/benchmarking/astro_metrics_summary1012.csv")
slide_data <- read.csv("/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/outputvisium/spaceranger_metrics_summary1012.csv")  # smart 数据

# 合并并转换格式
df_sil <- bind_rows(
  data.frame(down = astro_data$down_prob, data_type = "astro", value = astro_data$silhouette),
  data.frame(down = slide_data$down_prob, data_type = "spaceranger", value = slide_data$silhouette)
)

df_cal <- bind_rows(
  data.frame(down = astro_data$down_prob, data_type = "astro", value = astro_data$CH),
  data.frame(down = slide_data$down_prob, data_type = "spaceranger", value = slide_data$CH)
)

df_dav <- bind_rows(
  data.frame(down = astro_data$down_prob, data_type = "astro", value = astro_data$DB),
  data.frame(down = slide_data$down_prob, data_type = "spaceranger", value = slide_data$DB)
)

# 写出 TSV
write.table(df_sil, "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/benchmarking/df_silhouette1012.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_cal, "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/benchmarking/df_calinski1012.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_dav, "/vast/palmer/scratch/jun_lu/dz287/colab/yh742/spaceranger/cscc/benchmarking/df_davies1012.tsv", sep = "\t", row.names = FALSE, quote = FALSE)