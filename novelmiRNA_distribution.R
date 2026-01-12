# ===================
# INSTALL AND LOAD LIBRARIES
# ===================

library(dplyr)
library(tidyr)
library(ggplot2)

# ===================
# LOAD DATA
# ===================
project <- "SRA_bmarrow"
# FORWARD STRAND
fwd <- read.delim(file.path(paste0("C:/Users/javie/OneDrive/Escritorio/Project in Molecular Life Sciences/", "Analysis/", "AGGGCTGAACCTAGAATT_greps_", project, ".txt")),
                  header = FALSE)
colnames(fwd) <- c('counts', 'seq', 'sample')
# Sum reads from the same sample 
fwd_sum <- fwd %>%
  group_by(sample) %>%
  summarise(total = sum(counts))
fwd_sum$strand <- 'fwd'

# REV STRAND
rev <- read.delim(file.path(paste0("C:/Users/javie/OneDrive/Escritorio/Project in Molecular Life Sciences/", "Analysis/", "AATTCTAGGTTCAGCCCT_greps_", project, ".txt")),
                  header = FALSE)
colnames(rev) <- c('counts', 'seq', 'sample')
rev_sum <- rev %>%
  group_by(sample) %>%
  summarise(total = sum(counts))
rev_sum$strand <- 'rev'

# ===================
# PLOT THE DISTRIBUTIONS
# ===================

df <- rbind(fwd_sum, rev_sum)
if (project == "SRA_bmarrow"){
  project_long = "SRA bone marrow"
} else {
  project_long = project
}

count_distr <- df %>%
  ggplot(aes(x = total, fill = strand)) +
  geom_histogram(position = 'identity',  alpha = 0.4, bins = 30) +
  theme_minimal() +
  labs(
    title = paste0('Count distribution in ', project_long,' data of potential novel miRNA'),
    x = 'Counts',
    y = 'Frequency'
  ) +
  scale_x_continuous(limits = c(0, max(df$total) + 1))

count_distr
ggsave(paste0("novel_miRNA_count_distr_", project, ".png"), plot = count_distr,
       dpi = 600, device = "png", width = 8, height = 8)
