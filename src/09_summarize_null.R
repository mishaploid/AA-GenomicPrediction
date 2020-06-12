################################################################################
### AA-GP: process and export summary of null distribution
### Sarah Turner-Hissong
### 07 March 2019

library(tidyverse)
library(data.table)
# library(quantreg)

# summarize random gene groups (number of SNPs & number of genes)

feature_sizes <- read_csv("reports/random_subsets_summary.csv")
head(feature_sizes)

# import log likelihood for GBLUP model
load("reports/gblup.RData")

## list files for reml output of each random gene group
files <- list.files(path = "models/null_h2", pattern = ".reml$",
                    full.names = TRUE, recursive = TRUE)

lr_null <- rbindlist(lapply(files, read.table, header = TRUE),
                      idcol = "filename", use.names = TRUE) %>%
  mutate(filename = files[filename]) %>%
  separate(filename, sep = "/|[.]",
           into = c("dir", "source", "pathway", "trait", "metric", "model"), remove = FALSE) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(null_llik = as.numeric(as.character(X2))) %>%
  left_join(., gblup_llik, by = "trait") %>%
  select(trait, pathway, Num_Kinships, null_llik, gblup_llik) %>%
  left_join(., feature_sizes, by = "pathway") %>%
  mutate(group_size = cut(size, breaks = seq(0, 50000, by = 10000),
                          include.lowest = TRUE, dig.lab = 10),
         LR = 2 * (null_llik - gblup_llik))

# lr_null <- rbindlist(lapply(files, read.table, header = TRUE),
#                       idcol = "filename", use.names = TRUE) %>%
#   mutate(filename = files[filename]) %>%
#   separate(filename, sep = "/|[.]",
#            into = c("dir", "source", "pathway", "trait", "metric", "model"), remove = FALSE) %>%
#   filter(Num_Kinships %in% "Alt_Likelihood") %>%
#   # mutate(null_llik = as.numeric(as.character(X2))) %>%
#   # left_join(., gblup_llik, by = "trait") %>%
#   # select(trait, pathway, Num_Kinships, null_llik, gblup_llik) %>%
#   rename("LRT_Stat" = X2) %>%
#   select(trait, pathway, LRT_Stat) %>%
#   left_join(., feature_sizes, by = "pathway") %>%
#   mutate(group_size = cut(size, breaks = seq(0, 50000, by = 10000),
#                           include.lowest = TRUE, dig.lab = 10))

head(lr_null)

## read in variance component estimates
files <- list.files(path = "models/null_h2", pattern = ".share$",
                    full.names = TRUE, recursive = TRUE)

h2_null <- rbindlist(lapply(files, read.table, header = TRUE),
                      idcol = "filename", use.names = TRUE) %>%
  mutate(filename = files[filename]) %>%
  separate(filename, sep = "/|[.]", into = c("dir", "source", "pathway", "trait", "metric", "model"), remove = FALSE) %>%
  left_join(., feature_sizes, by = "pathway") %>%
  select(trait, pathway, size, n_genes, Component, Share, SD)

head(h2_null)

## merge!
null_dist <- left_join(lr_null, h2_null, by = c("trait", "pathway")) %>%
  write_csv(., "reports/lr_null_results.csv")
