################################################################################
### AA-GP: process and export summary of null distribution for each pathway
### Sarah Turner-Hissong
### 12 June 2020
################################################################################
# This script summarizes random gene groups for each pathway tested
# There are 1000 gene groups per pathway
# This serves as an empirical null distribution for each pathway to test
#   the significance of the MultiBLUP model
# output includes the following for each gene group:
#   the number of unique genes
#   the number of unique snps
#   the likelihood ratio test statistic from MultiBLUP
#   the proportion of heritability explained

library(tidyverse)
library(data.table)

################################################################################
# Step 1: summarise random gene groups

# use command line arguments
args = commandArgs(trailingOnly = TRUE)

# filepath
path <- paste0("data/processed/random_sets_pathways/", args[1])

# list files for random gene groups
random_files <- data.frame(
  filename = list.files(path = path,
                        pattern = ".txt$",
                        full.names = TRUE,
                        recursive = TRUE)
)

random_files <- data.frame(filename =
                         system(paste("find", path, "| grep .txt", sep = " "),
                         intern = TRUE)
                       )

# summarize number of snps and genes in each random group

random_summary <- random_files %>%
  mutate(filename = as.character(filename),
         # pathway id
         pathway = word(filename, start = 4, sep = "/"),
         # random group id (1-1000)
         random_group = word(filename, start = 5, sep = "/"),
         # cleaning up a bit.. this is kind of messy
         random_group = gsub(".txt", "", random_group),
         random_group = gsub("null_", "c_", random_group),
         # read in file contents into a list
         file_contents = map(filename,
                             ~ fread(., drop = 1))) %>%
#                             ~ read.table(.))) %>%
  # unnest list into a data frame
  unnest() %>%
  group_by(pathway, random_group) %>%
  # count number of unique genes and snps in each random gene group
  summarize(n_genes = n_distinct(ensembl_gene_id),
            n_snps = n_distinct(snp_id))

head(random_summary)

################################################################################
# Step 2: compute likelihood ratio test statistic

# import log likelihood for GBLUP model
# this is used to compute the likelihood ratio statistic
load("reports/gblup.RData")

# filepath to model output from REML
path <- paste0("models/null_h2_pathways/", args[1])

## list files for reml output of each random gene group

# list files using 'find' OS command
# this option is *much* faster than 'list.files'
lr_files <- data.frame(filename =
                         system(paste("find", path, "| grep .reml", sep = " "),
                         intern = TRUE)
                       )

lr_null <- lr_files %>%
  mutate(filename = as.character(filename),
         # pathway id
         pathway = word(filename, start = 3, sep = "/"),
         # random group id (1-1000)
         random_group = word(filename, start = 4, sep = "/"),
         # traid id
         trait = word(filename, start = 5, sep = "/|[.]"),
         # read in file contents into a list
         file_contents = map(filename,
                             ~ fread(.,
                                     header = FALSE,
                                     skip = 10,
                                     nrows = 1))) %>%
  # unnest list into a data frame
  unnest() %>%
#   filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(null_llik = as.numeric(as.character(V2))) %>%
  left_join(., gblup_llik, by = "trait") %>%
  select(trait, pathway, random_group, null_llik, gblup_llik) %>%
  mutate(LR = 2 * (null_llik - gblup_llik))

head(lr_null)

################################################################################
## Step 3: summarize proportion of variance explained

## list files for proportion of variance explained
# uses same path as LR output

h2_files <- data.frame(filename =
                         system(paste("find", path, "| grep .share", sep = " "),
                         intern = TRUE)
                       )


# h2_files <- data.frame(filename =
#                          list.files(path = path,
#                                     pattern = ".share$",
#                                     full.names = TRUE,
#                                     recursive = TRUE)
# )

h2_null <- h2_files %>%
  mutate(filename = as.character(filename),
         # pathway id
         pathway = word(filename, start = 3, sep = "/"),
         # random group id (1-1000)
         random_group = word(filename, start = 4, sep = "/"),
         # traid id
         trait = word(filename, start = 5, sep = "/|[.]"),
         # read in file contents into a list
         file_contents = map(filename,
                             ~ fread(., header = TRUE))) %>%
  # unnest list into a data frame
  unnest() %>%
  select(trait, pathway, random_group, Component, Share, SD) %>%
  left_join(., random_summary, by = c("pathway", "random_group"))

head(h2_null)

################################################################################
## Step 4: merge and export results!

output <- paste0("reports/null_dist_results_pathways/", args[1], "_null.csv")

null_dist <- left_join(lr_null, h2_null,
                       by = c("trait", "pathway", "random_group")) %>%
  write_csv(., output)
