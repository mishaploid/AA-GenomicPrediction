################################################################################
## AA-GP: summarize GBLUP results
## Sarah Turner-Hissong
## Updated 18 March 2019

library(tidyverse)
library(broom)

# use command line args
args <- commandArgs(trailingOnly = TRUE)

# phenotypes
pheno <- read.table(paste0(args[1]), header = TRUE) %>%
  gather(key = "trait",
         value = "phenotype",
         -FID, -IID)

# heritability
gblup_h2 <- tibble(file = list.files(path = "models/gblup_h2",
                                     pattern = ".reml$",
                                     full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "trait", "method2", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE, skip = 13)) %>%
  unnest(data) %>%
  filter(Component %in% "Her_ALL") %>%
  select(trait, Heritability, Her_SD)

head(gblup_h2)

# log likelihood
gblup_llik <- tibble(file = list.files(path = "models/gblup_h2",
                                     pattern = ".reml$",
                                     full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "trait", "metric", "fill", "fill2"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(gblup_llik = as.numeric(as.character(X1))) %>%
  select(trait, gblup_llik)

head(gblup_llik)

# individuals used for test set
cross_val <- tibble(filename = list.files(path = "data/processed/cross_validation",
                    pattern = ".test", full.names = TRUE)) %>%
             mutate(cvnum = word(filename, start = 4, end = 4, sep = "/|[.]"),
                    fold = word(filename, start = 5, end = 5, sep = "/|[.]"),
                    fold = gsub("test", "", fold)) %>%
             mutate(file_contents = map(filename, ~read.table(., header = FALSE))) %>%
             unnest(cols = c("file_contents")) %>%
             rename("FID" = "V1",
                    "IID" = "V2")

# predictive accuracy
gblup_pa <- tibble(file = list.files(path = "models/gblup",
                                     pattern = ".pred$",
                                     full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "trait", "cvnum", "fold", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  left_join(pheno, ., by = c("FID" = "ID1", "IID" = "ID2", "trait")) %>%
  inner_join(., cross_val)
  # left_join(., training_coef, by = c("trait", "cvnum", "fold"))

# training_coef <- tibble(file = list.files(path = "models/gblup",
#                                           pattern = ".coeff$",
#                                           full.names = TRUE)) %>%
#   separate(file, sep = "/|[.]", into = c("source", "method", "trait", "cvnum", "fold", "fill"), remove = FALSE) %>%
#   mutate(data = lapply(file, read.table, header = TRUE)) %>%
#   unnest(data) %>%
#   filter(Component %in% "Intercept") %>%
#   select(trait, cvnum, fold, Effect)
#
# head(training_coef)
#
# gblup_pa <- tibble(file = list.files(path = "models/gblup",
#                                      pattern = ".profile$",
#                                      full.names = TRUE)) %>%
#   separate(file, sep = "/|[.]", into = c("source", "method", "trait", "cvnum", "fold", "fill"), remove = FALSE) %>%
#   mutate(data = lapply(file, read.table, header = TRUE)) %>%
#   unnest(data) %>%
#   left_join(., training_coef, by = c("trait", "cvnum", "fold"))

head(gblup_pa)

save(gblup_h2, gblup_llik, gblup_pa, file = "reports/gblup.RData")
