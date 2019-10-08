################################################################################
## AA-GP: summarize GBLUP results
## Sarah Turner-Hissong
## Updated 18 March 2019

library(tidyverse)
library(broom)

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

# predictive accuracy

training_coef <- tibble(file = list.files(path = "models/gblup",
                                          pattern = ".coeff$",
                                          full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "trait", "cvnum", "fold", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  select(trait, cvnum, fold, Effect)

head(training_coef)

gblup_pa <- tibble(file = list.files(path = "models/gblup",
                                     pattern = ".profile$",
                                     full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "trait", "cvnum", "fold", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  left_join(., training_coef, by = c("trait", "cvnum", "fold"))

head(gblup_pa)

save(gblup_h2, gblup_llik, gblup_pa, file = "reports/gblup.RData")
