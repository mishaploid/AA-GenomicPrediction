################################################################################
## AA-GP: summarize GBLUP results
## Sarah Turner-Hissong
## Updated 18 March 2019

library(tidyverse)

# heritability
multiblup_h2 <- tibble(file = list.files(path = "models/multiblup_h2",
                                     pattern = ".reml$",
                                     full.names = TRUE, recursive = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "method2", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE, skip = 13)) %>%
  unnest(data) %>%
  filter(Component %in% c("Her_K1", "Her_K2")) %>%
  select(trait, pathway, Component, Heritability, Her_SD)

head(multiblup_h2)

# likelihood ratio
multiblup_lr <- tibble(file = list.files(path = "models/multiblup_h2",
                                     pattern = ".reml$",
                                     full.names = TRUE, recursive = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "method2", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(multiblup_llik = as.numeric(X2)) %>%
  select(trait, pathway, multiblup_llik)

head(multiblup_lr)

# predictive ability
multiblup_pa <- tibble(file = list.files(path = "models/multiblup",
                                         pattern = ".profile$",
                                         full.names = TRUE, recursive = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "cvnum", "fold"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  group_by(pathway, cvnum, fold, trait) %>%
  mutate(multiblup_cor = cor(Phenotype, Profile1, use = "complete.obs")) %>%
  summarise_at(., "multiblup_cor", .funs = mean) %>%
  ungroup() %>%
  select(trait, pathway, cvnum, fold, multiblup_cor)

head(multiblup_pa)

save(multiblup_h2, multiblup_lr, multiblup_pa, file = "reports/multiblup.RData")
