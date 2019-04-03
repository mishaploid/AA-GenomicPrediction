################################################################################
## AA-GP: summarize GBLUP results
## Sarah Turner-Hissong
## Updated 18 March 2019

library(tidyverse)
library(broom)

# read in amino acid trait ids
# aa_names <- data.frame(trait = colnames(read.table("data/processed/pheno_file",
#                                                    header = TRUE))) %>%
#   filter(!trait %in% c("trait", "FID", "IID")) %>%
#   mutate(id = as.character(1:length(trait)))

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
  mutate(gblup_llik = as.numeric(X1)) %>%
  select(trait, gblup_llik)

head(gblup_llik)

# predictive accuracy

gblup_pa <- tibble(file = list.files(path = "models/gblup",
                                     pattern = ".profile$",
                                     full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "trait", "cvnum", "fold", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  nest(-trait, -cvnum, -fold) %>%
  mutate(fit = map(data, ~ lm(Profile1 ~ Phenotype, data = .)), results = map(fit, augment)) %>%
  unnest(results) %>%
  group_by(trait, cvnum, fold) %>%
  summarize(gblup_cor = cor(Profile1, Phenotype, use = "complete.obs"),
    gblup_mse = sum((.resid/(1-.hat))^2)/length(.resid)) %>%
  select(trait, cvnum, fold, gblup_cor, gblup_mse) %>%
  ungroup()
  # unnest(data) %>%
  # group_by(trait, cvnum, fold) %>%
  # mutate(cor = cor(Phenotype, Profile1, use = "complete.obs")) %>%
  # summarise_at(., "cor", .funs = mean) %>%
  # ungroup() %>%
  # select(trait, cvnum, fold, cor)

head(gblup_pa)

save(gblup_h2, gblup_llik, gblup_pa, file = "reports/gblup.RData")
