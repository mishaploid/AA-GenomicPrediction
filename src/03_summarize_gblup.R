################################################################################
## AA-GP: summarize GBLUP results
## Sarah Turner-Hissong
## Updated 18 March 2019

library(tidyverse)

# read in amino acid trait ids
# aa_names <- data.frame(trait = colnames(read.table("data/processed/pheno_file",
#                                                    header = TRUE))) %>%
#   filter(!trait %in% c("trait", "FID", "IID")) %>%
#   mutate(id = as.character(1:length(trait)))

# heritability
gblup_h2 <- tibble(file = list.files(path = "models/gblup_h2",
                                     pattern = ".reml$",
                                     full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "metric", "trait", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE, skip = 13)) %>%
  unnest(data) %>%
  filter(Component %in% "Her_ALL") %>%
  select(trait, Heritability, Her_SD)

# log likelihood
gblup_llik <- tibble(file = list.files(path = "models/gblup_h2",
                                     pattern = ".reml$",
                                     full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "metric", "trait", "fill", "fill2"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(gblup_llik = as.numeric(X1)) %>%
  select(trait, gblup_llik)


# predictive accuracy

gblup_pa <- tibble(file = list.files(path = "models/gblup",
                                     pattern = ".profile$",
                                     full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "trait", "cvnum", "fold", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  group_by(trait, cvnum, fold) %>%
  mutate(cor = cor(Phenotype, Profile1, use = "complete.obs")) %>%
  summarise_at(., "cor", .funs = mean) %>%
  ungroup() %>%
  select(trait, cvnum, fold, cor)

save(gblup_h2, gblup_llik, gblup_pa, file = "reports/gblup.RData")
