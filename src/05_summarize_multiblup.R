################################################################################
## AA-GP: summarize multiBLUP results
## Sarah Turner-Hissong
## Updated 18 March 2019

library(tidyverse)
library(broom)

# heritability
multiblup_h2 <- tibble(file = list.files(path = "models/multiblup_h2",
                                     pattern = ".share$",
                                     full.names = TRUE, recursive = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "method2", "fill"), remove = FALSE) %>%
  # mutate(data = lapply(file, read.table, header = TRUE, skip = 13)) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  # filter(Component %in% c("Her_K1", "Her_K2")) %>%
  # select(trait, pathway, Component, Heritability, Her_SD)
  select(trait, pathway, Component, Share, SD) %>%
  # correct trait names
  mutate(trait = gsub("Corr", "", trait),
         trait = gsub("_asp", "", trait))

head(multiblup_h2)

# likelihood ratio
multiblup_llik <- tibble(file = list.files(path = "models/multiblup_h2",
                                     pattern = ".reml$",
                                     full.names = TRUE, recursive = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "method2", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(multiblup_llik = as.numeric(X2)) %>%
  select(trait, pathway, multiblup_llik) %>%
  # correct trait names
  mutate(trait = gsub("Corr", "", trait),
         trait = gsub("_asp", "", trait))

head(multiblup_llik)

# predictive ability
training_coef <- tibble(file = list.files(path = "models/gblup",
                                          pattern = ".coeff$",
                                          full.names = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "trait", "cvnum", "fold", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  select(trait, cvnum, fold, Effect) %>%
  # correct trait names
  mutate(trait = gsub("Corr", "", trait),
         trait = gsub("_asp", "", trait))

head(training_coef)

multiblup_pa <- tibble(file = list.files(path = "models/multiblup",
                                         pattern = ".profile$",
                                         full.names = TRUE, recursive = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "cvnum", "fold"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  # correct trait names
  mutate(trait = gsub("Corr", "", trait),
         trait = gsub("_asp", "", trait)) %>%
  left_join(., training_coef, by = c("trait", "cvnum", "fold")) 

  # nest(-trait, -pathway, -cvnum, -fold) %>%
  # mutate(fit = map(data, ~ lm(Profile1 ~ Phenotype, data = .)), results = map(fit, augment)) %>%
  # unnest(results) %>%
  # group_by(trait, pathway, cvnum, fold) %>%
  # summarize(multiblup_cor = cor(Profile1, Phenotype, use = "complete.obs"),
  #   multiblup_mse = sum((.resid/(1-.hat))^2)/length(.resid)) %>%
  # select(trait, pathway, cvnum, fold, multiblup_cor, multiblup_mse) %>%
  # ungroup()
#   group_by(pathway, cvnum, fold, trait) %>%
#   summarize(multiblup_cor = cor(Phenotype, Profile1, use = "complete.obs"),
#          multiblup_mse = mse(Phenotype, Profile1),
#          multiblup_rmse = rmse(Phenotype, Profile1)) %>%
# #  mutate(multiblup_cor = cor(Phenotype, Profile1, use = "complete.obs")) %>%
# #  summarise_at(., "multiblup_cor", .funs = mean) %>%
#   ungroup() %>%
#   select(trait, pathway, cvnum, fold, multiblup_cor, multiblup_mse, multiblup_rmse)

head(multiblup_pa)

# number of snps/genes in each pathway
files <- list.files(path = "data/processed/pathways", pattern = "*.txt", full.names = TRUE,
                    recursive = TRUE)

pathways <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("dir", "dir2", "dir3", "pathway", "pathway2"),
           remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE),
         data = lapply(data, mutate_all, as.character)) %>%
  unnest(data) %>%
  group_by(pathway) %>%
  summarise(size = length(unique(snp_id)),
            n_genes = length(unique(ensembl_gene_id.x))) %>%
  select(pathway, n_genes, size) %>%
  ungroup() %>%
  write_csv(., "reports/pathways_summary.csv")

head(pathways)

save(multiblup_h2, multiblup_llik, multiblup_pa, pathways, file = "reports/multiblup.RData")
