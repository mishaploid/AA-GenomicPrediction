################################################################################
## AA-GP: summarize multiBLUP results
## Sarah Turner-Hissong
## Updated 18 March 2019

library(tidyverse)
library(broom)

# phenotypes
pheno <- read.table("data/raw/pheno_file", header = TRUE) %>%
  gather(key = "trait",
         value = "phenotype",
         -FID, -IID)

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
  select(trait, pathway, Component, Share, SD)

head(multiblup_h2)

# likelihood ratio
multiblup_llik <- tibble(file = list.files(path = "models/multiblup_h2",
                                     pattern = ".reml$",
                                     full.names = TRUE, recursive = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "method2", "fill"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(multiblup_llik = as.numeric(as.character(X2))) %>%
  select(trait, pathway, multiblup_llik)

head(multiblup_llik)

# predictive ability

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
multiblup_pa <- tibble(file = list.files(path = "models/multiblup",
                                     pattern = ".pred$",
                                     full.names = TRUE,
                                     recursive = TRUE)) %>%
  separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "cvnum", "fold"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  left_join(pheno, ., by = c("FID" = "ID1", "IID" = "ID2", "trait")) %>%
  inner_join(., cross_val)

# training_coef <- tibble(file = list.files(path = "models/multiblup",
#                                           pattern = ".coeff$",
#                                           full.names = TRUE, recursive = TRUE)) %>%
#   separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "cvnum", "fold", "fill"), remove = FALSE) %>%
#   mutate(data = lapply(file, read.table, header = TRUE)) %>%
#   unnest(data) %>%
#   filter(Component %in% "Intercept") %>%
#   select(trait, pathway, cvnum, fold, Effect)
#
# head(training_coef)
#
# multiblup_pa <- tibble(file = list.files(path = "models/multiblup",
#                                          pattern = ".profile$",
#                                          full.names = TRUE, recursive = TRUE)) %>%
#   separate(file, sep = "/|[.]", into = c("source", "method", "pathway", "trait", "cvnum", "fold"), remove = FALSE) %>%
#   mutate(data = lapply(file, read.table, header = TRUE)) %>%
#   unnest(data) %>%
#   left_join(., training_coef, by = c("pathway", "trait", "cvnum", "fold"))

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
