################################################################################
## AA-GP: summarize GBLUP results
## Sarah Turner-Hissong
## Updated 18 March 2019

library(tidyverse)

load("reports/gblup.RData")

# read in amino acid trait ids
aa_names <- data.frame(trait = colnames(read.table("data/processed/pheno_file",
                                                   header = TRUE)[c(-1,-2)])) %>%
  mutate(id = as.factor(seq(1, n(), by = 1)))

## group into amino acid families (based on precursors)
aa_fam <- list(
  aspartate = data.frame(trait =
                           aa_names[grepl("asp|ile|lys|thr|met|Asp", aa_names$trait) &
                                      !grepl("GluFamCorr|SerFam", aa_names$trait), 1],
                         aa_fam = "aspartate"),
  glutamate = data.frame(trait =
                           aa_names[grepl("arg|glu|gln|pro|his|GluFam", aa_names$trait) &
                                      !grepl("SerFam", aa_names$trait),1],
                         aa_fam = "glutamate"),
  serine = data.frame(trait =
                        aa_names[grepl("gly|ser|SerFam", aa_names$trait), 1],
                      aa_fam = "serine"),
  pyruvate = data.frame(trait =
                          aa_names[grepl("ala|val|leu|BCAA|PyrFam", aa_names$trait) &
                                           !grepl("GluFamCorr|ile|SerFam", aa_names$trait), 1],
                        aa_fam = "pyruvate"),
  shikimate = data.frame(trait =
                           aa_names[grepl("phe|tyr|trp|ShikFam", aa_names$trait), 1],
                         aa_fam = "shikimate"),
  total = data.frame(trait = "Total", aa_fam = "total")
) %>%
  map_df(I) %>%
  unnest()

# group into categories by type of measurement (absolute, relative, ratio)
aa_cat <- list(absolute =
                 data.frame(trait = aa_names[!grepl("_",aa_names$trait), 1],
                            aa_cat = "absolute"),
  relative = data.frame(trait = aa_names[grepl("_t", aa_names$trait), 1], aa_cat = "relative"),
               family = data.frame(trait = aa_names[grepl("_|_shik|_pyr", aa_names$trait) &
                                                      !grepl("_t", aa_names$trait),1], aa_cat = "family")) %>%
  map_df(I)


files <- list.files("models/multiblup", pattern = ".profile$",
                    full.names = TRUE, recursive = TRUE)

multiblup <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("base", "method", "pathway", "ids"), remove = FALSE) %>%
  separate(ids, into = c("cv", "cvnum", "fold", "trait", "extension")) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
         #trait = gsub(".profile", "", trait)) %>%
  unnest(data) %>%
  group_by(pathway, cvnum, fold, trait) %>%
  mutate(cor = cor(Phenotype, Profile1, use = "complete.obs")) %>%
  summarise_at(., "cor", .funs = mean) %>%
  left_join(., aa_names, by = c("trait" = "id")) %>%
  left_join(., gblup_pa, by = c("trait.y" = "trait", "cvnum", "fold")) %>%
  select(pathway, cvnum, fold, trait.y, cor.x, cor.y) %>%
  rename(trait = trait.y,
         multiBLUP = cor.x,
         GBLUP = cor.y) %>%
  mutate(dif = multiBLUP - GBLUP)

save(multiblup, file = "reports/multiblup.RData")
