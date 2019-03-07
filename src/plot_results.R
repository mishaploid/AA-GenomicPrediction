library(dplyr)
library(tidyverse)
library(gridExtra)
library(data.table)
library(quantreg)

################ PATHWAY RESULTS

# load null distribution
load("data/processed/null_dist2.RData")

null_dist <- left_join(null_dist, feature_sizes, by = "pathway") %>%
  select(aa, pathway, size.y, Heritability, LR) %>%
  rename(size = size.y)

# split aa traits into categories
aa_cats <- list(aspartate = data.frame(aa = aa_names[grepl("asp|ile|lys|thr|met|Asp", aa_names$aa) &
                                                       !grepl("GluFamCorr", aa_names$aa), 1], aa_cat = "aspartate"),
                glutamate = data.frame(aa = aa_names[grepl("arg|glu|gln|pro|GluFam", aa_names$aa) &
                                         !grepl("his", aa_names$aa), 1], aa_cat = "glutamate"),
                serine = data.frame(aa = aa_names[grepl("gly|ser|SerFam", aa_names$aa), 1], aa_cat = "serine"),
                pyruvate = data.frame(aa = aa_names[grepl("ala|val|leu|BCAA|PyrFam", aa_names$aa) &
                                                      !grepl("GluFamCorr", aa_names$aa), 1], aa_cat = "pyruvate"),
                shikimate = data.frame(aa = aa_names[grepl("phe|tyr|trp|ShikFam", aa_names$aa), 1], aa_cat = "shikimate"),
                histidine = data.frame(aa = aa_names[grepl("his", aa_names$aa), 1], aa_cat = "histidine")) %>%
  map_df(I)

# list file paths of reml output for pathways
files <- list.files(path = "models/reml_h2", pattern = "h2_\\d+.reml$", full.names = TRUE, recursive = TRUE)

h2 <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("source", "sub", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE, skip = 13)) %>%
  unnest(data) %>%
  mutate(id = str_replace(id, "reml_h2_", ""),
         id = str_extract(id, "\\d+")) %>%
  left_join(., aa_names, by = "id") %>%
  select(aa, pathway, Component, Heritability, Her_SD) %>%
  filter(Component %in% c("Her_K1", "Her_K2", "Her_ALL")) %>% #& Heritability >= 0) %>%
  mutate(Component = factor(Component, levels = c("Her_K1", "Her_K2", "Her_ALL")))

h3 <- h2 %>% 
  select(aa, pathway, Component, Heritability) %>%
  spread(key = Component, value = Heritability) %>%
  mutate(order = Her_K1) %>%
  gather(key = Component, value = Heritability, Her_K1, Her_K2, Her_ALL) %>%
  select(aa, pathway, Component, order) %>%
  left_join(h2, ., by = c("aa", "pathway", "Component")) %>%
  mutate(Component = factor(Component, levels = c("Her_K1", "Her_K2", "Her_ALL")),
         pathway = factor(pathway)) %>%
  left_join(., aa_cats, by = "aa")


########### determine significant pathways
# first, need size of specific pathways
files <- list.files("pathways", pattern = "list1", recursive = TRUE, full.names = TRUE)

pathways <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("source", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = FALSE)) %>%
  unnest(data) %>%
  count(pathway)

write.table(pathways, "reports/snp_count_pathways.txt", row.names = FALSE)


## need to run predict.rqss for each amino acid and each pathway

lr_sig <- function(x, y) {
  return(tryCatch(predict.rqss(x, data.frame(size = y)),
                  error = function(e) NA))
}

lr_fit <- NULL

for (i in unique(null_dist$aa)) {
  lr_fit[[i]] <- rqss(as.numeric(LR) ~ qss(size, constraint = "I"), 
                      data = na.omit(null_dist[null_dist$aa==i,]), tau = 0.95, lambda = 500)
}

h2_fit <- NULL

for (i in unique(null_dist$aa)) {
  h2_fit[[i]] <- rqss(as.numeric(Heritability) ~ qss(size, constraint = "I"), 
                      data = null_dist[null_dist$aa==i,], tau = 0.95, lambda = 500)
}

lr_test <- NULL
h2_test <- NULL

for(i in 1:nrow(pathways)) {
  lr_test[[i]] <- lapply(lr_fit, lr_sig, pathways[[i,2]])
  h2_test[[i]] <- lapply(h2_fit, lr_sig, pathways[[i,2]])
}

names(lr_test) <- unique(pathways$pathway)
names(h2_test) <- unique(pathways$pathway)

# likelihood ratio for pathways
files <- list.files(path = "models/reml_h2", pattern = "h2_\\d+.reml$", 
                    full.names = TRUE, recursive = TRUE)

lr <- tibble(file = files) %>%
  separate(file, sep = "/", 
           into = c("source", "sub", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(id = str_replace(id, "reml_h2_", ""),
         id = str_replace(id, ".reml", ""),
         X2 = as.numeric(X2)) %>%
  left_join(., aa_names, by = "id") %>%
  left_join(., gblup, by = "id") %>%
  select(aa, id, pathway, Num_Kinships, X2, X1) %>%
  mutate(LR = 2 * (X2 - X1)) %>%
  select(pathway, aa, LR)

# determine cutoffs based on null distribution

lr_cutoff <- rbindlist(lr_test, idcol = "pathway") %>%
  gather(key = aa, value = lr_95, -pathway)

h2_cutoff <- rbindlist(h2_test, idcol = "pathway") %>%
  gather(key = aa, value = h2_95, -pathway)


h4 <- h3 %>%
  filter(Component == "Her_K1") %>%
  left_join(., lr, by = c("pathway", "aa")) %>%
  left_join(., h2_cutoff, by = c("pathway", "aa")) %>%
  left_join(., lr_cutoff, by = c("pathway", "aa")) %>%
  select(aa, pathway, aa_cat, Heritability, Her_SD, LR, h2_95, lr_95) %>%
  mutate(h2_pass = Heritability > h2_95,
         lr_pass = LR > lr_95) %>%
  filter(Heritability > 0, 
         h2_pass == TRUE, 
         lr_pass == TRUE) 

# split pathways into categories
primary <- filter(pathways[,1], grepl("aa_|glycolysis|tca_cycle", pathway))
protein <- filter(pathways[,1], grepl("protein|degradation_", pathway)) 
secondary <- filter(pathways[,1], !grepl("aa_|glycolysis|tca|protein|degradation", pathway)) 

# define categories for plotting
cats <- list(primary = primary, protein = protein, secondary = secondary) %>%
  map_df(I, .id = "category")

h4 <- left_join(h4, cats, by = "pathway")
#   fct_relevel(pathway, "AAS_subset", "aa_transport", "AspartateFam", "BCAAFam", "bcat_genes", "GlutamineFam", 
#               "Shikimate", "glycolysis", "tca_cycle",)
  

ggplot(h4, aes(aa, Heritability, col = pathway)) +
  geom_point(position = position_dodge(width = 0.3)) +
#  scale_shape_manual(values = 1:length(unique(h4$pathway))) +
  geom_errorbar(aes(ymin = Heritability - Her_SD, ymax = Heritability + Her_SD), width = 0.1,
                position = position_dodge(width = 0.3)) +
  facet_grid(aa_cat ~ category) + # scales = "free_x", ncol = 10) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


h4$category <- factor(h4$category, levels = unique(h4$category))

ggplot(h4, aes(aa, reorder(pathway, -as.numeric(category)), col = aa_cat, size = Heritability)) +
  geom_point() +
  scale_colour_discrete(name = "Family") + 
  facet_grid(category ~ aa_cat, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing.y = unit(0.5, "lines"),
        panel.grid.major = element_line(colour = 'lightgray', size = 0.25),
        text = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  labs(y = "", x = "Amino Acid") 

ggsave("reports/figures/heritability_by_family_v2.png", height = 4, width = 8)


# chisq test for LR 
foo <- h3 %>%
  filter(Component == "Her_K1") %>%
  left_join(., lr, by = c("pathway", "aa")) %>%
  left_join(., h2_cutoff, by = c("pathway", "aa")) %>%
  left_join(., lr_cutoff, by = c("pathway", "aa")) %>%
  dplyr::select(aa, pathway, aa_cat, Heritability, Her_SD, LR, h2_95, lr_95) %>%
  mutate(h2_pass = Heritability > h2_95,
         lr_pass = LR > lr_95) %>%
  group_by(aa) %>%
  mutate(chisq = pchisq(LR, df = 1, lower.tail = FALSE),
            x_fdr = p.adjust(chisq, method = "fdr")) 

bar <- foo %>%
  filter(aa == "asp_t") %>%
  mutate(chisq = pchisq(LR, df = 1, lower.tail = FALSE),
         x_fdr = p.adjust(chisq, method = "fdr")) 

### PLOT SCRATCH SPACE 

h5 <- split(h4, h4$aa_cat)
p <- NULL

for (i in names(h5)){
p[[i]] <- ggplot(h5[[i]], aes(aa, pathway, col = category, size = Heritability)) +
  geom_point() +
  facet_wrap(. ~ category) +
  ggtitle(paste(i))
}
p[[2]]

absolute <- ggplot(h5[[1]], aes(aa, Heritability, col = pathway, shape = pathway)) +
  geom_point(position = position_dodge(width = 0.3)) +
  scale_shape_manual(values = 1:length(unique(h4$pathway))) +
  geom_errorbar(aes(ymin = Heritability - Her_SD, ymax = Heritability + Her_SD), width = 0.1,
                position = position_dodge(width = 0.3)) +
  facet_grid(. ~ category.y, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

absolute

family <- ggplot(h5[[2]], aes(aa, Heritability, col = pathway, shape = pathway)) +
  geom_point(position = position_dodge(width = 0.3)) +
  scale_shape_manual(values = 1:length(unique(h4$pathway))) +
  geom_errorbar(aes(ymin = Heritability - Her_SD, ymax = Heritability + Her_SD), width = 0.1,
                position = position_dodge(width = 0.3)) +
  facet_grid(. ~ category, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

family


relative <- ggplot(h5[[3]], aes(aa, Heritability, col = pathway, shape = pathway)) +
  geom_point(position = position_dodge(width = 0.3)) +
  scale_shape_manual(values = 1:length(unique(h4$pathway))) +
  geom_errorbar(aes(ymin = Heritability - Her_SD, ymax = Heritability + Her_SD), width = 0.1,
                position = position_dodge(width = 0.3)) +
  facet_grid(. ~ category.y, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

relative


legend_plot <- ggplot(h4, aes(aa, Heritability, col = pathway, shape = pathway)) +
  geom_point(position = position_dodge(width = 0.3)) +
  scale_shape_manual(values = 1:length(unique(h4$pathway))) +
  theme_bw() +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))

legend <- get_legend(legend_plot)

library(cowplot)

pgrid <- plot_grid(absolute, relative, family, ncol = 1)
plot_grid(pgrid, legend, ncol = 2, rel_widths = c(1, .3), labels = levels(h4$category.y))

ggsave("figures/multiblup_results.png", height = 6, width = 10)




















# p <- list()

# for (i in levels(h3$pathway)) {
#   tmp <- h3 %>%
#     filter(pathway %in% i) 
#   head(tmp)
#   p[[i]] <- ggplot(tmp, aes(reorder(aa, -order), Heritability, colour = Component)) +
#     geom_point(position = position_dodge(width = 0.3)) +    
#     geom_errorbar(aes(ymin = Heritability - Her_SD, ymax = Heritability + Her_SD), width = 0.1,
#                   position = position_dodge(width = 0.3)) + 
#     facet_wrap(. ~ category, ncol = 1, scales = "free_x") +
#     geom_hline(aes(yintercept = 0)) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#     ggtitle(paste0(i))
#   ggsave(paste0("figures/", i, ".png"))
# }


