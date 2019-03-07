library(dplyr)
library(tidyverse)
library(gridExtra)
library(data.table)
library(quantreg)

## read in amino acid names + numeric codes 
aa_names <- data.frame(aa = colnames(read.table("data/processed/pheno_file", header = TRUE)[3:56]), id = as.character(1:54))

## read in size of each feature group from control sets 

trait <- "his"
trait_pattern <- paste0(trait, "_\\d+")

feature_sizes <- tibble(file = list.files(path = "../AA-GenomicPrediction/gfblup_null", pattern = trait_pattern, full.names = TRUE)) %>%
  separate(file, sep = "/", into = c("dir", "source", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  mutate(pathway = str_replace(id, "his_", ""),
         pathway = str_replace(pathway, ".txt", ""),
         pathway = paste0("c_", pathway)) %>%
  select(pathway, size)

# read in gblup results 
gblup <- tibble(file = list.files(path = "models/reml_old_results/old/gblup",
                                  pattern = "h2.\\d+.reml$", full.names = TRUE)) %>%
  separate(file, sep = "/", into = c("source", "sub", "dir", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(id = str_replace(id, "reml_h2.", ""),
         id = str_replace(id, ".reml", ""),
         X1 = as.numeric(X1)) %>%
  select(id, X1)

files <- list.files(path = "reml_null", pattern = "h2_\\d+.reml$", 
                    full.names = TRUE, recursive = TRUE)

# calculate likelihood ratio for null distribution

lr_null <- tibble(file = files) %>%
  separate(file, sep = "/", 
           into = c("source", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(id = str_replace(id, "reml_h2_", ""),
         id = str_replace(id, ".reml", ""),
         X2 = as.numeric(X2)) %>%
  left_join(., aa_names, by = "id") %>%
  left_join(., gblup, by = "id") %>%
  select(aa, id, pathway, Num_Kinships, X2, X1) %>%
  left_join(., feature_sizes, by = "pathway") %>%
  mutate(group_size = cut(size, breaks = seq(0, 130000, by = 10000), 
                          include.lowest = TRUE, dig.lab = 10)) %>%
  mutate(LR = 2 * (X2 - X1)) %>%
  filter(LR > 0)

# check distribution of likelihood ratio
ggplot(lr_null, aes(LR)) + 
  geom_histogram(bins = 50) + 
  facet_wrap(. ~ aa, ncol = 10, scales = 'free_x')

ggsave("figures/lr_histogram.png", height = 6, width = 15)

# compare to chisq with 1 df

foo <- lr_null %>% filter(size <= 50000)
ggplot(foo, aes(sample = LR, colour = group_size)) + 
  stat_qq(distribution = qchisq, dparams = list(df = 1)) +
  facet_wrap(. ~ aa, ncol = 8, scales = 'free_y') +
  geom_abline(slope = 1, intercept = 0) 

ggsave("figures/lr_qqchisq.png", height = 15, width = 15)

## 95% LR threshold for plotting

lr_fit <- NULL
coefs <- NULL

for (i in unique(lr_null$aa)) {
  lr_fit[[i]] <- rqss(as.numeric(LR) ~ qss(size, constraint = "I"), 
                      data = lr_null[lr_null$aa==i,], tau = 0.95, lambda = 500)
  coefs[[i]] <- data.frame(lr_fit[[i]]$qss$size$xyz) %>%
    mutate(lr_95 = X2 + lr_fit[[i]]$coef[1],
           size = X1) %>%
    select(size, lr_95)
}

lr_threshold <- rbindlist(coefs, idcol = "aa")


lr_null <- lr_null %>%
  left_join(., lr_threshold, by = c("aa", "size"))

# 
# # extract threshold for plotting
# coefs <- data.frame(lr_fit$qss$size$xyz) %>%
#   mutate(lr_95 = X2 + lr_fit$coef[1], # add intercept to coefficients
#          size = X1) %>%
#   select(size, lr_95)

# read in genomic heritability for null 
h2_null <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("source", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE, skip = 13)) %>%
  unnest(data) %>%
  mutate(id = str_replace(id, "reml_h2_", ""),
         id = str_extract(id, "\\d+")) %>%
  left_join(., aa_names, by = "id") %>%
  select(aa, pathway, Component, Heritability, Her_SD) %>%
  filter(Component %in% c("Her_K1")) 

# h2_null <- h2_null %>% 
#   left_join(., feature_sizes, by = "pathway") 

# plotdf <- left_join(h2_null, lr_null, by = c("aa", "pathway")) %>%
#   mutate(col = LR >= lr_95) %>%
#   select(aa, pathway, size, Component, Heritability, Her_SD, col)

plotdf <- left_join(h2_null, lr_null, by = c("aa", "pathway")) %>%
#  mutate(col = LR >= lr_95) %>%
  left_join(., feature_sizes, by = "pathway") %>%
  select(aa, pathway, size.y, Component, Heritability, Her_SD, LR) %>%
  rename(size = size.y)

null_dist <- plotdf %>% 
  select(aa, pathway, size, Heritability, LR)

plotdf <- plotdf %>%
  filter(Heritability >= 0,
         size <= 50000)

ggplot(plotdf, aes(size, Heritability, colour = col)) +
  geom_point(position = position_dodge(width = 0.3)) +    
  #   geom_errorbar(aes(ymin = Heritability - Her_SD, ymax = Heritability + Her_SD), width = 0.1,
  #                 position = position_dodge(width = 0.3)) +
  geom_quantile(method = "rqss", lambda = 500, quantiles = c(0.95),
                formula = as.formula(y ~ qss(x, constraint = "I")),
                aes(linetype = factor(..quantile..)), colour = "black") + 
  facet_wrap(. ~ aa, ncol = 10, scales = "free_x") +
  geom_hline(aes(yintercept = 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("figures/null_dist_all.png", height = 16, width = 16)

save(aa_names, null_dist, gblup, feature_sizes, file = "data/processed/null_dist2.RData")

# ggplot(cs, aes(size.x, h2, colour = col, size = col)) +
#   geom_point() +
#   scale_size_manual(values = c(1, 1.75), guide = FALSE) +  
#   scale_colour_manual(name = "", values = c('gray85', 'darkcyan'), 
#                       labels = c(expression('LR < LR'[95]), expression('LR > LR'[95]))) + 
#   geom_quantile(method = "rqss", lambda = 500, quantiles = c(0.95, 0.50),
#                 formula = as.formula(y ~ qss(x, constraint = "I")),
#                 aes(linetype = factor(..quantile..)), size = 1.5, colour = "darkcyan") + 
#   scale_linetype_manual(name = expression('H'[set]^2 * ' percentile'), values = c(2, 1)) +
#   geom_line(aes(size, naive), colour = "black", size = 1, linetype = "dotdash") +
#   guides(linetype = guide_legend(override.aes = list(size = 1))) + 
# #  ggtitle(paste(trait)) +
#   xlab("Group size") + 
#   ylab(expression('H'[set]^2)) + 
#   theme_bw() +
#   theme(legend.key.width = unit(1, "cm"))
# 