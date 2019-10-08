################################################################################
## AA-GP: PCA for SNP data
## Sarah Turner-Hissong
## Updated 06 October 2019

library(tidyverse)
library(adegenet)
library(ggfortify)
library(factoextra)
library(snpStats)

g <- read.plink(bed = "models/pca/ld_pruned.bed",
                bim = "models/pca/ld_pruned.bim",
                fam = "models/pca/ld_pruned.fam")

str(g)

geno <- as(g$genotypes, "numeric")
geno[1:5,1:5]

pca <- prcomp(geno)

pca_plot <- autoplot(pca, x = 1, y = 2)

summary(pca)

scree <- fviz_screeplot(pca, addlabels = TRUE, geom = "line", ncp = 20)

screeloadings <- as.data.frame(pca$x) %>%
  mutate(IID = as.numeric(rownames(pca$x))) %>%
  select(IID, everything())

# read in phenotype data
pheno <- read.table("data/raw/pheno_file",
                    header = TRUE)

# create input data for linear regression
lminput <- left_join(pheno, screeloadings[,c(1:7)], by = "IID")

reformulate(names(lminput[-c(1:ncol(pheno))]), names(lminput[3]))

# regress phenotype data against principal components
# use 'apply' to perform over each column (phenotype)
# function exports residuals from linear regression against the first
# six principal components from lminput
# then use cbind to add individual id information in PLINK format
pheno_pc <- apply(lminput[,3:ncol(pheno)], 2, function(e)
  residuals(lm(e ~ ., data = as.data.frame(e = pheno[,e], pca1$li[,c(1:6)]), na.action = na.exclude))) %>%
  as.data.frame() %>%
  cbind(pheno[,1:2], .)

write_delim(pheno_pc, "~/Documents/AA-GenomicPrediction/data/raw/pheno_file_pcs")
test <- read.table("~/Documents/AA-GenomicPrediction/data/raw/pheno_file_pcs", header = TRUE)

rownames(test) <- pheno$IID
