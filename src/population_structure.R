library(tidyverse)
library(adegenet)
library(ggfortify)
library(factoextra)
library(snpStats)

g <- read.plink(bed = "models/pca/ld_pruned.bed", bim = "models/pca/ld_pruned.bim",
                fam = "models/pca/ld_pruned.fam")

str(g)

geno <- as(g$genotypes, "numeric")
geno[1:5,1:5]

pca <- prcomp(geno)

pca_plot <- autoplot(pca, x = 1, y = 2)
pca_plot

summary(pca)

scree <- fviz_screeplot(pca, addlabels = TRUE, geom = "line", ncp = 20)
scree

screeloadings <- as.data.frame(pca$x) %>%
  mutate(IID = as.numeric(rownames(pca$x))) %>%
  select(IID, everything()) 

### regress phenotype data against principal components

pheno <- read.table("~/Documents/AA-GenomicPrediction/data/raw/pheno_file",
                    header = TRUE)

lminput <- left_join(pheno, screeloadings[,c(1:7)], by = "IID")

reformulate(names(lminput[-c(1:ncol(pheno))]), names(lminput[3]))

test <- apply(lminput[,3:ncol(pheno)], 2, function(e)
  residuals(lm(e ~ ., data = as.data.frame(e = pheno[,e], pca1$li[,c(1:6)]), na.action = na.exclude))) %>%
  as.data.frame() %>%
  cbind(pheno[,1:2], .)

write_delim(test, "~/Documents/AA-GenomicPrediction/data/raw/pheno_file_pcs")
test <- read.table("~/Documents/AA-GenomicPrediction/data/raw/pheno_file_pcs", header = TRUE)

rownames(test) <- pheno$IID

library(tidyverse)

one <- pheno[,-c(1:2)] %>%
  as.data.frame() %>%
  gather(key = trait, value = value) %>%
  ggplot(aes(trait, value)) +
  geom_boxplot()

two <- test[,-c(1:2)] %>%
  as.data.frame() %>%
  gather(key = trait, value = value) %>%
  ggplot(aes(trait, value)) +
  geom_boxplot()

library(cowplot)

plot_grid(one, two)


# g <- read.PLINK("~/Documents/AA-GenomicPrediction/models/pca/ld_pruned.raw")
# 
# grp <- find.clusters(g, max.n.clust = 40)
# 
# dapc.g <- dapc(g, grp$grp, n.pca = 200)
# 
# scatter(dapc.g)
# 
# xval <- xvalDapc(g, grp$grp, n.pca.max = 300, training.set = 0.9, result = "groupMean", 
#                  center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 30, xval.plot = TRUE)
# xval

# nice tutorial here https://speciationgenomics.github.io/pca/

# pca1 <- dudi.pca(g, scale = FALSE, scannf = FALSE, nf = 50)
# 
# dudi_scree <- screeplot(pca1)
# 
# s.label(pca1$li, xax = 1, yax = 2)
# 
# pca2 <- prcomp(g)
