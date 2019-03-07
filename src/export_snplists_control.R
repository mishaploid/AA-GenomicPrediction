## quick script to get snp lists
# use command line arguments
# input files located in data/processed/random_sets
library(tidyverse)

# all genomic SNPs
all <- read.table("data/processed/input_nomissing.bim")[,2]

# bcat genes
bcat_genes <- read.table("data/processed/pathways/bcat_genes/bcat_genes.txt", header = TRUE)

## control sets
control <- list.files(path = "data/processed/random_sets", pattern = "*.txt", full.names = TRUE)

for (i in control) {
  tmp <- read.table(i) %>%
     distinct(., SNP_ID) # remove duplicates
  print(head(tmp))
  # write list 1 (snps in feature set)
  dir <- paste0("data/processed/random_sets/c_", parse_number(i))
  print(dir)
  dir.create(dir)
  write.table(tmp, paste0(dir, "/list1"), quote = FALSE,
                          col.names = FALSE, row.names = FALSE)
  # write list 2 (remaining genomic snps)
  tmp2 <- all[!all %in% tmp$SNP_ID]
  write.table(tmp2, paste0(dir, "/list2"), quote = FALSE,
                           col.names = FALSE, row.names = FALSE)
  ## exclude BCAT
  dir <- paste0("data/processed/random_sets_nobcat/c_", parse_number(i))
  print(dir)
  dir.create(dir)
  nobcat <- anti_join(tmp, bcat_genes, by = "SNP_ID")
  write.table(nobcat, paste0(dir, "/nobcat1"), quote = FALSE,
                             col.names = FALSE, row.names = FALSE)
  nobcat2 <- all[!all %in% nobcat$SNP_ID]
  write.table(nobcat2, paste0(dir, "/nobcat2"), quote = FALSE,
                              col.names = FALSE, row.names = FALSE)
}
