## quick script to get snp lists
# use command line arguments
# usage: Rscript export_snplists.R path/to/file.txt

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
print(path)

dir <- paste0("data/processed/pathways/", unlist(strsplit(path, "/"))[4])
dir.create(dir)

# all genomic SNPs
all <- read.table("data/processed/input_nomissing.bim")[,2]

# SNPs in feature set
list1 <- read.table(path, header = TRUE)[,2]
write.table(list1, paste0(dir, "/list1"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# SNPs not in feature set
list2 <- all[!all %in% list1]
write.table(list2, paste0(dir, "/list2"), quote = FALSE,
            row.names = FALSE, col.names = FALSE)
