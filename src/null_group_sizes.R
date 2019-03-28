
# uniform distribution of markers up to 50000
gene_group_sizes <- round(runif(5000, 1, 50000))
hist(gene_group_sizes)
write.table(gene_group_sizes, "data/interim/null_group_sizes.txt", row.names = FALSE)

# write sampling file (control_sampling.txt) for sbatch
sampling <- as.data.frame(seq(1, length(gene_group_sizes), by = 50))
colnames(sampling) <- "start"
sampling$end <- sampling$start + 49

write.table(sampling, "data/interim/control_sampling.txt", row.names = FALSE, col.names = FALSE)
