library(biomaRt)

mart <- useMart(biomart = "plants_mart", host = "plants.ensembl.org", dataset ="athaliana_eg_gene")

foo <- data_frame(chr = rep(1:5), length = c(34.96, 22.04, 25.50, 20.86, 31.27))

files <- list.files("pathways", pattern = ".txt", recursive = TRUE, full.names = TRUE)

pathways <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("source", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, fread, select = c(1:5))) %>%
  unnest(data)

gene_list <- getBM(attributes = c("ensembl_gene_id", "start_position", "end_position"),
                   filters = "ensembl_gene_id",
                   values = pathways$ensembl_gene_id.y,
                   mart = mart)

test <- left_join(pathways, gene_list, by = c("ensembl_gene_id.y" = "ensembl_gene_id")) %>%
  dplyr::select(pathway, Chromosome, ensembl_gene_id.y, start_position, end_position) %>%
  rename(ensembl_gene_id = ensembl_gene_id.y) %>%
  group_by(pathway) %>%
  distinct(pathway, ensembl_gene_id, .keep_all = TRUE)

test2 <- test %>% dplyr::select(pathway, ensembl_gene_id) %>% split(.$pathway) %>%
  lapply(function(x) x[(names(x) %in% "ensembl_gene_id")]) 

# ggplot(foo, aes(y = chr)) +
#  geom_segment(data = foo, aes(x = 0, xend = length, y = chr, yend = chr)) +
# geom_point(data = test, aes(start_position/1e6, Chromosome, col = pathway), size = .1,
#            position = "dodge")

# pathway categories
primary <- filter(test[,1], grepl("aa_|glycolysis|tca_cycle", pathway))
protein <- filter(test[,1], grepl("protein|degradation_", pathway)) 
secondary <- filter(test[,1], !grepl("aa_|glycolysis|tca|protein|degradation", pathway)) 

cats <- list(primary = unique(primary), protein = unique(protein), secondary = unique(secondary)) %>%
  map_df(I, .id = "category")

plot_df <- left_join(test, cats, by = "pathway")

ggplot(plot_df, aes(start_position/1e6, pathway, col = pathway)) +
  geom_point(shape = 3) +
  facet_wrap(Chromosome ~ category, ncol = 3, scales = "free_y") +
  theme(legend.position = "none") +
  ylab("") + xlab("")
ggsave("reports/figures/snp_locations.png", height = 10, width = 12)