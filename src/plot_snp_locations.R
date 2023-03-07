### pull corresponding SNP id and plot locations 

library(biomaRt)
library(tidyverse)
library(data.table)

setwd("~/Documents/AA-GenomicPrediction/")

mart <- useMart(biomart = "plants_mart", host = "plants.ensembl.org", dataset ="athaliana_eg_gene")

foo <- data_frame(chromosome = rep(1:5), length = c(34.96, 22.04, 25.50, 20.86, 31.27))

files <- list.files("data/processed/pathways", pattern = ".txt", recursive = TRUE, full.names = TRUE)

pathways <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("source", "dir", "id", "pathway"), remove = FALSE) %>%
  mutate(data = lapply(file, fread, select = c(1:5))) %>%
  unnest(data) %>%
  filter(!pathway %in% c("aa_synthesis.txt", "aa_degradation.txt")) %>%
  mutate(pathway = gsub(".txt", "", pathway),
         pathway = fct_recode(pathway,
                             aa_degradation = "aa_deg_bcat_haplo",
                             aa_synthesis = "aa_syn_nobcat_haplo"),
         category = ifelse(
           grepl("degradation_|protein_|protease", pathway),
           "protein",
           ifelse(
             grepl(
               "flavonoids|isoprenoids|N_containing|phenylpropanoids|S_containing",
               pathway
             ),
             "specialized",
             ifelse(
               grepl("aa_degradation|aa_synthesis|aa_transport", pathway),
               "amino acid",
               "core"
             )
           )
         ) ,
         category = factor(
           category,
           levels = c("amino acid", "core", "specialized", "protein")
         )
  )


pathways <- pathways %>%
  complete(chromosome, nesting(pathway, category)) %>% 
  left_join(., foo, by = "chromosome") %>% 
  mutate(chromosome = paste0("chr", chromosome))


library(RColorBrewer)

ggplot(pathways, aes(x = position/1e6,
                     y = pathway,
                     color = category,
                     alpha = 0.1)) + 
  geom_segment(aes(x = 0,
                   xend = length,
                   y = pathway,
                   yend = pathway),
               inherit.aes = FALSE,
               color = "lightgray",
               size = 1) +
  geom_point(position = position_jitter(height = 0.2),
             size = 1,
             alpha = 0.1) +
  scale_color_brewer(palette = "Dark2") + 
  facet_grid(category ~ chromosome,
             scales = "free",
             space = "free") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  xlab("Position (Mb)")

# ggsave("reports/figures/FigS3.pdf", height = 5, width = 8)

### cover 

pathways %>%
  # filter(!pathway %in% "degradation_ubiquitin") %>% 
ggplot(., aes(x = position/1e6,
                     y = chromosome,
                     color = category,
                     alpha = 0.1)) + 
  geom_segment(aes(x = 0,
                   xend = length,
                   y = chromosome,
                   yend = chromosome),
               inherit.aes = FALSE,
               color = "gray50",
               size = 1) +
  geom_point(aes(fill = category),
             stat = "identity",
             position = position_jitterdodge(dodge.width = 0.7,
                                             jitter.width = 0.4),
             size = 0.1,
             alpha = 0.2) +
  scale_color_brewer(palette = "Dark2") + 
  theme_void() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "gray5",
                                        colour = "gray5",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "gray5"), 
        plot.margin = unit(c(5,1,1,1), "cm")) +
  xlab("Position (Mb)") +
  scale_x_reverse() + 
  coord_flip() 

ggsave("manuscript/aagp_cover_submission_full.png",
       width = 8.375,
       height = 10.875,
       units = "in", 
       dpi = 300)
