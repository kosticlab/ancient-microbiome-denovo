library(tidyverse)

#Name and location of the result plot
plot_name <- "barplot.pdf"
#Path to the kraken folders run on different thresholds
kraken15_dir <- "kraken15"
kraken30_dir <- "kraken03"
kraken60_dir <- "kraken06"
kraken90_dir <- "kraken09"
#File containing the ids of all parasites in the database
parasite_ids_file <- "parasite_ids.tsv"
#Samples to extract and visualize
sample_ids <- c("AW107", "AW108", "AW110A", "UT30.3", "UT43.2", "Zape1", "Zape2", "Zape3", "1026.1.4", "1043.3.1", "3567.1.1")



read_kraken_reprot_helper <- function(file){
    df <- read_tsv(file,  col_names = FALSE)
    colnames(df) <- c("percentage", "subtree_count", "node_count", "level", "taxonomic_id", "taxonomic_name")
    df$sample <- basename(file)
    return(df)
}

read_kraken_report <- function(directory){
    files <- list.files(directory, full.names = TRUE)
    df <- map_df(files, read_kraken_reprot_helper)
    df <- df %>% dplyr::select(-percentage)
    return(df)
}

parasite_ids <- read_tsv(parasite_ids_file, col_names = FALSE)

df15 <- read_kraken_report(kraken15_dir)
df30 <- read_kraken_report(kraken30_dir)
df60 <- read_kraken_report(kraken60_dir)
df90 <- read_kraken_report(kraken90_dir)

df15 <- df15 %>% filter(taxonomic_id %in% parasite_ids$X1) %>% filter(sample %in% sample_ids) %>% mutate(threshold = 0.15)
df30 <- df30 %>% filter(taxonomic_id %in% parasite_ids$X1) %>% filter(sample %in% sample_ids) %>% mutate(threshold = 0.30)
df60 <- df60 %>% filter(taxonomic_id %in% parasite_ids$X1) %>% filter(sample %in% sample_ids) %>% mutate(threshold = 0.60)
df90 <- df90 %>% filter(taxonomic_id %in% parasite_ids$X1) %>% filter(sample %in% sample_ids) %>% mutate(threshold = 0.90)

at_least_1000 <- df30 %>% filter(subtree_count >= 1000) %>% select(taxonomic_name) %>% distinct()

df15 <- df15 %>% filter(taxonomic_name %in% at_least_1000$taxonomic_name)
df30 <- df30 %>% filter(taxonomic_name %in% at_least_1000$taxonomic_name)
df60 <- df60 %>% filter(taxonomic_name %in% at_least_1000$taxonomic_name)
df90 <- df90 %>% filter(taxonomic_name %in% at_least_1000$taxonomic_name)

df <- bind_rows(df15,df30,df60,df90)

farben <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")

ggplot(data = df, aes(taxonomic_name, subtree_count, group = sample)) +
    geom_col(aes(fill = fct_rev(as.factor(threshold))), position="dodge") + theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90),
        panel.background = element_blank(),
        strip.background = element_rect(fill="white")) +
    geom_hline(yintercept = 0) + facet_wrap(~sample) +
    scale_fill_manual(name = "threshold",values = farben)+
    xlab("") + ylab("read count") + ggsave(plot_name)



