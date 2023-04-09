library(tidyverse)
library(vegan)

asv_table = read_csv("raw_data/asv_table.csv")

# parse raw sample_id with external python script
raw_sample_ids = colnames(asv_table[,2:ncol(asv_table)])
write_lines(raw_sample_ids, "output/raw_sample_ids.txt")
system2("/opt/miniconda3/bin/python", "parse_raw_sample_ids.py")

# rename sample_id in asv table
parsed_metadata = read_tsv("data/parsed_sample_metadata.tsv")
new_sample_ids = parsed_metadata$sample_id
names(new_sample_ids) = parsed_metadata$raw_sample_id
colnames(asv_table) = c("ASV_ID", new_sample_ids[raw_sample_ids])

# drop outlier ASVs in JX5P1-2
JX5 = select(asv_table, ASV_ID, starts_with("JX5"))
JX5 = JX5[rowSums(JX5[,2:ncol(JX5)])>0,]
JX5_ca = cca(t(JX5[,2:ncol(JX5)])>0)
outlier_id = JX5$ASV_ID[which(JX5_ca$CA$v[,1] > 2)]
#asv_table[which(asv_table$ASV_ID %in% outlier_id), "JX5P1-2"] = 0

S2 = select(asv_table, ASV_ID, starts_with("S2"))
S2 = S2[rowSums(S2[,2:ncol(S2)])>0,]
S2_ca = cca(t(S2[,2:ncol(S2)])>0)
outlier_id = S2$ASV_ID[which(S2_ca$CA$v[,1] > 2)]
#asv_table[which(asv_table$ASV_ID %in% outlier_id), "S2P2-3"] = 0

# Filtering ASV table
#   drop out ASV with low abundance (<5)
#   drop samples at S5, S6, S7 in Argentina
#   drop two outliers JX5P1-2 and S2P2-3 (too many ASV singletons)
asv_table = asv_table %>%
    select(-starts_with("S5"), -starts_with("S6"), -starts_with("S7")) %>%
    #select(-`JX5P1-2`, -`S2P2-3`) %>%
    mutate(rs = rowSums(.[,2:ncol(.)]),
           rsi = rowSums(.[,2:ncol(.)]>0)) %>%
    filter(rs > 5) %>%
    select(-rs, -rsi)

write_csv(asv_table, "data/asv_table_renamed.csv")

# load taxonomy table and FUNGuild annotations
tax_table = read_csv("raw_data/taxonomy_table.csv")
tax_table = tax_table %>%
    select(ASV_ID = qseqid,
           taxid = staxid,
           kingdom:species) %>%
    filter(ASV_ID %in% asv_table$ASV_ID)

funguild_annotations = read_tsv("raw_data/taxonomy_table_for_FUNGuild.guilds_matched.txt", na=c("", "NA", "NULL"))
funguild_annotations = funguild_annotations %>% 
    select(ASV_ID=`Feature ID`,
           trophic_mode=`Trophic Mode`,
           guild=Guild,
           morphology=`Growth Morphology`,
           trait=Trait) %>%
    filter(ASV_ID %in% asv_table$ASV_ID)

asv_annotations = left_join(tax_table, funguild_annotations, by="ASV_ID")
write_lines(unique(asv_annotations$taxid), "data/taxids.txt")

write_csv(asv_annotations, "data/asv_annotations.csv")
