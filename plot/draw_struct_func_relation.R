library(tidyverse)
library(reshape2)
library(vegan)
library(cowplot)
library(SoDA)
library(ggrastr)

total_traits = read_csv("data/total_traits_final.csv")
asv_table_t = read.csv("data/asv_table_final.csv")
asv_annotations = read_csv("data/asv_annotations.csv")
asv_table = read_csv("data/asv_table_renamed.csv")
rownames(asv_table_t) = total_traits$sample_id

# remove outliers
drop = which(total_traits$sample_id %in% c("JX5P1-2", "S2P2-3"))
total_traits = total_traits[-drop,]
asv_table_t = asv_table_t[-drop,]
asv_table_t = asv_table_t[,colSums(asv_table_t)>0]
asv_table = select(asv_table, -all_of(c("JX5P1-2", "S2P2-3")))
asv_table = asv_table[rowSums(asv_table[,2:ncol(asv_table)]) > 0,]

group_by_troph = function(group_level) {
    aa = asv_table
    aa[2:ncol(aa)] = decostand(select(aa, -ASV_ID), method="pa", MARGIN=2)
    
    table_by_taxa = left_join(aa, asv_annotations, by="ASV_ID") %>%
        # grouping and summarizing
        group_by(across(all_of(group_level))) %>%
        summarise(across(matches("[AHGDJXS]{1,2}\\d{1,2}P\\d-\\d"), sum)) %>%
        ungroup()
    
    sample_order = total_traits %>%
        distinct(sample_id, lat, lon, origin) %>% 
        arrange(origin, lat, sample_id) %>%
        `$`(sample_id)
    
    table_by_taxa_ordered = table_by_taxa %>%
        mutate(rs = rowSums(.[,2:ncol(.)])) %>%
        arrange(desc(rs)) %>%
        select(-rs) %>%
        relocate(any_of(c(group_level,sample_order)))
    
    return(table_by_taxa_ordered)
}

# trophic guilds dissimilarity
tm_table = group_by_troph("trophic_mode")

tm_tags = tm_table$trophic_mode
tm_tags[which(is.na(tm_tags))] = "unkown"
tm_table = tm_table[2:ncol(tm_table)]
tm_table = t(tm_table)
colnames(tm_table) = tm_tags
tm_table = as.data.frame(tm_table)
tm_table = select(tm_table, -unkown)

dist_tro = as.matrix(vegdist(tm_table, method="jaccard"))

dist_tro[lower.tri(dist_tro, diag=TRUE)] = NA
dist_tro_df = dist_tro %>%
    melt() %>%
    as_tibble() %>%
    filter(!is.na(value)) 

# community dissimilarity
att = as.data.frame(asv_table_t)
rownames(att) = total_traits$sample_id
dist_com = vegdist(att, method = "jaccard", binary="T") %>%
    as.matrix()
dist_com[lower.tri(dist_com, diag=TRUE)] = NA
dist_com_df = dist_com %>%
    reshape2::melt() %>%
    as_tibble() %>%
    filter(!is.na(value)) %>%
    left_join(select(total_traits, Var1 = sample_id, Lat1 = lat, Lon1 = lon, Origin1=origin ,Cluster1=cluster)) %>%
    left_join(select(total_traits, Var2 = sample_id, Lat2 = lat, Lon2 = lon, Origin2=origin, Cluster2=cluster)) %>%
    mutate(Origin=Origin1==Origin2, Cluster=Cluster1==Cluster2)

dist_df = left_join(select(dist_tro_df, Var1, Var2, tro=value),
              select(dist_com_df, Var1, Var2, com=value)) %>%
    left_join(select(total_traits, Var1 = sample_id, Lat1 = lat, Lon1 = lon, Origin1=origin ,Cluster1=cluster)) %>%
    left_join(select(total_traits, Var2 = sample_id, Lat2 = lat, Lon2 = lon, Origin2=origin, Cluster2=cluster)) %>%
    mutate(Origin=Origin1==Origin2, Cluster=Cluster1==Cluster2) %>%
    mutate(spatial_scale = ifelse(Cluster, 1, ifelse(Origin, 2, 3)))

ggplot(aes(x=com, y=tro, color=as.factor(spatial_scale)), data=dist_df) +
    rasterise(geom_point(shape=1, alpha=0.8), dpi=300) + 
    xlab("Difference in community composition") +
    ylab("Difference in trophic composition") +
    theme_bw()

units = "cm"
width = 13
aspect_ratio = 0.6
ggsave("output/struct_func_relation.pdf", width=width, height=aspect_ratio*width, units=units)

# summaries trophic composition
ds = left_join(select(asv_annotations, ASV_ID, guild), asv_table) %>%
    pivot_longer(cols=3:ncol(.), names_to = "sample_id", values_to = "abundance") %>%
    left_join(select(total_traits, sample_id, origin, cluster)) %>%
    mutate(all=1) %>%
    filter(abundance > 0)

ds %>% distinct(ASV_ID, guild, all) %>% count(guild, all) %>% write_csv("output/trophic_composition_all.csv")
ds %>% distinct(ASV_ID, guild, origin) %>% count(guild, origin)%>% pivot_wider(names_from=origin, values_from=n) %>% write_csv("output/trophic_composition_origin.csv")
ds %>% distinct(ASV_ID, guild, cluster) %>% count(guild, cluster) %>% pivot_wider(names_from=cluster, values_from=n) %>% write_csv("output/trophic_composition_cluster.csv")
