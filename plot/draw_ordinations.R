library(tidyverse)
library(cowplot)
library(ggsci)
library(vegan)

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

asv_table_t_pa = decostand(asv_table_t, method="pa")

# Check length of the first axis
decorana(asv_table_t_pa) # DCA1 > 4, use unimodal method
decorana(asv_table_t)    # DCA1 < 3, use linear method

ca_pa = cca(asv_table_t_pa)
pcoa = cmdscale(vegdist(asv_table_t, method="jaccard", binary=T), eig=T)
#pcoa = cmdscale(vegdist(asv_table_t, method="bray"), eig=T)

ord_table = tibble(sample_id = rownames(asv_table_t),
                   CA1 = ca_pa$CA$u[,1],
                   CA2 = ca_pa$CA$u[,2],
                   PCoA1 = pcoa$points[,1],
                   PCoA2 = pcoa$points[,2]) %>%
    left_join(distinct(total_traits, sample_id, site, origin, lat, cluster), 
              by="sample_id")

units = "cm"
width = 11.4
aspect_ratio = 0.68

# Main plot
ggplot(ord_table, aes(x=PCoA1, y=PCoA2, color=as.factor(cluster))) +
    geom_point() +
    xlab(paste0("PCoA1 ", round(pcoa$eig/sum(pcoa$eig)*100, 1)[1], "%")) +
    ylab(paste0("PCoA2 ", round(pcoa$eig/sum(pcoa$eig)*100, 1)[2], "%")) +
    theme_bw() + 
    scale_color_d3()
ggsave("output/PCoA_pa_cluster.pdf", width=width, height=aspect_ratio*width, units=units)

# Auxiliary plots
ggplot(ord_table, aes(x=CA1, y=CA2, color=as.factor(cluster))) +
    geom_point(alpha=.75) +
    theme_bw() + 
    scale_color_d3()
ggsave("output/CA_pa_cluster.pdf", width=width, height=aspect_ratio*width, units=units)

ggplot(ord_table, aes(x=CA1, y=CA2, color=lat, shape=origin)) +
    geom_point() +
    scale_color_gsea() +
    ggtitle("CA with presence-absence data")
ggsave("output/CA_pa_lat.pdf", width=6.15, height=4.5)

ggplot(ord_table, aes(x=PCoA1, y=PCoA2, color=lat, shape=origin)) +
    geom_point() +
    ggtitle("PCoA analysis with Jaccard distance") +
    scale_color_gsea() +
    xlab(paste0("PCoA1 ", round(pcoa$eig/sum(pcoa$eig)*100, 1)[1], "%")) +
    ylab(paste0("PCoA2 ", round(pcoa$eig/sum(pcoa$eig)*100, 1)[2], "%"))
ggsave("output/PCoA_lat.pdf", width=6.15, height=4.5)