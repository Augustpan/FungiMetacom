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

#### Alpha diversity: Richness ####
g_richness = ggplot() +
    geom_point(aes(x=lat, y=richness, color=origin), data=total_traits) +
    stat_smooth(aes(x=lat, y=richness), data=filter(total_traits, origin=="AR"), 
                method="lm", formula=y~x+I(x^2), color="blue") +
    stat_smooth(aes(x=lat, y=richness), data=filter(total_traits, origin=="CN"), 
                method="lm", formula=y~x, color="red") +
    xlab("Latitude") +
    ylab("Fungal richness") +
    scale_color_manual(values=c("AR"="blue", "CN"="red")) +
    theme_bw()

g_richness_ar = ggplot() +
    geom_point(aes(x=lat, y=richness, color=origin), data=filter(total_traits, origin=="AR")) +
    stat_smooth(aes(x=lat, y=richness), data=filter(total_traits, origin=="AR"), 
                method="lm", formula=y~x+I(x^2), color="blue") +
    xlab("Latitude") +
    ylab("Fungal richness") +
    scale_color_manual(values=c("AR"="blue", "CN"="red")) +
    theme_bw()

g_richness_cn = ggplot() +
    geom_point(aes(x=lat, y=richness, color=origin), data=filter(total_traits, origin=="CN")) +
    stat_smooth(aes(x=lat, y=richness), data=filter(total_traits, origin=="CN"), 
                method="lm", formula=y~x, color="red") +
    xlab("Latitude") +
    ylab("Fungal richness") +
    scale_color_manual(values=c("AR"="blue", "CN"="red")) +
    theme_bw()

width = 11.4
aspect_ratio = 2.2
plot_grid(g_richness, g_richness_ar, g_richness_cn, nrow=3, ncol=1, align = "hv")
ggsave("output/richness_lat.pdf", width=width, height=aspect_ratio*width, units="cm")

#### Alpha diversity: Shannon and Evenness ####
g_shannon = ggplot(total_traits, aes(x=lat, y=shannon, color=origin)) +
    geom_point() +
    geom_smooth(method="gam", formula=y~s(x, k=5, bs="cr")) +
    xlab("Latitude") +
    ylab("Shannon index") +
    scale_color_manual(values=c("AR"="blue", "CN"="red")) +
    theme_bw()

g_evenness = ggplot(total_traits, aes(x=lat, y=evenness, color=origin)) +
    geom_point() +
    geom_smooth(method="gam", formula=y~s(x, k=5, bs="cr")) +
    xlab("Latitude") +
    ylab("Pielou's evenness") +
    scale_color_manual(values=c("AR"="blue", "CN"="red")) +
    theme_bw()

units = "cm"
width = 11.4
aspect_ratio = 2.2
plot_grid(g_richness, g_shannon, g_evenness, nrow=3, ncol=1, align = "hv")
ggsave("output/alpha_diversity_lat.pdf", width=width, height=aspect_ratio*width, units=units)

#### TODO: family, genus richness ####

# byfamily = group_by_taxa("family", F, 1000)
# byfamily = select(byfamily, -family) %>% t()
# total_traits$richness_family = rowSums(byfamily>0)

# byfamily = group_by_taxa("genus", F, 1000)
# byfamily = select(byfamily, -genus) %>% t()
# total_traits$richness_genus = rowSums(byfamily>0)

# byfamily = group_by_taxa("order", F, 1000)
# byfamily = select(byfamily, -order) %>% t()
# total_traits$richness_order = rowSums(byfamily>0)

#### Beta diversity, by compartment ####
d = betadiver(asv_table_t, "w") %>%
    as.matrix() %>%
    as.data.frame()
d[lower.tri(d, diag=TRUE)] = NA
d$R = rownames(d) 
dd = d %>% 
    pivot_longer(cols=1:(ncol(.)-1), names_to="C") %>% 
    na.omit() %>%
    left_join(distinct(total_traits,R=sample_id,cluster_R=cluster)) %>%
    left_join(distinct(total_traits,C=sample_id,cluster_C=cluster)) %>%
    mutate(cluster = ifelse(cluster_R==cluster_C, cluster_R, NA)) %>%
    na.omit() %>%
    left_join(distinct(total_traits,cluster,origin))

ggplot(aes(x=factor(cluster, levels=c(3,1,2,7,4,5,6)), y=value, fill=as.factor(cluster)), data=dd) +
    geom_boxplot(alpha=0.8) + 
    scale_fill_d3() + 
    theme_bw()

ggsave("output/beta_diversity.pdf", width=11.4, height=8, units="cm")