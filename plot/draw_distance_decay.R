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

# candidate explainary variables
host_predictor = c(
    "leafherb", "leafbio", 
    "SSL", "SLA", "BI",
    "trichome", "saponins", "lignin")
clim_predictor = c("cpc1", "cpc2", "cpc3")
soil_predictor = c(
    "soil_water", "soil_PH", 
    "soil_SOM", "soil_CN", 
    "soil_NH", "soil_NO", "soil_AP")
all_predictor = c(clim_predictor, host_predictor, soil_predictor)

# dataframes of predictor variables
scaled_data = mutate(total_traits, across(all_of(c(soil_predictor, host_predictor)), ~scale(.x)[,1])) %>%
    mutate(origin = as.factor(origin),
           site = as.factor(site),
           cluster = as.factor(cluster))
envir_mat = select(scaled_data, all_of(all_predictor))

# environment dissimilarity
dist_env = as.matrix(dist(envir_mat))
rownames(dist_env) = scaled_data$sample_id
colnames(dist_env) = scaled_data$sample_id
dist_env[lower.tri(dist_env, diag=TRUE)] = NA
dist_env_df = dist_env %>%
    melt() %>%
    as_tibble() %>%
    filter(!is.na(value)) %>%
    left_join(select(total_traits, Var1 = sample_id, Lat1 = lat, Lon1 = lon, Origin1=origin, Cluster1=cluster)) %>%
    left_join(select(total_traits, Var2 = sample_id, Lat2 = lat, Lon2 = lon, Origin2=origin, Cluster2=cluster)) %>%
    mutate(Origin=Origin1==Origin2, Cluster=Cluster1==Cluster2)

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

# geographic distance
XY1 = geoXY(dist_env_df$Lat1, dist_env_df$Lon1, 20, 58, unit=1000)
XY2 = geoXY(dist_env_df$Lat2, dist_env_df$Lon2, 20, 58, unit=1000)
dist_env_df = dist_env_df %>%
    mutate(X1 = XY1[,1], Y1=XY1[,2], X2=XY2[,1], Y2=XY2[,2]) %>%
    mutate(GDist = sqrt((X2-X1)^2+(Y2-Y1)^2)) %>%
    filter(Origin==TRUE, value != 0)

XY1 = geoXY(dist_com_df$Lat1, dist_com_df$Lon1, 20, 58, unit=1000)
XY2 = geoXY(dist_com_df$Lat2, dist_com_df$Lon2, 20, 58, unit=1000)
dist_com_df = dist_com_df %>%
    mutate(X1 = XY1[,1], Y1=XY1[,2], X2=XY2[,1], Y2=XY2[,2]) %>%
    mutate(GDist = sqrt((X2-X1)^2+(Y2-Y1)^2)) %>%
    filter(Origin==TRUE, value != 0)

# plots
g_env = ggplot(dist_env_df, aes(x=GDist, y=value, color=Cluster)) +
    rasterise(geom_point(alpha=.3), dpi = 300, scale = 1)+
    geom_smooth(method="lm") +
    xlab("Distance (km)") +
    ylab("Difference in environment") +
    theme_bw()

g_com = ggplot(dist_com_df, aes(x=GDist, y=value, color=Cluster)) +
    rasterise(geom_point(alpha=.3), dpi = 300, scale = 1)+
    geom_smooth(method="lm") +
    xlab("Distance (km)") +
    ylab("Difference in community composition") +
    theme_bw()


plot_grid(g_com, g_env, nrow=2, ncol=1, align = "hv", axis="lr")

units = "cm"
width = 11.4
aspect_ratio = 1.5
ggsave("output/distance_decay.pdf", width=width, height=aspect_ratio*width, units=units)