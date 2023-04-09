library(rnaturalearth)
library(scatterpie)
library(tidyverse)
library(ggspatial)
library(cowplot)
library(ggsci)
library(sf)

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

d = total_traits %>%
    count(site, lat, lon, origin, cluster) %>% 
    pivot_wider(names_from = "cluster", values_from="n", values_fill=0)

my_pal = pal_d3("category10")(7)
names(my_pal) = paste0(1:7)

# draw China
country_list = c('china', 'philippines')
sf_tw = ne_countries(country = 'taiwan', returnclass="sf", scale="large")

sf_states = ne_states(country = country_list, returnclass="sf")
sf_countries = ne_countries(country = country_list, returnclass="sf", scale="large")

lat = total_traits$lat
lon = total_traits$lon
site = total_traits$site
region = as.factor(total_traits$cluster)
latlon = data.frame(lon=lon, lat=lat, site, region)
sf_site = st_as_sf(latlon, coords = c("lon", "lat"), crs = 4326)

fg_fill = "#e1e1e1"
bg_fill = "#ffffff"

chn = ggplot() +
    geom_sf(data = sf_states, aes(fill=iso_a2), size=0.1) +
    geom_sf(data = sf_countries, aes(fill=iso_a2), size=0.3, alpha=0) +
    geom_sf(data = sf_tw, fill=fg_fill, size=0.3) +
    #geom_sf(data = sf_site, aes(color=region), colour="black",pch=21, size=3.5) +
    geom_scatterpie(aes(x=lon, y=lat, group=site, r=0.35), cols=paste(1:7), data=d, alpha=1) + 
    #annotation_scale(location = "br", width_hint = 0.4) +
    coord_sf(xlim = c(110, 120), ylim = c(20, 37), label_axes = list(bottom="E", right = "N")) +
    theme_bw() + 
    theme(panel.grid = element_blank(),axis.text=element_text(size=12,face = "bold")) + 
    scale_fill_manual(values = c("CN" = fg_fill, "TW"=fg_fill, "PN"=bg_fill, my_pal)) + 
    scale_color_manual(values = my_pal) +
    scale_x_continuous(breaks = c(111, 115, 119)) +
    guides(fill="none", color="none")

# draw Argentina
country_list = c('argentina', 'uruguay', 'brazil', 'paraguay', 'bolivia')

sf_states = ne_states(country = country_list, returnclass="sf")
sf_countries = ne_countries(country = country_list, returnclass="sf", scale="large")

lat = -total_traits$lat
lon = -total_traits$lon
site = total_traits$site
region = as.factor(total_traits$cluster)
latlon = data.frame(lon=lon, lat=lat, site, region)
sf_site = st_as_sf(latlon, coords = c("lon", "lat"), crs = 4326)

fg_fill = "#e1e1e1"
bg_fill = "#ffffff"

arg = ggplot() +
    geom_sf(data = sf_states, aes(fill=iso_a2), size=0.1) +
    geom_sf(data = sf_countries, aes(fill=iso_a2), size=0.3, alpha=0) +
    #geom_sf(data = sf_site, aes(color=region), size = 3.5) +
    geom_scatterpie(aes(x=-lon, y=-lat, group=site, r=0.35), cols=paste(1:7), data=d, alpha=1) + 
    #annotation_scale(location = "br", width_hint = 0.5) +
    coord_sf(xlim = c(-64, -54), ylim = c(-40, -23)) +
    theme_bw() + 
    theme(panel.grid = element_blank(), axis.text=element_text(size=12,face = "bold")) + 
    scale_fill_manual(values = c("AR" = fg_fill, "UY"=bg_fill, "BR"=bg_fill, "PY"=bg_fill, "BO"=bg_fill, my_pal)) + 
    scale_color_manual(values = my_pal) +
    scale_x_continuous(breaks = c(-63, -59, -55)) +
    guides(fill="none", color="none")

map_merged = plot_grid(arg, chn, align = "hv", axis="lr")

ggsave(plot = map_merged, 
       filename = "output/map.pdf", 
       height = 8.25, 
       units = "cm")