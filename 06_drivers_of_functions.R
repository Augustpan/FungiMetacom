library(tidyverse)
library(vegan)
library(ggvegan)
library(ggsci)
library(cowplot)

# set random seed
set.seed(999)

# allow some vegan functions to utilize multi cores
options(mc.cores = 12)

# Load data
# To keep number of samples balanced between sites (15 samples per site),
#   we do not remove the two outliers (sample_id: JX5P1-2 and S2P2-3).
#   Instead, we use imputed data to replace them.
#   A balanced design is required for restricted permutation test in vegan.
#   Imputation was done in "01_detaclean.R"
total_traits = read_csv("data/total_traits_final.csv")
asv_table_t = read.csv("data/asv_table_final.csv")
rownames(asv_table_t) = total_traits$sample_id

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

if (T) {
    ind_nat = which(scaled_data$origin == "AR")
    scaled_data = scaled_data[ind_nat,]
    asv_table_t = asv_table_t[ind_nat,]
    asv_table_t = asv_table_t[,colSums(asv_table_t)>0]
}

envir = select(scaled_data, all_of(c(clim_predictor, soil_predictor)))
host = select(scaled_data, all_of(c(host_predictor)))
mems = select(scaled_data, starts_with("dbMEM"))
meta = select(scaled_data, origin, site, cluster, lat, lon)

centroids = select(scaled_data, centroid_x, centroid_y)
trophic_scores = select(scaled_data, patho_s, sapro_s, symbio_s)

mod = rda(trophic_scores ~ cpc2 + soil_NH + soil_SOM + dbMEM.1 + dbMEM.2, data=scaled_data)

# cn

#mod = rda(trophic_scores ~ cpc2 + cpc3 + cpc1 + soil_PH + SLA + leafbio + dbMEM.8 + dbMEM.6 + dbMEM.7 + dbMEM.11 + dbMEM.10 + dbMEM.9, data=scaled_data)

#mod = dbrda(vegdist(asv_table_t, method = "jaccard", binary = T) ~ cpc1 + cpc2 + cpc3 + soil_SOM + trichome + saponins +
#                SLA + leafbio + BI + SSL + dbMEM.6 + dbMEM.1 + dbMEM.2 + dbMEM.7 + dbMEM.8 + dbMEM.9, data=scaled_data)

st = scores(mod)$site
sp = scores(mod)$species
wa = scores(mod)$biplot * 1.6

my_pal = pal_d3("category10")(7)
names(my_pal) = paste0(1:7)

ggplot() +
    geom_point(aes(x=st[,1], y=st[,2], color=scaled_data$cluster)) +
    
    # Constrains
    geom_segment(aes(x=0, y=0, xend=wa[,1], yend=wa[,2]), color="red",
                 arrow=arrow(length=unit(0.15, "cm"))) +
    geom_text(aes(x=wa[,1],y=wa[,2],label=rownames(wa)),
              position='dodge') + 

    geom_vline(xintercept=0, linetype="dotted") + 
    geom_hline(yintercept=0, linetype="dotted") +
    
    xlab("RDA 1") + 
    ylab("RDA 2") +
    guides(color="none") + 
    scale_color_manual(values=my_pal) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size=1.5))

ggplot() +
    geom_point(aes(x=st[,1], y=st[,2], color=scaled_data$cluster)) +
    
    # Constrains
    geom_segment(aes(x=0, y=0, xend=wa[,1], yend=wa[,2]), color="red",
                 arrow=arrow(length=unit(0.15, "cm"))) +
    geom_text(aes(x=wa[,1],y=wa[,2],label=rownames(wa)),
              position='dodge') + 
    
    # Species Scores
    geom_segment(aes(x=0, y=0, xend=sp[,1], yend=sp[,2]), color="blue",
                 arrow=arrow(length=unit(0.15, "cm"))) +
    geom_text(aes(x=sp[,1],y=sp[,2],label=rownames(sp)),
              position='dodge') + 
    
    geom_vline(xintercept=0, linetype="dotted") + 
    geom_hline(yintercept=0, linetype="dotted") +
    xlab("RDA 1") + 
    ylab("RDA 2") +
    guides(color="none") + 
    scale_color_manual(values=my_pal) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size=1.5)) +
    ylim(c(-1.2,0.5)) # AR
