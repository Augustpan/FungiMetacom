library(tidyverse)
library(vegan)
library(mice)
library(SoDA)
library(adespatial)
library(pracma)

# set random seed
set.seed(999)

#### Define utility functions ####
create.dbmem = function (coord = NULL, D.mat = NULL, nsites, randtest=T) {
    if (is.null(coord) & is.null(D.mat)) 
        stop("Geographic information must be provided in 'coord' or in 'D.mat'")
    if (is.null(D.mat)) 
        D.mat <- dist(coord)
    D.mat <- as.matrix(D.mat)
    if (!is.null(coord)) {
        n <- nrow(coord)
    } else {
        n <- nrow(D.mat)
    }
    if (sum(nsites) != n) 
        stop("Vector nsites does not sum to nrow(coord) or nrow(D.mat)")
    if (min(nsites) == 1) 
        stop("At least one group contains a single site")
    out <- matrix(0, n, n)
    end <- 0
    end.mem <- 0
    for (k in 1:length(nsites)) {
        start <- end + 1
        end <- end + nsites[k]
        tmp <- as.dist(D.mat[start:end, start:end])
        res <- dbmem(tmp, MEM.autocor = "positive")
        if (randtest) {
            test.mem = moran.randtest(res, attributes(res)$listw)
            tb_res = as_tibble(res)
            sel_res = select(tb_res, test.mem$names[test.mem$pvalue<0.05])
        }
        dbMEM <- as.matrix(sel_res)
        n.mem <- ncol(dbMEM)
        out[start:end, (end.mem + 1):(end.mem + n.mem)] <- dbMEM
        end.mem <- end.mem + n.mem
    }
    out <- out[, 1:end.mem]
    if (is.null(rownames(coord))) {
        rownames(out) <- rownames(out, do.NULL = FALSE, prefix = "Site.")
    } else {
        rownames(out) <- rownames(coord)
    }
    colnames(out) <- colnames(out, do.NULL = FALSE, prefix = "dbMEM.")
    out
}
calc_trophic_scores = function(x) {
    
    patho_score = function (x){
        switch(x, 
               "Pathotroph"=1, 
               "Pathotroph-Symbiotroph"=1/2,
               "Pathotroph-Saprotroph"=1/2,
               "Pathotroph-Saprotroph-Symbiotroph"=1/3,
               0)
    }
    sapro_score = function (x){
        switch(x, 
               "Saprotroph"=1, 
               "Saprotroph-Symbiotroph"=1/2,
               "Pathotroph-Saprotroph"=1/2,
               "Pathotroph-Saprotroph-Symbiotroph"=1/3,
               0)
    }
    symbio_score = function (x){
        switch(x, 
               "Symbiotroph"=1, 
               "Pathotroph-Symbiotroph"=1/2,
               "Saprotroph-Symbiotroph"=1/2,
               "Pathotroph-Saprotroph-Symbiotroph"=1/3,
               0)
    }
    
    patho_score = Vectorize(patho_score)
    sapro_score = Vectorize(sapro_score)
    symbio_score = Vectorize(symbio_score)
    
    pa = patho_score(x)
    sa = sapro_score(x)
    sy = symbio_score(x)
    
    return(tibble(patho=pa, sapro=sa, symbio=sy))
}

#### Load raw data ####
plant_traits = read_csv("data/plant_traits.csv")
asv_table = read_csv("data/asv_table_renamed.csv")
asv_annotations = read_csv("data/asv_annotations.csv")
compartment_mapping = read_csv("data/compartment_mapping.csv")
compartment_mapping = distinct(compartment_mapping, site, compartment)

#### Alpha diversity ####
asv_table_t = t(select(asv_table,-ASV_ID))
fungal_diversity = tibble(
    sample_id = rownames(asv_table_t),
    richness = specnumber(asv_table_t),
    shannon = diversity(asv_table_t, "shannon"),
    evenness = shannon / log(richness)
)
total_traits = left_join(fungal_diversity, plant_traits, by="sample_id")

#### Clustering ####
run_optimize_clusters = FALSE
if (run_optimize_clusters) {
    cluster_func = function(x, k) {
        d = vegdist(x, method="jaccard", binary="T")
        return(hcut(d, k))
    }
    require(factoextra)
    fviz_nbclust(asv_table_t, cluster_func, method="gap_stat")
    ggsave("output/optimal_cluster_gap_stat.pdf", width=4.9, height=3.45, units="in")
    
    fviz_nbclust(asv_table_t, cluster_func, method="wss")
    ggsave("output/optimal_cluster_wss.pdf", width=4.9, height=3.45, units="in")
}

optimal_num_clusters = 7
clusters = asv_table_t %>%
    vegdist(method = "jaccard", binary = T) %>%
    hclust(method = "ward.D2") %>%
    cutree(k = optimal_num_clusters)
total_traits$cluster = clusters

#### Trophic composition ####
tm = left_join(asv_table, distinct(asv_annotations, ASV_ID, trophic_mode), by="ASV_ID")
ts = calc_trophic_scores(tm$trophic_mode)
tms = tm %>% mutate(across(2:(ncol(tm)-1), ~ifelse(.x>0, 1, 0))) %>%group_by(trophic_mode) %>% summarise(across(2:(ncol(tm)-1), sum))
ttm = tms$trophic_mode
ttm[which(is.na(ttm))] = "Unknowntroph"
tdf = t(tms[,2:ncol(tms)])
sid = rownames(tdf)
tdf = as_tibble(tdf)
colnames(tdf) = ttm
tdf$sample_id = sid
total_traits = left_join(total_traits, tdf)
asv_table_pa = select(asv_table, -ASV_ID) %>% decostand(method="pa", MARGIN=2)
pa = (asv_table_pa * ts$patho) %>% 
    as_tibble() %>%
    colSums()
sa = (asv_table_pa * ts$sapro) %>% 
    as_tibble() %>%
    colSums()
sy = (asv_table_pa * ts$symbio) %>% 
    as_tibble() %>%
    colSums()

trophic_table = tibble(
    sample_id = colnames(asv_table_pa),
    patho = pa,
    sapro = sa,
    symbio = sy
)

x = trophic_table[,2:4] %>% scale()
center = sweep(x, 2, apply(x, 2, min),'-')
R = apply(x, 2, max) - apply(x, 2, min)
x_star = sweep(center, 2, R, "/") %>% as.data.frame()

theta = (c(1,2,3)-1) * 2*pi/3
xm = sweep(x_star, 2, cos(theta), "*")
ym = sweep(x_star, 2, sin(theta), "*")

centroids = t(apply(cbind(xm, ym), 1, function(x) poly_center(x[1:3], x[4:6])))
area = apply(cbind(xm, ym), 1, function(x) polyarea(x[1:3], x[4:6]))

trophic_table = trophic_table %>%
    mutate(patho_s = x_star[["patho"]], 
           sapro_s = x_star[["sapro"]],
           symbio_s = x_star[["symbio"]],
           centroid_x = centroids[,1],
           centroid_y = centroids[,2],
           area = area) 
write_csv(trophic_table, "data/trophic_table.csv")

#### Assign compartment ####

total_traits = total_traits %>%
    left_join(compartment_mapping, by="site") %>%
    left_join(trophic_table, by="sample_id")

#### Calculate climate PCs ####
bioclim_vars = select(total_traits, V1:V19)
clim_pc = prcomp(scale(bioclim_vars))
total_traits$cpc1 = clim_pc$x[,1]
total_traits$cpc2 = clim_pc$x[,2]
total_traits$cpc3 = clim_pc$x[,3]

#### Generate dbMEM ####
ordered = arrange(total_traits, origin)
cartesian_coord = geoXY(ordered$lat, ordered$lon, unit=1000) # unit=km
num_nat = sum(total_traits$origin == "AR")
num_int = sum(total_traits$origin == "CN")
dbMEMs = create.dbmem(coord=cartesian_coord, nsites=c(num_nat, num_int), randtest=T)
dbMEMs.nat = colnames(dbMEMs[,which(dbMEMs[1,]!=0)])
dbMEMs.int = colnames(dbMEMs[,which(dbMEMs[nrow(dbMEMs),]!=0)])
total_traits = as_tibble(cbind(ordered, dbMEMs))

#### Data imputation ####
total_traits = arrange(total_traits, sample_id)
asv_table_t = t(select(asv_table, -ASV_ID)) %>%
    as.data.frame() %>%
    mutate(sample_id = rownames(.)) %>%
    arrange(sample_id) %>%
    select(-sample_id)

host_predictor = c(
    "leafherb", "leafbio", 
    "SSL", "SLA", "BI",
    "trichome", "saponins", "lignin")
clim_predictor = c("cpc1", "cpc2", "cpc3")
soil_predictor = c(
    "soil_water", "soil_PH", 
    "soil_SOM", "soil_CN", 
    "soil_NH", "soil_NO", "soil_AP")
spatial_predictor = c(dbMEMs.nat, dbMEMs.int)
all_predictor = c(clim_predictor, host_predictor, soil_predictor)

d = select(total_traits, all_of(all_predictor))
temp_data = mice(d, m=5, maxit=50, meth='pmm', seed=999, printFlag=F)
completed_data = complete(temp_data, 1)
total_traits[,all_predictor] = completed_data

#### Gamma richness of each cluster ####
gamma_richness = NULL
for (i in 1:7) {
    rowidx = which(total_traits$cluster == i)
    gamma_richness[i] = sum(colSums(asv_table_t[rowidx,]) > 0)
}
nsamples = count(total_traits, cluster, name="cluster_n")

total_traits = total_traits %>%
    left_join(tibble(cluster=1:7, gamma_richness=gamma_richness), by="cluster") %>%
    left_join(nsamples, by="cluster")

write_csv(total_traits, "data/total_traits_final.csv")
write_csv(asv_table_t, "data/asv_table_final.csv")