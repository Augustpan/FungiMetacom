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

envir = select(scaled_data, all_of(c(clim_predictor, soil_predictor)))
host = select(scaled_data, all_of(c(host_predictor)))
mems = select(scaled_data, starts_with("dbMEM"))
meta = select(scaled_data, origin, site, cluster, lat, lon)

if (T) {
    perm = how(nperm = 1000, 
               within = Within(mirror=F, constant=F, type="free"), 
               plots = Plots(strata = meta$site, mirror=F, type="free"))
    trophic_scores = select(scaled_data, ends_with("troph"))
    mod0 = rda(trophic_scores ~ 1, data=envir)
    mod1 = rda(trophic_scores ~ ., data=envir)
    mods = ordiR2step(mod0, 
                      scope = formula(mod1), 
                      direction = "forward",
                      #permutations = perm, 
                      trace = FALSE)
    envir_vars = attr(mods$terms, "term.labels")
    
    mod0 = rda(trophic_scores ~ 1, data=host)
    mod1 = rda(trophic_scores ~ ., data=host)
    mods = ordiR2step(mod0, 
                      scope = formula(mod1), 
                      direction = "forward",
                      #permutations = perm, 
                      trace = FALSE)
    host_vars = attr(mods$terms, "term.labels")
    
    mod0 = rda(trophic_scores ~ 1, data=mems)
    mod1 = rda(trophic_scores ~ ., data=mems)
    mods = ordiR2step(mod0, 
                      scope = formula(mod1), 
                      direction = "forward",
                      #permutations = perm, 
                      trace = FALSE)
    space_vars = attr(mods$terms, "term.labels")
    
    mod_trop_all = rda(trophic_scores ~ ., select(scaled_data, all_of(c(envir_vars, host_vars, space_vars))))
}

if (T) {
    ind_nat = which(scaled_data$origin == "CN")
    scaled_data = scaled_data[ind_nat,]
    asv_table_t = asv_table_t[ind_nat,]
    asv_table_t = asv_table_t[,colSums(asv_table_t)>0]
    envir = select(scaled_data, all_of(c(clim_predictor, soil_predictor)))
    host = select(scaled_data, all_of(c(host_predictor)))
    mems = select(scaled_data, starts_with("dbMEM"))
    meta = select(scaled_data, origin, site, cluster, lat, lon)
    
    perm = how(nperm = 1000, 
               within = Within(mirror=F, constant=F, type="free"), 
               plots = Plots(strata = meta$site, mirror=F, type="free"))
    trophic_scores = select(scaled_data, ends_with("troph"))
    mod0 = rda(trophic_scores ~ 1, data=envir)
    mod1 = rda(trophic_scores ~ ., data=envir)
    mods = ordiR2step(mod0, 
                      scope = formula(mod1), 
                      direction = "forward",
                      #permutations = perm, 
                      trace = FALSE)
    envir_vars = attr(mods$terms, "term.labels")
    
    mod0 = rda(trophic_scores ~ 1, data=host)
    mod1 = rda(trophic_scores ~ ., data=host)
    mods = ordiR2step(mod0, 
                      scope = formula(mod1), 
                      direction = "forward",
                      #permutations = perm, 
                      trace = FALSE)
    host_vars = attr(mods$terms, "term.labels")
    
    mod0 = rda(trophic_scores ~ 1, data=mems)
    mod1 = rda(trophic_scores ~ ., data=mems)
    mods = ordiR2step(mod0, 
                      scope = formula(mod1), 
                      direction = "forward",
                      #permutations = perm, 
                      trace = FALSE)
    space_vars = attr(mods$terms, "term.labels")
    
    mod_trop_cn = rda(trophic_scores ~ ., select(scaled_data, all_of(c(envir_vars, host_vars, space_vars))))
    
    my_pal = pal_d3("category10")(7)
    names(my_pal) = paste0(1:7)
    st = scores(mod_trop_cn)$sites
    sp = scores(mod_trop_cn)$species*2
    wa = scores(mod_trop_cn)$biplot * 10
    ggplot() +
        geom_point(aes(x=st[,1], y=st[,2], color=as.factor(scaled_data$cluster))) + 
        geom_segment(aes(x=0, y=0, xend=wa[,1], yend=wa[,2]), arrow=arrow(length = unit(0.3, "cm")), color="red", alpha=0.7) + 
        geom_segment(aes(x=0, y=0, xend=sp[,1], yend=sp[,2]), arrow=arrow(length = unit(0.3, "cm")), color="blue", alpha=0.7) + 
        geom_text(aes(x=sp[,1], y=sp[,2], label=rownames(sp))) + 
        geom_text(aes(x=wa[,1], y=wa[,2], label=rownames(wa))) + 
        scale_color_manual(values=my_pal)
}

if (T) {
    ind_nat = which(scaled_data$origin == "AR")
    scaled_data = scaled_data[ind_nat,]
    asv_table_t = asv_table_t[ind_nat,]
    asv_table_t = asv_table_t[,colSums(asv_table_t)>0]
    envir = select(scaled_data, all_of(c(clim_predictor, soil_predictor)))
    host = select(scaled_data, all_of(c(host_predictor)))
    mems = select(scaled_data, starts_with("dbMEM"))
    meta = select(scaled_data, origin, site, cluster, lat, lon)
    
    perm = how(nperm = 1000, 
               within = Within(mirror=F, constant=F, type="free"), 
               plots = Plots(strata = meta$site, mirror=F, type="free"))
    trophic_scores = select(scaled_data, ends_with("troph"), -Unknowntroph)
    mod0 = rda(trophic_scores ~ 1, data=envir)
    mod1 = rda(trophic_scores ~ ., data=envir)
    mods = ordiR2step(mod0, 
                      scope = formula(mod1), 
                      direction = "forward",
                      #permutations = perm, 
                      trace = FALSE)
    envir_vars = attr(mods$terms, "term.labels")
    
    mod0 = rda(trophic_scores ~ 1, data=host)
    mod1 = rda(trophic_scores ~ ., data=host)
    mods = ordiR2step(mod0, 
                      scope = formula(mod1), 
                      direction = "forward",
                      #permutations = perm, 
                      trace = FALSE)
    host_vars = attr(mods$terms, "term.labels")
    
    mod0 = rda(trophic_scores ~ 1, data=mems)
    mod1 = rda(trophic_scores ~ ., data=mems)
    mods = ordiR2step(mod0, 
                      scope = formula(mod1), 
                      direction = "forward",
                      #permutations = perm, 
                      trace = FALSE)
    space_vars = attr(mods$terms, "term.labels")
    
    mod_trop_ar = rda(trophic_scores ~ ., select(scaled_data, all_of(c(envir_vars, host_vars, space_vars))))
    
    my_pal = pal_d3("category10")(7)
    names(my_pal) = paste0(1:7)
    st = scores(mod_trop_ar)$sites
    sp = scores(mod_trop_ar)$species*2
    wa = scores(mod_trop_ar)$biplot * 10
    ggplot() +
        geom_point(aes(x=st[,1], y=st[,2], color=as.factor(scaled_data$cluster))) + 
        geom_segment(aes(x=0, y=0, xend=wa[,1], yend=wa[,2]), arrow=arrow(length = unit(0.3, "cm")), color="red", alpha=0.7) + 
        geom_segment(aes(x=0, y=0, xend=sp[,1], yend=sp[,2]), arrow=arrow(length = unit(0.3, "cm")), color="blue", alpha=0.7) + 
        geom_text(aes(x=sp[,1], y=sp[,2], label=rownames(sp))) + 
        geom_text(aes(x=wa[,1], y=wa[,2], label=rownames(wa))) + 
        ylim(c(-10,10))+
        scale_color_manual(values=my_pal)
}

load_saved = F
if (load_saved == T) {
    load("output/selected_models.RData")
} else {
    # overall
    if (T) {
        perm = how(nperm = 1000, 
                   within = Within(mirror=F, constant=F, type="none"), 
                   plots = Plots(strata = meta$site, mirror=F, type="free"))
        
        mod0 = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ 1, data=envir)
        mod1 = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ ., data=envir, distance="jaccard", binary=T)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          permutations = perm, 
                          trace = FALSE)
        envir_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ 1, data=host)
        mod1 = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ ., data=host)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          permutations = perm, 
                          trace = FALSE)
        host_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ 1, data=mems)
        mod1 = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ ., data=mems)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          permutations = perm, 
                          trace = FALSE)
        space_vars = attr(mods$terms, "term.labels")
        
        mems_sel = select(scaled_data, all_of(space_vars)) %>% as.matrix()
        
        mod_all = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ .,
                        data=select(scaled_data, all_of(c(envir_vars, host_vars, space_vars))))
    }
    
    # plot of overall
    if (T) {
        st1 = scores(mod_all)$sites[,1]
        st2 = scores(mod_all)$sites[,2]
        
        wscale= 1.85
        wa1 = scores(mod_all)$biplot[,1]*wscale
        wa2 = scores(mod_all)$biplot[,2]*wscale
        sp1 = scores(mod_all)$species[,1]
        sp2 = scores(mod_all)$species[,2]
        
        ggplot() + 
            geom_point(aes(x=st1, y=st2, color=as.factor(cluster)), alpha=0.85, data=scaled_data) +
            annotate("segment", x=0,y=0,xend=wa1,yend=wa2,
                     arrow=arrow(length = unit(0.3, "cm"))) +
            geom_text(aes(x=wa1,y=wa2,label=names(wa1))) +
            theme_bw() + 
            scale_color_d3()
        
        units = "cm"
        width = 11.4
        aspect_ratio = 0.68
        ggsave("output/dbrda_mod_all.pdf", width=width, height=aspect_ratio*width, units=units)
    }
    
    mod_origin = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ origin, data=scaled_data)
    mod_cluster = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ cluster + Condition(origin), data=scaled_data)
    mod_site = dbrda(vegdist(asv_table_t, method="jaccard", binary=T) ~ site + Condition(cluster), data=scaled_data)
    
    fit_board_scale_model = T
    if (fit_board_scale_model) {
        perm = how(nperm = 1000, 
                   within = Within(mirror=F, constant=F, type="none"), 
                   plots = Plots(strata = meta$site, mirror=F, type="free"))
        
        mod0 = dbrda(fitted(mod_origin) ~ 1, data=envir)
        mod1 = dbrda(fitted(mod_origin) ~ ., data=envir)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          permutations = perm, 
                          trace = FALSE)
        envir_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(fitted(mod_origin) ~ 1, data=host)
        mod1 = dbrda(fitted(mod_origin) ~ ., data=host)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          permutations = perm, 
                          trace = FALSE)
        host_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(fitted(mod_origin) ~ 1, data=mems)
        mod1 = dbrda(fitted(mod_origin) ~ ., data=mems)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          permutations = perm, 
                          trace = FALSE)
        space_vars = attr(mods$terms, "term.labels")
        
        mod_bs = dbrda(fitted(mod_origin) ~ ., data=select(scaled_data, all_of(c(envir_vars, host_vars, space_vars))))
    }
    
    fit_mid_scale_model = T
    if (fit_mid_scale_model) {
        perm = how(nperm = 1000, 
                   within = Within(mirror=F, constant=F, type="none"), 
                   plots = Plots(strata = meta$site, mirror=F, type="free"))
        
        mod0 = dbrda(fitted(mod_cluster) ~ 1, data=envir)
        mod1 = dbrda(fitted(mod_cluster) ~ ., data=envir)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          permutations = perm, 
                          trace = FALSE)
        envir_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(fitted(mod_cluster) ~ 1, data=host)
        mod1 = dbrda(fitted(mod_cluster) ~ ., data=host)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          permutations = perm, 
                          trace = FALSE)
        host_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(fitted(mod_cluster) ~ 1, data=mems)
        mod1 = dbrda(fitted(mod_cluster) ~ ., data=mems)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          permutations = perm, 
                          trace = FALSE)
        space_vars = attr(mods$terms, "term.labels")
        
        vp = varpart(fitted(mod_cluster),
                     select(scaled_data, all_of(envir_vars)) %>% as.matrix(),
                     select(scaled_data, all_of(host_vars)) %>% as.matrix(),
                     select(scaled_data, all_of(space_vars)) %>% as.matrix())
        tiff("output/varpart_cluster.tif", height=400, width=500)
        plot(vp, Xnames=c("Envir.", "Host", "Space"), bg=2:4)
        title(paste0("between clusters", " Total R2adj = ", RsquareAdj(mod_cluster)$adj.r.squared))
        dev.off()
        
        mod_ms = dbrda(fitted(mod_cluster) ~ ., data=select(scaled_data, all_of(c(envir_vars, host_vars, space_vars))))
        write_lines(c("Envir:", envir_vars, "Host:", host_vars, "Space:", space_vars),
                    "selected_board_scale.txt")
    }
    
    fit_fine_scale_model = T
    if (fit_fine_scale_model) {
        perm = how(nperm = 1000, 
                   blocks = meta$cluster)
        
        mod0 = dbrda(fitted(mod_site) ~ 1, data=envir)
        mod1 = dbrda(fitted(mod_site) ~ ., data=envir)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          trace = FALSE)
        envir_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(fitted(mod_site) ~ 1, data=host)
        mod1 = dbrda(fitted(mod_site) ~ ., data=host)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          trace = FALSE)
        host_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(fitted(mod_site) ~ 1, data=mems)
        mod1 = dbrda(fitted(mod_site) ~ ., data=mems)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          trace = FALSE)
        space_vars = attr(mods$terms, "term.labels")
        
        vp = varpart(fitted(mod_site),
                     select(scaled_data, all_of(envir_vars)) %>% as.matrix(),
                     select(scaled_data, all_of(host_vars)) %>% as.matrix(),
                     select(scaled_data, all_of(space_vars)) %>% as.matrix())
        
        tiff("output/varpart_within_cluster.tif", height=400, width=500)
        plot(vp, Xnames=c("Envir.", "Host", "Space"), bg=2:4)
        title(paste0("within clusters", " Total R2adj = ", RsquareAdj(mod_site)$adj.r.squared))
        dev.off()
        
        mod_fs = dbrda(fitted(mod_site) ~ ., data=select(scaled_data, all_of(c(envir_vars, host_vars, space_vars))))
        write_lines(c("Envir:", envir_vars, "Host:", host_vars, "Space:", space_vars),
                    "selected_mid_scale.txt")
    }
    
    fit_resid_scale_model = TRUE
    if (fit_resid_scale_model) {
        perm = how(nperm = 1000, 
                   blocks = meta$cluster)
        
        mod0 = dbrda(resid(mod_cluster) ~ 1, data=envir)
        mod1 = dbrda(resid(mod_cluster) ~ ., data=envir)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          trace = FALSE)
        envir_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(resid(mod_cluster) ~ 1, data=host)
        mod1 = dbrda(resid(mod_cluster) ~ ., data=host)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          trace = FALSE)
        host_vars = attr(mods$terms, "term.labels")
        
        mod0 = dbrda(resid(mod_cluster) ~ 1, data=mems)
        mod1 = dbrda(resid(mod_cluster) ~ ., data=mems)
        mods = ordiR2step(mod0, 
                          scope = formula(mod1), 
                          direction = "forward",
                          trace = FALSE)
        space_vars = attr(mods$terms, "term.labels")
        
        vp = varpart(resid(mod_cluster),
                     select(scaled_data, all_of(envir_vars)) %>% as.matrix(),
                     select(scaled_data, all_of(host_vars)) %>% as.matrix(),
                     select(scaled_data, all_of(space_vars)) %>% as.matrix())
        
        tiff("output/varpart_within_cluster.tif", height=400, width=500)
        plot(vp, Xnames=c("Envir.", "Host", "Space"), bg=2:4)
        title(paste0("within clusters", " Total R2adj = ", RsquareAdj(mod_site)$adj.r.squared))
        dev.off()
        
        mod_re = dbrda(resid(mod_cluster) ~ ., data=select(scaled_data, all_of(c(envir_vars, host_vars, space_vars))))
        write_lines(c("Envir:", envir_vars, "Host:", host_vars, "Space:", space_vars),
                    "selected_mid_scale.txt")
    }
    
    save(mod_all, mod_origin, mod_cluster, mod_site, mod_bs, mod_ms, mod_fs, mod_re,
         file="output/selected_models.RData")
}



