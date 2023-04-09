library(tidyverse)
library(vegan)
library(ggsci)
library(ComplexHeatmap)

toggle_EMS_test = F

asv_table_t = read.csv("data/asv_table_final.csv")
total_traits = read_csv("data/total_traits_final.csv")
rownames(asv_table_t) = total_traits$sample_id

ra = function(x) {
    ca = cca(x > 0)
    st = ca$CA$u[,1] # site score
    sp = ca$CA$v[,1] # species score
    st_ord = order(st) # site order
    sp_ord = order(sp) # species order
    x_ordered = (x[st_ord, sp_ord] > 0)*1
    return(x_ordered)
}


asv_table_no_singletons = asv_table_t[,colSums(asv_table_t > 0)>1]
asv_table_no_singletons_ordered = ra(asv_table_no_singletons)
asv_table_ordered = ra(asv_table_t)
write_csv(as_tibble(asv_table_ordered), "output/asv_table_ordered.csv")
write_csv(as_tibble(asv_table_no_singletons_ordered), "output/asv_table_no_singletons_ordered.csv")

# run EMS tests
#   EMS analysis on full ASV table may take several hours
if (toggle_EMS_test) {
    system2("/usr/local/bin/julia", "run_EMS_analysis.jl", "output/asv_table_ordered.csv")
    system2("/usr/local/bin/julia", "run_EMS_analysis.jl", "output/asv_table_no_singletons_ordered.csv")
}


cCN = asv_table_no_singletons_ordered[39:165,]
cAR = asv_table_no_singletons_ordered[211:333,]
cCN = cCN[,colSums(cCN)>0]
cAR = cAR[,colSums(cAR)>0]
cpCN = ra(cCN)
cpAR = ra(cAR)

my_pal = pal_d3("category10")(7)
names(my_pal) = paste0(1:7)

A = asv_table_no_singletons_ordered
B = cpCN
C = cpAR
metainfo = left_join(
    tibble(sample_id = rownames(A)),
    select(total_traits, sample_id, origin, cluster, site, lat),
    by="sample_id"
) %>% mutate(cluster = as.factor(cluster))

l_anno = rowAnnotation(
    cluster = as.factor(metainfo$cluster),
    col = list(origin = c("AR"="blue","CN"="red"), cluster = my_pal)
)

A[1:38,] %>% as_tibble() %>% write_csv("output/asv_table_no_singletons_ordered_a1.csv")
A[39:165,] %>% as_tibble() %>% write_csv("output/asv_table_no_singletons_ordered_a2.csv")
A[166:210,] %>% as_tibble() %>% write_csv("output/asv_table_no_singletons_ordered_a3.csv")
A[211:333,] %>% as_tibble() %>% write_csv("output/asv_table_no_singletons_ordered_a4.csv")
A[334:nrow(A),] %>% as_tibble() %>% write_csv("output/asv_table_no_singletons_ordered_a5.csv")
B[1:42,] %>% as_tibble() %>% write_csv("output/asv_table_no_singletons_ordered_b1.csv")
B[42:nrow(B),] %>% as_tibble() %>% write_csv("output/asv_table_no_singletons_ordered_b2.csv")
C[1:63,] %>% as_tibble() %>% write_csv("output/asv_table_no_singletons_ordered_c1.csv")
C[64:nrow(C),] %>% as_tibble() %>% write_csv("output/asv_table_no_singletons_ordered_c2.csv")


# compartment annotations
row_sp = 1:nrow(A)
row_sp[1:38] = "1"
row_sp[39:165] = "2"
row_sp[166:210] = "3"
row_sp[211:333] = "4"
row_sp[334:length(row_sp)] = "5"

row_sp = 1:nrow(A)
row_sp[1:210] = "1"
row_sp[210:length(row_sp)] = "2"
tiff("output/compartment_A.tif", width=800, height=800)
heatmap_A = A %>%
    apply(2, as.factor) %>%
    as.matrix() %>%
    Heatmap(left_annotation=l_anno, 
            col = c("white", "black"),
            row_split = row_sp,
            row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), 
            border = TRUE,
            border_gp = gpar(col = "black", lwd=5),
            cluster_columns=F, 
            cluster_rows=F, 
            show_row_names = F,
            column_labels = rep("", ncol(.)),
            color_space="RGB",
            use_raster=T,
            raster_quality = 10)
draw(heatmap_A)
dev.off()

metainfo = left_join(
    tibble(sample_id = rownames(B)),
    select(total_traits, sample_id, origin, cluster, site, lat),
    by="sample_id"
) %>% mutate(cluster = as.factor(cluster))

l_anno = rowAnnotation(
    cluster = as.factor(metainfo$cluster),
    col = list(origin = c("AR"="blue","CN"="red"), cluster = my_pal)
)

# compartment annotations
row_sp = 1:nrow(B)
row_sp[1:42] = "2"
row_sp[42:length(row_sp)] = "1"

tiff("output/compartment_B.tif", width=800, height=800)
write_csv(as.data.frame(B), "B.csv")
heatmap_B = B %>%
    apply(2, as.factor) %>%
    as.matrix() %>%
    Heatmap(left_annotation=l_anno, 
            col = c("white", "black"),
            #row_split = row_sp,
            row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE,
            border_gp = gpar(col = "black", lwd=5),
            cluster_columns=F, 
            cluster_rows=F, 
            show_row_names = F,
            column_labels = rep("", ncol(.)),
            color_space="RGB",
            use_raster=T,
            raster_quality = 10)
draw(heatmap_B)
dev.off()

metainfo = left_join(
    tibble(sample_id = rownames(C)),
    select(total_traits, sample_id, origin, cluster, site, lat),
    by="sample_id"
) %>% mutate(cluster = as.factor(cluster))

l_anno = rowAnnotation(
    cluster = as.factor(metainfo$cluster),
    col = list(origin = c("AR"="blue","CN"="red"), cluster = my_pal)
)

# compartment annotations
row_sp = 1:nrow(C)
row_sp[1:63] = "2"
row_sp[64:length(row_sp)] = "1"

tiff("output/compartment_C.tif", width=800, height=800)
write_csv(as.data.frame(C[2:nrow(C),]), "C.csv")
heatmap_C = C %>%
    apply(2, as.factor) %>%
    as.matrix() %>%
    Heatmap(left_annotation=l_anno, 
            col = c("white", "black"),
            #row_split = row_sp,
            row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), 
            border_gp = gpar(col = "black", lwd=5),
            border = TRUE,
            cluster_columns=F, 
            cluster_rows=F, 
            show_row_names = F,
            column_labels = rep("", ncol(.)),
            color_space="RGB",
            use_raster=F)
draw(heatmap_C)
dev.off()
