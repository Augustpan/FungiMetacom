library(tidyverse)
library(ggsci)

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

#### Taxonomic composition ####

group_by_taxa = function(group_level, plot=T, max_legend=10, by_cluster=T) {
    
    table_by_taxa = left_join(asv_table, asv_annotations, by="ASV_ID") %>%
        # grouping and summarizing
        group_by(across(all_of(group_level))) %>%
        summarise(across(matches("[AHGDJXS]{1,2}\\d{1,2}P\\d-\\d"), sum)) %>%
        ungroup()
    
    table_by_taxa[,2:ncol(table_by_taxa)] = decostand(select(table_by_taxa, -all_of(group_level)), method="total", MARGIN=2)
    
    sample_order = total_traits %>%
        distinct(sample_id, lat, lon, origin) %>% 
        arrange(origin, lat, sample_id) %>%
        `$`(sample_id)
    
    table_by_taxa_ordered = table_by_taxa %>%
        mutate(rs = rowSums(.[,2:ncol(.)])) %>%
        arrange(desc(rs)) %>%
        select(-rs) %>%
        relocate(any_of(c(group_level,sample_order)))
    
    if (nrow(table_by_taxa_ordered) > max_legend) {
        table_by_taxa_ordered[max_legend:nrow(table_by_taxa_ordered),1] = "Others"
        table_by_taxa_ordered = table_by_taxa_ordered %>%
            group_by(across(all_of(group_level))) %>%
            summarise(across(matches("[AHGDJXS]{1,2}\\d{1,2}P\\d-\\d"), sum)) %>%
            ungroup()
    }
    #write_csv(table_by_taxa_ordered, paste0("output/table_by_", group_level, ".csv"))
    
    taxa_order = table_by_taxa_ordered %>% 
        filter(if_all(group_level, ~.x != "Others")) %>%
        mutate(rs = rowSums(.[,2:ncol(.)])) %>%
        arrange(desc(rs)) %>%
        select(-rs) %>%
        `[[`(group_level)
    
    taxa_order = c(taxa_order, "Others")
    
    if (plot) {
        # pivot_longer for ggplot drawing
        table_by_taxa_long = pivot_longer(table_by_taxa_ordered, 
                                          cols = matches("[AHGDJXS]{1,2}\\d{1,2}P\\d-\\d"), 
                                          names_to = "sample_id", 
                                          values_to = "abundance") %>%
            left_join(distinct(total_traits, sample_id, origin, site, cluster), by="sample_id")
        
        fill_colors = c("#000000", "#B15928",
                        "#33A02C", "#B2DF8A", 
                        "#1F78B4", "#FB9A99", 
                        "#FDBF6F", "#E3A1AC", 
                        "#6A3D9A", "#CAB2D6",
                        "#FF7F00", "#FFFF99", 
                        "#6ACE13")
        if (by_cluster) {
            g = ggplot(table_by_taxa_long, 
                       aes(x = factor(sample_id, levels=sample_order), 
                           y = abundance, 
                           fill = factor(.data[[group_level]], levels=taxa_order))) + 
                geom_bar(stat="identity", width=1) + 
                facet_wrap(~factor(cluster, levels=c(3,1,2,7,4,5,6)), scales="free_x") +
                xlab("Sample") +
                ylab("Relative abundance") +
                labs(fill=group_level) +
                theme(axis.ticks.x = element_blank(),
                      axis.text.x = element_blank(),
                      panel.grid = element_blank(),
                      panel.background = element_blank()) + 
                scale_fill_manual(values=rev(fill_colors))
            ggsave(filename=paste0("output/composition_", group_level, "_by_cluster" , ".pdf"), 
                   plot=g, width=7, height=3.3, units="in")
        } else {
            g = ggplot(table_by_taxa_long, 
                       aes(x = factor(sample_id, levels=sample_order), 
                           y = abundance, 
                           fill = factor(.data[[group_level]], levels=taxa_order))) + 
                geom_bar(stat="identity", width=1) + 
                facet_wrap(~factor(origin, levels=c("AR", "CN")), scales="free_x") +
                xlab("Sample") +
                ylab("Relative abundance") +
                labs(fill=group_level) +
                theme(axis.ticks.x = element_blank(),
                      axis.text.x = element_blank(),
                      panel.grid = element_blank(),
                      panel.background = element_blank()) + 
                scale_fill_manual(values=rev(fill_colors))
            ggsave(filename=paste0("output/composition_", group_level, "_by_origin" , ".pdf"), 
                   plot=g, width=7, height=3.3, units="in")
        }
        
        
    }
    
    return(table_by_taxa_ordered)
}

group_by_taxa("phylum", plot=T, max_legend=20, T)
group_by_taxa("class", plot=T, max_legend=11, T)
group_by_taxa("trophic_mode", plot=T, max_legend=10, T)

group_by_taxa("phylum", plot=T, max_legend=20, F)
group_by_taxa("class", plot=T, max_legend=11, F)
group_by_taxa("trophic_mode", plot=T, max_legend=10, F)