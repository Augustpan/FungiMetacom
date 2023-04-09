library(tidyverse)
library(vegan)
library(cowplot)
library(ggsci)
library(pracma)

total_traits = read_csv("data/total_traits_final.csv")
asv_table_t = read_csv("data/asv_table_final.csv")

# remove outliers
drop = which(total_traits$sample_id %in% c("JX5P1-2", "S2P2-3"))
total_traits = total_traits[-drop,]
asv_table_t = asv_table_t[-drop,]
asv_table_t = asv_table_t[,colSums(asv_table_t)>0]

qplot(x=factor(cluster,levels=c(3,1,2,7,4,5,6)), 
      y=Pathotroph, 
      fill=factor(cluster,levels=c(3,1,2,7,4,5,6)), 
      data=total_traits, geom="boxplot") +
    scale_fill_d3()

qplot(x=factor(cluster,levels=c(3,1,2,7,4,5,6)), 
      y=Saprotroph, 
      fill=factor(cluster,levels=c(3,1,2,7,4,5,6)), 
      data=total_traits, geom="boxplot") +
    scale_fill_d3()

smm = total_traits %>% 
    group_by(cluster) %>% 
    summarise(pa = mean(patho_s), 
              sa = mean(sapro_s), 
              sy = mean(symbio_s), 
              cx = mean(centroid_x), 
              cy = mean(centroid_y), 
              ar = mean(area)) %>%
    mutate(pax = pa*cos(0),
           pay = pa*sin(0),
           sax = sa*cos(2*pi/3),
           say = sa*sin(2*pi/3),
           syx = sy*cos(4*pi/3),
           syy = sy*sin(4*pi/3))

glist = list()
for (i in 1:7) {
    axis_scale = 0.5
    glist[[i]] = ggplot(data=smm[i,]) +
        annotate("polygon", 
                 x=c(smm[i,]$pax,smm[i,]$sax,smm[i,]$syx), 
                 y=c(smm[i,]$pay,smm[i,]$say,smm[i,]$syy), 
                 fill=my_pal[i], alpha=.8) + 
        geom_point(aes(x=cx,y=cy)) + 
        geom_segment(aes(x=0,y=0,xend=cos(0) * axis_scale, yend=sin(0) * axis_scale), 
                     arrow = arrow(length = unit(0.3, "cm"))) +
        geom_segment(aes(x=0,y=0,xend=cos(2*pi/3) * axis_scale, yend=sin(2*pi/3) * axis_scale),
                     arrow = arrow(length = unit(0.3, "cm"))) +
        geom_segment(aes(x=0,y=0,xend=cos(4*pi/3) * axis_scale, yend=sin(4*pi/3) * axis_scale),
                     arrow = arrow(length = unit(0.3, "cm"))) +
        annotate("path",
                 x=0.4*cos(seq(0,2*pi,length.out=100)),
                 y=0.4*sin(seq(0,2*pi,length.out=100)),
                 linetype = "dashed",
                 size = 0.2) +
        annotate("path",
                 x=0.3*cos(seq(0,2*pi,length.out=100)),
                 y=0.3*sin(seq(0,2*pi,length.out=100)),
                 linetype = "dashed",
                 size = 0.2) +
        annotate("path",
                 x=0.2*cos(seq(0,2*pi,length.out=100)),
                 y=0.2*sin(seq(0,2*pi,length.out=100)),
                 linetype = "dashed",
                 size = 0.2) + 
        annotate("path",
                 x=0.1*cos(seq(0,2*pi,length.out=100)),
                 y=0.1*sin(seq(0,2*pi,length.out=100)),
                 linetype = "dashed",
                 size = 0.2) +
        annotate("text",
                 x = cos(0) * axis_scale - 0.08,
                 y = sin(0) * axis_scale - 0.08,
                 label = "Pa",
                 hjust = "left") +
        annotate("text",
                 x = cos(2*pi/3) * axis_scale + 0.04,
                 y = sin(2*pi/3) * axis_scale + 0.03,
                 label = "Sa",
                 hjust = "left") +
        annotate("text",
                 x = cos(4*pi/3) * axis_scale + 0.04,
                 y = sin(4*pi/3) * axis_scale - 0.03,
                 label = "Sy",
                 hjust = "left") +
        annotate("path",
                 x=0.55*cos(seq(0,2*pi,length.out=100)),
                 y=0.55*sin(seq(0,2*pi,length.out=100)),
                 color=NA) +
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              panel.background = element_blank())
    ggsave(paste0("output/trophic_triangles_", i, ".pdf"), width=4.5, height=4.5, units="cm")
}

plot_grid(glist[[1]],glist[[2]],glist[[3]],glist[[4]],glist[[5]],glist[[6]],glist[[7]])
ggsave("output/trophic_triangles.pdf", width=10, height=10, units="in")

qplot(x=sqrt(area), y=richness, data=total_traits) + 
    geom_smooth(method="lm") + 
    ggpubr::stat_cor()
ggsave("output/triangle_area_vs_richness.pdf", width=5, height=5, units="in")


gsa = qplot(x=lat, y=sapro_s, data=total_traits) + facet_wrap(~origin) + geom_smooth()
gpa = qplot(x=lat, y=patho_s, data=total_traits) + facet_wrap(~origin) + geom_smooth()
gsy = qplot(x=lat, y=symbio_s, data=total_traits) + facet_wrap(~origin) + geom_smooth()
plot_grid(gsa, gpa, gsy, nrow=3)
ggsave("output/trophic_scores_lat.pdf", width=5, height=8, units="in")

gcx = qplot(x=lat, y=centroid_x, data=total_traits) + facet_wrap(~origin) + geom_smooth()
gcy = qplot(x=lat, y=centroid_y, data=total_traits) + facet_wrap(~origin) + geom_smooth()
plot_grid(gcx, gcy, nrow=2)
ggsave("output/trophic_centroids_lat.pdf", width=5, height=5.5, units="in")

axis_scale = 0.1
qplot(x=centroid_x, y=centroid_y, color=lat, data=total_traits) + facet_wrap(~origin) + 
    geom_segment(aes(x=0,y=0,xend=cos(0) * axis_scale, yend=sin(0) * axis_scale), 
                 arrow = arrow(length = unit(0.3, "cm")), color="black") +
    geom_segment(aes(x=0,y=0,xend=cos(2*pi/3) * axis_scale, yend=sin(2*pi/3) * axis_scale),
                 arrow = arrow(length = unit(0.3, "cm")), color="black") +
    geom_segment(aes(x=0,y=0,xend=cos(4*pi/3) * axis_scale, yend=sin(4*pi/3) * axis_scale),
                 arrow = arrow(length = unit(0.3, "cm")), color="black")

axis_scale = 0.1
ggplot(smm) + 
    geom_segment(aes(x=0, y=0, xend=cx, yend=cy, color=as.factor(cluster)), 
                 arrow=arrow(length = unit(0.3, "cm"))) + 
    geom_segment(aes(x=0,y=0,xend=cos(0) * axis_scale, yend=sin(0) * axis_scale), 
                 arrow = arrow(length = unit(0.3, "cm")), color="black") +
    geom_segment(aes(x=0,y=0,xend=cos(2*pi/3) * axis_scale, yend=sin(2*pi/3) * axis_scale),
                 arrow = arrow(length = unit(0.3, "cm")), color="black") +
    geom_segment(aes(x=0,y=0,xend=cos(4*pi/3) * axis_scale, yend=sin(4*pi/3) * axis_scale),
                 arrow = arrow(length = unit(0.3, "cm")), color="black") +
    scale_color_d3()

qplot(x=centroid_x, y=centroid_y, color=lat, shape=as.factor(cluster), data=total_traits) + stat_ellipse() + scale_color_gsea() + facet_wrap(~origin) + geom_hline(yintercept=0, size=0.2) + geom_vline(xintercept=0, size=0.2) + 
    geom_segment(aes(x=0,y=0,xend=cos(0) * axis_scale, yend=sin(0) * axis_scale), 
                 arrow = arrow(length = unit(0.3, "cm")), color="black", size=0.01) +
    geom_segment(aes(x=0,y=0,xend=cos(2*pi/3) * axis_scale, yend=sin(2*pi/3) * axis_scale),
                 arrow = arrow(length = unit(0.3, "cm")), color="black", size=0.01) +
    geom_segment(aes(x=0,y=0,xend=cos(4*pi/3) * axis_scale, yend=sin(4*pi/3) * axis_scale),
                 arrow = arrow(length = unit(0.3, "cm")), color="black", size=0.01)

g_sapro = ggplot() +
    geom_point(aes(x=lat, y=Saprotroph, color=origin), data=total_traits) +
    stat_smooth(aes(x=lat, y=Saprotroph), data=filter(total_traits, origin=="AR"), 
                method="lm", formula=y~x+I(x^2), color="blue") +
    stat_smooth(aes(x=lat, y=Saprotroph), data=filter(total_traits, origin=="CN"), 
                method="lm", formula=y~x, color="red") +
    scale_color_manual(values=c("AR"="blue", "CN"="red")) +
    theme_bw()

g_patho = ggplot() +
    geom_point(aes(x=lat, y=Pathotroph, color=origin), data=total_traits) +
    stat_smooth(aes(x=lat, y=Pathotroph), data=filter(total_traits, origin=="AR"), 
                method="lm", formula=y~x+I(x^2), color="blue") +
    stat_smooth(aes(x=lat, y=Pathotroph), data=filter(total_traits, origin=="CN"), 
                method="lm", formula=y~x, color="red") +
    scale_color_manual(values=c("AR"="blue", "CN"="red")) +
    theme_bw()

g_symbio = ggplot() +
    geom_point(aes(x=lat, y=Symbiotroph, color=origin), data=total_traits) +
    stat_smooth(aes(x=lat, y=Symbiotroph), data=filter(total_traits, origin=="AR"), 
                method="lm", formula=y~x+I(x^2), color="blue") +
    stat_smooth(aes(x=lat, y=Symbiotroph), data=filter(total_traits, origin=="CN"), 
                method="lm", formula=y~x+I(x^2), color="red") +
    scale_color_manual(values=c("AR"="blue", "CN"="red")) +
    theme_bw()

width = 11.4
aspect_ratio = 2.2
plot_grid(g_sapro, g_patho, g_symbio, nrow=3, ncol=1, align = "hv")
ggsave("output/troph_lat.pdf", width=width, height=aspect_ratio*width, units="cm")