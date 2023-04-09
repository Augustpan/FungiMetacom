library(tidyverse)

total_traits = read_csv("data/total_traits_final.csv")

tmp = total_traits %>%
  select(origin, lat, cpc1, cpc2, cpc3) %>%
  pivot_longer(cpc1:cpc3, names_to="CPC_names", values_to="CPC_values")

ggplot(data=tmp) +
  geom_point(aes(x=lat, y=CPC_values, color=origin)) +
  facet_grid(~CPC_names) +
  ylab("Climate PCs") +
  xlab("latitude") +
  theme_bw()

tmp = total_traits %>%
  select(origin, lat, V1:V19) %>%
  pivot_longer(V1:V19, names_to="Bioclim_var_names", values_to="Bioclim_var_values")

ggplot(data=tmp) +
  geom_point(aes(x=lat, y=Bioclim_var_values, color=origin)) +
  facet_wrap(~Bioclim_var_names, scales="free_y") +
  ylab("Bioclimatic variables") +
  xlab("latitude") +
  theme_bw()

tmp = total_traits %>%
  select(origin, lat, starts_with("soil_")) %>%
  pivot_longer(starts_with("soil_"), names_to="var_names", values_to="var_values")

ggplot(data=tmp) +
  geom_point(aes(x=lat, y=var_values, color=origin)) +
  facet_wrap(~var_names, scales="free_y") +
  ylab("Soil properties") +
  xlab("latitude") +
  theme_bw()

tmp = total_traits %>%
  select(origin, lat, leafherb:CN) %>%
  pivot_longer(leafherb:CN, names_to="var_names", values_to="var_values")

ggplot(data=tmp) +
  geom_point(aes(x=lat, y=var_values, color=origin)) +
  facet_wrap(~var_names, scales="free_y") +
  ylab("Plant traits") +
  xlab("latitude") +
  theme_bw()
