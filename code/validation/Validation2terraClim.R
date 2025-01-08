library(terra)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggpubr)

### First Validate Historical data by resampling topoterra to 4km

hist_years <- 1961:2022
hist_residual_maps <- vector("list", length(hist_years))
hist_residual_list <- foreach(i = hist_years, .packages = c("terra", "tidyverse")) %do% {
  # load topoterra
  topoterra <- rast(paste0("data/merged_output/topoterra_hist_", i, ".tif"))
  # convert tmax and tmin to C
  topoterra[["tmax"]] <- topoterra[["tmax"]]
  topoterra[["tmin"]] <- topoterra[["tmin"]]

  # load terraclim
  terraclim_file <- list.files("data/tile_templates", pattern = paste0("hist_", i, ".tif$"), full.names = TRUE)
  terraclim <- rast(terraclim_file)
  ##  clean terraclim like topoterra
  # AET
  terraclim[["aet"]][terraclim[["aet"]] < 0] <- 0
  terraclim[["aet"]][terraclim[["aet"]] > 1000] <- 1000
  terraclim[["aet"]] <- round(terraclim[["aet"]], 0)
  # DEF
  terraclim[["def"]][terraclim[["def"]] < 0] <- 0
  terraclim[["def"]][terraclim[["def"]] > 2500] <- 2500
  terraclim[["def"]] <- round(terraclim[["def"]], 0)
  # TMAX
  terraclim[["tmax"]][terraclim[["tmax"]] < -5] <- -5 # This happens to not be needed, but just for completeness
  terraclim[["tmax"]][terraclim[["tmax"]] > 50] <- 50
  terraclim[["tmax"]] <- round(terraclim[["tmax"]], 2)
  # TMIN
  terraclim[["tmin"]][terraclim[["tmin"]] < -25] <- -25
  terraclim[["tmin"]][terraclim[["tmin"]] > 20] <- 20
  terraclim[["tmin"]] <- round(terraclim[["tmin"]], 2)

  # resample topoterra to 4km
  topoterra_4km <- resample(topoterra, terraclim, threads = TRUE)

  # build topoterra vs terraclim residuals
  residuals <- terraclim - topoterra_4km

  hist_residual_maps[[i - 1960]] <- residuals
  # plot histograms of residuals using ggplot, facet by variable name
  df <- residuals %>%
    as.data.frame() %>%
    pivot_longer(cols = everything())

  residual_histograms <- df %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 100) +
    facet_wrap(~name, scales = "free") +
    ggtitle(paste0(i, " TopoTerra vs TerraClim Historical Residuals")) +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(strip.text = element_text(size = 14))

  final_list <- list(df = df, residual_histograms = residual_histograms)
}
# save list
save(hist_residual_list, file = "validation/hist_residual_list.RData")

# stack by variable, average then recombined
vars <- c("aet", "def", "tmax", "tmin")
hist_residual_maps_avg <- vector("list")
for (var in vars) {
  hist_residual_maps_var <- hist_residual_maps %>%
    map(~ .[[var]]) %>%
    rast() %>%
    mean()

  names(hist_residual_maps_var) <- var
  hist_residual_maps_avg[[var]] <- hist_residual_maps_var
}
hist_residual_maps_avg <- rast(hist_residual_maps_avg)

# plot average residuals histogram
hist_avg_residual_plot <- hist_residual_maps_avg %>%
  as.data.frame() %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~name, scales = "free") +
  ggtitle("Residual from TopoTerra to Terraclimate") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 14))
ggsave("validation/figures/hist_avg_residual_plot.jpg", hist_avg_residual_plot, width = 10, height = 10, units = "in")
writeRaster(hist_residual_maps_avg, "validation/figures/hist_residual_maps_avg.tif", overwrite = TRUE)

### now perform for 2c

c2_years <- 1985:2015
c2_residual_maps <- vector("list", length(c2_years))
c2_residual_list <- foreach(i = c2_years, .packages = c("terra", "tidyverse")) %do% {
  # load topoterra
  topoterra <- rast(paste0("data/merged_output/topoterra_2C_", i, ".tif"))
  # convert tmax and tmin to C
  topoterra[["tmax"]] <- topoterra[["tmax"]]
  topoterra[["tmin"]] <- topoterra[["tmin"]]

  # load terraclim
  terraclim_file <- list.files("data/tile_templates", pattern = paste0("2C_", i, ".tif$"), full.names = TRUE)
  terraclim <- rast(terraclim_file)
  ##  clean terraclim like topoterra
  # AET
  terraclim[["aet"]][terraclim[["aet"]] < 0] <- 0
  terraclim[["aet"]][terraclim[["aet"]] > 1000] <- 1000
  terraclim[["aet"]] <- round(terraclim[["aet"]], 0)
  # DEF
  terraclim[["def"]][terraclim[["def"]] < 0] <- 0
  terraclim[["def"]][terraclim[["def"]] > 2500] <- 2500
  terraclim[["def"]] <- round(terraclim[["def"]], 0)
  # TMAX
  terraclim[["tmax"]][terraclim[["tmax"]] < -5] <- -5 # This happens to not be needed, but just for completeness
  terraclim[["tmax"]][terraclim[["tmax"]] > 50] <- 50
  terraclim[["tmax"]] <- round(terraclim[["tmax"]], 2)
  # TMIN
  terraclim[["tmin"]][terraclim[["tmin"]] < -25] <- -25
  terraclim[["tmin"]][terraclim[["tmin"]] > 20] <- 20
  terraclim[["tmin"]] <- round(terraclim[["tmin"]], 2)

  # resample topoterra to 4km
  topoterra_4km <- resample(topoterra, terraclim, threads = TRUE)

  # build topoterra vs terraclim residuals
  residuals <- terraclim - topoterra_4km


  c2_residual_maps[[i - 1984]] <- residuals
  # plot histograms of residuals using ggplot, facet by variable name
  df <- residuals %>%
    as.data.frame() %>%
    pivot_longer(cols = everything())

  residual_histograms <- df %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 100) +
    facet_wrap(~name, scales = "free") +
    ggtitle(paste0(i, " TopoTerra vs TerraClim 2C Residuals")) +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(strip.text = element_text(size = 14))

  final_list <- list(df = df, residual_histograms = residual_histograms)
}

# save
save(c2_residual_list, file = "validation/2C_residual_list.RData")

# stack by variable, average then recombined
vars <- c("aet", "def", "tmax", "tmin")
c2_residual_maps_avg <- vector("list")
for (var in vars) {
  c2_residual_maps_var <- c2_residual_maps %>%
    map(~ .[[var]]) %>%
    rast() %>%
    mean()

  names(c2_residual_maps_var) <- var
  c2_residual_maps_avg[[var]] <- c2_residual_maps_var
}
c2_residual_maps_avg <- rast(c2_residual_maps_avg)

# plot average residuals histogram
c2_avg_residual_plot <- c2_residual_maps_avg %>%
  as.data.frame() %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~name, scales = "free") +
  ggtitle("Residual from TopoTerra to Terraclimate") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 14))
ggsave("validation/figures/2C_avg_residual_plot.jpg", c2_avg_residual_plot, width = 10, height = 10, units = "in")

writeRaster(c2_residual_maps_avg, "validation/figures/2C_residual_maps_avg.tif", overwrite = TRUE)

### now perform for 30 year normals- historical ----
# load topoterra normal
topoterra_normal <- rast("data/merged_output/topoterra_hist_1961-2022.tif")
# terraclim_normal
terraclim_normal <- rast("data/merged_output/terra_hist_1961-2022.tif")

# resample topoterra_normal to terraclim
topoterra_normal_4k <- resample(topoterra_normal, terraclim_normal, threads = T)

# build topoterra vs terraclim residuals
residuals_normal <- terraclim_normal - topoterra_normal_4k
# convert tmax and tmin to C
residuals_normal[["tmax"]] <- residuals_normal[["tmax"]]
residuals_normal[["tmin"]] <- residuals_normal[["tmin"]]



# write out
writeRaster(residuals_normal, "validation/figures/hist_normal_residuals.tif", overwrite = TRUE)

# plot histograms of residuals using ggplot, facet by variable name
df <- residuals_normal %>%
  as.data.frame() %>%
  pivot_longer(cols = everything())

residual_histograms <- df %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~name, scales = "free") +
  ggtitle("TopoTerra vs TerraClim Historical Normals Residuals") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 14))
ggsave("validation/figures/hist_normal_residuals.jpg", residual_histograms, width = 10, height = 10, units = "in")

historical_normal_list <- list(df = df, residual_histograms = residual_histograms)

### repeat for 2C ###
topoterra_normal <- rast("data/merged_output/topoterra_2C_1985-2015.tif")

# terraclim_normal cleaned to match topoterra

terraclim_normal <- rast("data/merged_output/terra_2C_1985-2015.tif")

# resample topoterra_normal to terraclim
topoterra_normal_4k <- resample(topoterra_normal, terraclim_normal, threads = T)

# build topoterra vs terraclim residuals
residuals_normal <- terraclim_normal - topoterra_normal_4k
# convert tmax and tmin to C
residuals_normal[["tmax"]] <- residuals_normal[["tmax"]]
residuals_normal[["tmin"]] <- residuals_normal[["tmin"]]

# write out
writeRaster(residuals_normal, "validation/figures/2C_normal_residuals.tif", overwrite = TRUE)

# plot histograms of residuals using ggplot, facet by variable name
df <- residuals_normal %>%
  as.data.frame() %>%
  pivot_longer(cols = everything())

residual_histograms <- df %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~name, scales = "free") +
  ggtitle("TopoTerra vs TerraClim 2C Normals Residuals") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 14))
ggsave("validation/figures/2C_normal_residuals.jpg", residual_histograms, width = 10, height = 10, units = "in")

c2_normal_list <- list(df = df, residual_histograms = residual_histograms)

# combine and save
normal_residual_list <- list(historical_normal_list = historical_normal_list, c2_normal_list = c2_normal_list)
save(normal_residual_list, file = "validation/normal_residual_list.RData")


### Major accuracy metrics
library(tidyverse)
library(viridisLite)
library(RColorBrewer)

# load residuals
load("validation/hist_residual_list.RData")
load("validation/2C_residual_list.RData")
load("validation/normal_residual_list.RData")

# calculate accuracy metrics

# historical
hist_residual_list_summary <- hist_residual_list %>%
  map(~ .[["df"]]) %>%
  bind_rows() %>%
  mutate(year = rep(1961:2022, each = nrow(hist_residual_list[[1]][["df"]]))) %>%
  # mutate(value = ifelse(name %in% c("tmax", "tmin"), value, value)) %>%
  group_by(name, year) %>%
  summarise(
    mean = mean(value, na.rm = T),
    median = median(value, na.rm = T),
    sd = sd(value, na.rm = T),
    min = min(value, na.rm = T),
    max = max(value, na.rm = T)
  ) %>%
  # rename to AET (mm), DEF (mm), TMax (°C) and Tmin (°C)
  mutate(name = case_when(
    name == "aet" ~ "AET (mm)",
    name == "def" ~ "DEF (mm)",
    name == "tmax" ~ "TMax (°C)",
    name == "tmin" ~ "Tmin (°C)"
  ))


# plot residuals over time by variable
historical_residuals_time <- hist_residual_list_summary %>%
  ggplot(aes(x = year, y = mean, color = name)) +
  facet_wrap(~name, scales = "free") +
  geom_line() +
  geom_ribbon(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd), alpha = 0.2) +
  ggtitle("Historical Residuals Over Time") +
  labs(x = "Year", y = "Residual") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    strip.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  ) +
  scale_color_brewer(palette = "Dark2")

# save
ggsave("validation/figures/historical_residuals_time.jpg", historical_residuals_time, width = 10, height = 10, units = "in")

# 2C

c2_residual_list_summary <- c2_residual_list %>%
  map(~ .[["df"]]) %>%
  bind_rows() %>%
  mutate(year = rep(1985:2015, each = nrow(c2_residual_list[[1]][["df"]]))) %>%
  mutate(value = ifelse(name %in% c("tmax", "tmin"), value, value)) %>%
  group_by(name, year) %>%
  summarise(
    mean = mean(value, na.rm = T),
    median = median(value, na.rm = T),
    sd = sd(value, na.rm = T),
    min = min(value, na.rm = T),
    max = max(value, na.rm = T)
  ) %>%
  # rename to AET (mm), DEF (mm), TMax (°C) and Tmin (°C)
  mutate(name = case_when(
    name == "aet" ~ "AET (mm)",
    name == "def" ~ "DEF (mm)",
    name == "tmax" ~ "TMax (°C)",
    name == "tmin" ~ "Tmin (°C)"
  ))
# plot residuals over time by variable

c2_residuals_time <- c2_residual_list_summary %>%
  ggplot(aes(x = year, y = mean, color = name)) +
  facet_wrap(~name, scales = "free") +
  geom_line() +
  geom_ribbon(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd), alpha = 0.2) +
  ggtitle("2C Residuals Over Time") +
  labs(x = "Year", y = "Residual") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    strip.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  ) +
  scale_color_brewer(palette = "Dark2")

# save

ggsave("validation/figures/2C_residuals_time.jpg", c2_residuals_time, width = 10, height = 10, units = "in")

# normal summary
normal_residual_list_summary <- normal_residual_list %>%
  map(~ .[["df"]]) %>%
  bind_rows() %>%
  mutate(dataset = rep(c("Historical", "2C"), each = nrow(normal_residual_list[[1]][["df"]]))) %>%
  group_by(name) %>%
  summarise(
    mean = mean(value, na.rm = T),
    median = median(value, na.rm = T),
    sd = sd(value, na.rm = T),
    min = min(value, na.rm = T),
    max = max(value, na.rm = T)
  ) %>%
  # rename to AET (mm), DEF (mm), TMax (°C) and Tmin (°C)
  mutate(name = case_when(
    name == "aet" ~ "AET (mm)",
    name == "def" ~ "DEF (mm)",
    name == "tmax" ~ "TMax (°C)",
    name == "tmin" ~ "Tmin (°C)"
  ))


# save summaries as one list
summary_list <- list(hist_residual_list_summary = hist_residual_list_summary, c2_residual_list_summary = c2_residual_list_summary, normal_residual_list_summary = normal_residual_list_summary)
save(summary_list, file = "validation/summary_list.RData")


# scatterplot of topoterra vs terraclim normals
topoterra_normal <- rast("data/merged_output/topoterra_hist_1961-2022.tif")
terraclim_normal <- rast("data/merged_output/terra_hist_1961-2022.tif")
# upsample topoterra to terraclim
topoterra_normal_4k <- resample(topoterra_normal, terraclim_normal, threads = T, method = "bilinear")

# plot
layer_names <- paste0(rep(c("aet", "def", "tmax", "tmin"), 2), rep(c("_topoterra", "_terraclim"), each = 4))

normal_scatter_df <- c(topoterra_normal_4k, terraclim_normal) %>%
  as.data.frame() %>%
  setNames(layer_names) %>%
  mutate(row_id = row_number()) %>%
  # build so that theres a terraclim and topoterra column for each variable, and a variable column listing the variable
  pivot_longer(cols = -row_id, names_to = "variable", values_to = "value") %>%
  separate(
    col = variable,
    into = c("variable", "source"),
    sep = "_"
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = value,
    id_cols = c(variable, row_id)
  ) %>%
  # change variable names to AET (mm), DEF (mm), TMax (°C) and Tmin (°C)
  mutate(variable = case_when(
    variable == "aet" ~ "AET (mm)",
    variable == "def" ~ "DEF (mm)",
    variable == "tmax" ~ "TMax (°C)",
    variable == "tmin" ~ "Tmin (°C)"
  )) %>%
  drop_na()
normal_scatter <- normal_scatter_df %>%
  ggplot(aes(x = terraclim, y = topoterra)) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  # add R^2adj equation
  stat_regline_equation(
    aes(label = paste(after_stat(adj.rr.label), sep = "~~~~")),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 5
  ) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("TopoTerra vs TerraClim Normals") +
  labs(x = "TerraClimate", y = "TopoTerra") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 14))
ggsave("validation/figures/topo_terra_hist_norm_scatter.jpg", normal_scatter, width = 10, height = 10, units = "in")

# repeat with 2C
topoterra_normal <- rast("data/merged_output/topoterra_2C_1985-2015.tif")
terraclim_normal <- rast("data/merged_output/terra_2C_1985-2015.tif")
# upsample topoterra to terraclim
topoterra_normal_4k <- resample(topoterra_normal, terraclim_normal, threads = T, method = "bilinear")

# plot
layer_names <- paste0(rep(c("aet", "def", "tmax", "tmin"), 2), rep(c("_topoterra", "_terraclim"), each = 4))

normal_scatter_df <- c(topoterra_normal_4k, terraclim_normal) %>%
  as.data.frame() %>%
  setNames(layer_names) %>%
  mutate(row_id = row_number()) %>%
  # build so that theres a terraclim and topoterra column for each variable, and a variable column listing the variable
  pivot_longer(cols = -row_id, names_to = "variable", values_to = "value") %>%
  separate(
    col = variable,
    into = c("variable", "source"),
    sep = "_"
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = value,
    id_cols = c(variable, row_id)
  ) %>%
  # change variable names to AET (mm), DEF (mm), TMax (°C) and Tmin (°C)
  mutate(variable = case_when(
    variable == "aet" ~ "AET (mm)",
    variable == "def" ~ "DEF (mm)",
    variable == "tmax" ~ "TMax (°C)",
    variable == "tmin" ~ "Tmin (°C)"
  )) %>%
  drop_na()
normal_scatter <- normal_scatter_df %>%
  ggplot(aes(x = terraclim, y = topoterra)) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  # add R^2adj equation
  stat_regline_equation(
    aes(label = paste(after_stat(adj.rr.label), sep = "~~~~")),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 5
  ) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("TopoTerra vs TerraClim 2C Normals") +
  labs(x = "TerraClimate", y = "TopoTerra") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 14))
ggsave("validation/figures/topo_terra_2C_norm_scatter.jpg", normal_scatter, width = 10, height = 10, units = "in")
