# This script compare the 3 methods performed by gids_glm_downscaling.R
library(tidyr)
library(ggplot2)

terraclim <- rast('data/cogs/def_hist_1981_2010_cog.tif') 
names(terraclim) <- 'terra_val'

tile <- 'grand'

tile_list <- list.files('output/def/', pattern = paste0('buff.*_',tile,'.tiff$'), full.names = T)
#grand_list <- list.files('output/def/', pattern = 'buff.*_grand.tiff$', full.names = T)

out_stack <- rast(tile_list) 
names(out_stack) <- gsub(paste0('_',tile,'.tiff$'),'',gsub('output/def/buff_','',tile_list))

sample_pts <- spatSample(out_stack, size = 200000, as.points = T)

terra_sample <- terra::extract(terraclim, sample_pts)

tile_sample <- terra::extract(out_stack, sample_pts) |> 
  pivot_longer(-ID, names_to = 'temp', values_to = 'pred_val') |> 
  separate(temp, into = c('buff','method'), sep = '_', extra = 'merge') |> 
  left_join(terra_sample, by = 'ID')

# Calculate fit
fit_fn <- function(data){
  # Fit the linear model
  lm_mod <- lm(terra_val ~ pred_val, data = data)
  # Save the results as a data.frame
  result <- data.frame('slope' = lm_mod$coefficients['pred_val'],
                       'r2' = round(summary(lm_mod)$r.squared, 3),
                       'intercept' = round(lm_mod$coefficients['(Intercept)']))
  # Return the data.frame
  return(result)
}

model_df <- tile_sample %>% 
  group_by(buff, method) %>% 
  group_modify(~ fit_fn(.x))

# Plot fits
ggplot(slice_sample(tile_sample, n = 5000)) +
  geom_point(aes(x = pred_val, y = terra_val), alpha = 0.5) +
  geom_abline(data = model_df, 
              aes(intercept = intercept, slope = slope), linewidth = 0.75, color = 'red') +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.75, color = 'blue') +
  geom_label(data = model_df,
             aes(x = 1700, y = 800, label = r2)) + #x = 100, y = 450  
  facet_grid(buff~method) +
  theme_bw(base_size = 14) +
  labs(x = 'Predicted 250-m value', y = 'TerraClimate 4-km at fine centroid')

ggsave(paste0('fits_',tile,'.png'), width = 8, height = 6, dpi = 300)

# HISTORGRAMS
ggplot(tile_sample) +
  geom_histogram(aes(x = pred_val), fill = 'black') +
  facet_grid(buff~method) +
  theme_bw(base_size = 14) +
  labs(x = 'Predicted 250-m value', y = 'Count')

ggsave(paste0('hist_',tile,'.png'), width = 8, height = 6, dpi = 300)







