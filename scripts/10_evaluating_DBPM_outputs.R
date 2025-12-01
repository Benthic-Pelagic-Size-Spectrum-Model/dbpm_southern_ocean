# Evaluating DBPM estimated catches against observations using the observation
# range adjusted (ORA) method

# Loading libraries -------------------------------------------------------
library(dplyr)
library(ggplot2)
library(arrow)
library(tidyr)
library(cowplot)
library(stringr)

# Define function to load catches data ------------------------------------
get_catch_dbpm_obs <- function(base_dir){
  catch_dbpm_nonspatial <- list.files(
    file.path(base_dir, "run_nonspatial/1deg"), full.names = T) |> 
    read_parquet() |> 
    select(c(time, year, total_catch)) |> 
    group_by(year) |> 
    summarise(catch_calib_org = mean(total_catch, na.rm = T)) 
  
  catch_dbpm_1deg <- list.files(
    file.path(base_dir, "gridded_dbpm_outputs/1deg"), 
    pattern = "mean_year_catch_dbpm_1", full.names = T) |> 
    read_parquet(col_select = !units) |> 
    rename(catch_1_org = mean_catch)
  
  catch_si_dbpm_1deg <- list.files(
    file.path(base_dir, "gridded_dbpm_outputs/1deg"),
    pattern = "mean_year_catch_dbpm_simask", full.names = T) |> 
    read_parquet(col_select = !units) |> 
    rename(catch_1_si = mean_catch)
  
  catch_si_dbpm_025deg <- list.files(
    file.path(base_dir, "gridded_dbpm_outputs/025deg"),
    pattern = "mean_year_catch_dbpm_simask", full.names = T) |> 
    read_parquet(col_select = !units) |> 
    rename(catch_025_si = mean_catch)
  
  catch_dbpm <- list.files(
    file.path(base_dir, "gridded_dbpm_outputs/025deg"),
    pattern = "mean_year_catch_dbpm_025", full.names = T) |> 
    read_parquet(col_select = !units) |> 
    rename(catch_025_org = mean_catch) |>
    left_join(catch_dbpm_1deg, by = "year") |>
    left_join(catch_si_dbpm_1deg, by = "year") |>
    left_join(catch_si_dbpm_025deg, by = "year") |>
    left_join(catch_dbpm_nonspatial, by = "year") |>
    filter(year >= 1961) |> 
    pivot_longer(cols = starts_with("catch"), names_to = c("res", "run"), 
                 names_prefix = "catch_", names_pattern = "(.*)_(.*)", 
                 values_to = "vals")
  
  catch_obs <- list.files(file.path(base_dir, "monthly_weighted/1deg"),
                          pattern = "^dbpm_clim", full.names = T) |> 
    read_parquet() |> 
    filter(year >= 1961) 
  
  reg_code <- catch_obs |> 
    distinct(region) |> 
    pull()
  
  catch_obs <-  catch_obs|> 
    group_by(year) |> 
    summarise(obs_ccamlr = mean(catch_ccamlr*1e6, na.rm = T),
              obs_pauly = mean(catch_pauly*1e6, na.rm = T),
              obs_watson = mean(catch_tonnes_area_m2*1e6, na.rm = T)) |> 
    rowwise() |>
    mutate(min_catch_density = min(obs_ccamlr, obs_pauly, obs_watson,
                                   na.rm = T),
           max_catch_density = max(obs_ccamlr, obs_pauly, obs_watson, 
                                   na.rm = T)) |> 
    pivot_longer(cols = starts_with("obs"), names_to = "source", 
                 values_to = "obs") |> 
    mutate(min_catch_density = case_when(is.infinite(min_catch_density) ~ NA, 
                                         T ~ min_catch_density),
           max_catch_density = case_when(is.infinite(max_catch_density) ~ NA, 
                                         T ~ max_catch_density))
  
  catch_dbpm <- catch_dbpm |>
    left_join(catch_obs, by = "year") |> 
    #Calculating pseudo observations
    mutate(pseudo = case_when(vals < min_catch_density ~ min_catch_density,
                              vals > max_catch_density ~ max_catch_density,
                              .default = vals)) |> 
    mutate(region = reg_code, .after = run)
  
  return(catch_dbpm)
}


# Loading data ------------------------------------------------------------
base_dirs <- dir("/g/data/vf71/la6889/dbpm_inputs", pattern = "^[e|w]",
                 full.names = T)

catch_dbpm_obs <- data.frame()

for(fp in base_dirs){
  reg_name <- basename(fp) |> 
    str_replace_all("_", " ") |> 
    str_to_title()
  da <- get_catch_dbpm_obs(fp) |> 
    mutate(reg_name = reg_name, .after = region)
  catch_dbpm_obs <- catch_dbpm_obs |> 
    bind_rows(da)
}

# Save results
catch_dbpm_obs |> 
  write_parquet("outputs/simulated_observed_catches_all.parquet")

# Plotting modelled vs observed catches  ----------------------------------
si_catch <- catch_dbpm_obs |> 
  #Select runs with sea ice mask
  filter(run == "si") |> 
  ggplot()+
  geom_ribbon(aes(x = year, ymin = min_catch_density, 
                  ymax = max_catch_density), fill = "#d8d6dd")+
  geom_line(aes(year, vals, color = res), linewidth = 0.9)+
  geom_point(aes(year, vals, color = res, shape = res), size = 1.5)+
  geom_line(aes(year, obs, color = source), linetype = "dashed")+
  geom_point(aes(year, obs, color = source, shape = source), size = 0.2)+
  facet_grid(region~., scales = "free_y")+
  scale_color_manual(values = c("#d7301f", "#fc8d59",  "#045a8d", "#2b8cbe", 
                                "#74a9cf"))+
  ylab(~paste("Mean annual catch per unit area (t*", km^-2, ")"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(family = "sans", size = 14),
        axis.text = element_text(family = "sans", size = 12), 
        strip.text = element_text(family = "sans", size = 12),
        legend.position = "none", 
        plot.margin = margin(5, 15, 5, 5, "pt"))


# Creating figure to get legend for different DBPM resolutions
fig_leg <- catch_dbpm_obs |> 
  #Select runs with sea ice mask and calibration (non-spatial) for sample region
  filter((run == "si") & region == "FAO 58") |> 
  ggplot()+
  geom_line(aes(year, vals, color = res), linewidth = 0.9)+
  geom_point(aes(year, vals, color = res, shape = res), size = 1.5)+
  scale_color_manual("DBPM simulated catches", 
                     values = c("#d7301f", "#fc8d59"),
                     labels = c("0.25°", "1°"))+
  scale_shape_discrete("DBPM simulated catches", 
                       labels = c("0.25°", "1°"))+
  theme_bw()+
  theme(legend.title = element_text(family = "sans", size = 12, face = "bold", 
                                    hjust = 0.5), legend.title.position = "top",
        legend.byrow = T, legend.direction = "horizontal", 
        legend.position = "top",
        legend.text = element_text(family = "sans", size = 12))

leg_res <- get_plot_component(fig_leg, "guide-box-top", return_all = T)

# Creating figure to get legend for different observed catches sources
fig_leg <- catch_dbpm_obs |> 
  #Select runs with sea ice mask and calibration (non-spatial) for sample region
  filter((run == "si") & region == "FAO 58") |> 
  ggplot()+
  geom_line(aes(year, obs, color = source), linetype = "dashed")+
  geom_point(aes(year, obs, color = source, shape = source), size = 0.2)+
  scale_color_manual("Observed catches", 
                     values = c("#045a8d", "#2b8cbe", "#74a9cf"),
                     labels = c("CCAMLR", "Pauly & Zeller", "Watson & Tidd"))+
  scale_shape_discrete("Observed catches",
                       labels = c("CCAMLR", "Pauly & Zeller", "Watson & Tidd"))+
  theme_bw()+
  theme(legend.title = element_text(family = "sans", size = 12, face = "bold", 
                                    hjust = 0.5), legend.title.position = "top",
        legend.byrow = T, legend.direction = "horizontal", 
        legend.position = "top",
        legend.text = element_text(family = "sans", size = 12))

leg_obs <- get_plot_component(fig_leg, "guide-box-top", return_all = T)


fig <- plot_grid(plot_grid(leg_res, leg_obs, ncol = 2, rel_widths = c(0.55, 1)), 
                 si_catch, nrow = 2, rel_heights = c(0.15, 1))

ggsave("outputs/composite_fig_catches_obs.png", fig, bg = "white", width = 7.5, 
       height = 5)


## Comparison original and sea ice mask runs ----
si_catch_both_runs <- catch_dbpm_obs |> 
  filter(res != "calib") |> 
  ggplot()+
  geom_ribbon(aes(x = year, ymin = min_catch_density, 
                  ymax = max_catch_density), fill = "#d8d6dd")+
  geom_line(aes(year, vals, color = res, linetype = run))+
  geom_point(aes(year, vals, color = res, shape = res), size = 1.5)+
  geom_line(aes(year, obs, color = source), linetype = "dashed")+
  geom_point(aes(year, obs, color = source, shape = source), size = 0.2)+
  facet_grid(region~., scales = "free_y")+
  scale_color_manual("", values = c("#d7301f", "#fc8d59", "#045a8d", "#2b8cbe",
                                    "#74a9cf"))+
  scale_linetype_manual("Runs", values = c(1, 3), 
                        labels = c("original", "sea ice mask"))+
  ylab(~paste("Mean annual catch per unit area (t*", km^-2, ")"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), 
        legend.position = "none",
        axis.title.y = element_text(family = "sans", size = 14),
        axis.text = element_text(family = "sans", size = 12),
        strip.text = element_text(family = "sans", size = 12),
        plot.margin = margin(5, 15, 5, 5, "pt"))

# Creating plots to get legends
# Observed catches
fig_leg <- catch_dbpm_obs |>
  filter(region == "FAO 58" & res == "1") |> 
  ggplot()+
  geom_line(aes(year, vals, linetype = run))+
  scale_linetype_manual("DBPM run type", values = c(1, 3), 
                        labels = c("Original", "Sea ice adjusted"))+
  theme_bw()+
  theme(legend.title = element_text(family = "sans", size = 12, face = "bold", 
                                    hjust = 0.5), legend.position = "top",
        legend.byrow = T, legend.title.position = "top",
        legend.text = element_text(family = "sans", size = 12))

leg_runs <- get_plot_component(fig_leg, "guide-box-top", return_all = T)

fig2 <- plot_grid(plot_grid(leg_runs, leg_res, ncol = 2), si_catch_both_runs, 
                  leg_obs, nrow = 3, rel_heights = c(0.15, 1, 0.15))

ggsave("outputs/composite_fig_catches_obs_comp_runs.png", fig2, bg = "white")


# Calculating ORA metrics -------------------------------------------------
ora_metrics <- catch_dbpm_obs |> 
  group_by(res, run, region) |>
  summarise(pear_cor = cor(pseudo, vals),
            mae = sum(abs(vals-pseudo))/n(),
            ri = exp(sqrt(sum(log(pseudo/vals)^2)/n())),
            mef = ((sum((pseudo-mean(pseudo))^2)-sum((vals-mean(vals))^2))/
                     sum((pseudo-mean(pseudo))^2)))

fout <- paste0("outputs/ora_metric_catches_all-res_all_regions_all_runs_1961",
               "-2010.parquet")
ora_metrics |> 
  write_parquet(fout)
