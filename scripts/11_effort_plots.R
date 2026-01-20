## Comparing effort: CCAMLR vs Novaglio
source("scripts/useful_functions.R")
library(arrow)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(grid)
library(gridExtra)
library(zoo)
library(purrr)

fao_region <- 48
# Name of region
name_region <- "weddell"
reg_name_plot <- "Atlantic sector"

base_dir <- file.path("/g/data/vf71/la6889/dbpm_inputs", name_region)

#Creating path to ocean inputs
forcing_folder <- file.path("/g/data/vf71/la6889/dbpm_inputs", name_region, 
                            "monthly_weighted/1deg")

eff <- read_parquet(file.path(forcing_folder, 
                              paste0("dbpm_clim-fish-inputs_fao-", fao_region,
                                     "_1841-2010.parquet"))) |> 
  filter(year >= 1961)


eff_ccamlr <- read_parquet("data/effort_per_year_by_area.parquet") |> 
  filter(asd_area_code == fao_region & year <= 2010) |> 
  select(year, Effort_kW_Days_yr)


new_eff_params <- eff |> 
  select(!nom_active_area_m2_relative) |> 
  left_join(eff_ccamlr, by = "year") |> 
  rename(nom_active_area_m2_relative = Effort_kW_Days_yr) |> 
  mutate(nom_active_area_m2_relative = 
           if_else(is.na(nom_active_area_m2_relative), 
                   0, nom_active_area_m2_relative))

fishing_params <- read_parquet(
  list.files(file.path(base_dir, "fishing_params",
                       paste0("best_fish_vals_fao-", fao_region, 
                              "_1deg")), 
             pattern = "1000.parquet", full.names = T)) |> 
  filter(cor >= 0) |> 
  arrange(desc(cor), rmse) |> 
  slice(1)

non_spatial_run <- run_model(fishing_params, new_eff_params)

new_run <- non_spatial_run |> 
  select(c(time, year, total_catch)) |> 
  group_by(year) |> 
  summarise(new_effort = mean(total_catch, na.rm = T))

est_catches <- catch_dbpm_obs |> 
  filter(res == "catch_nonspat") |> 
  select(year, vals) |> 
  rename(old_effort = vals) |> 
  left_join(new_run) |> 
  pivot_longer(!year, names_to = "run", values_to = "estimates")

obs_catches <- catch_dbpm_obs |> 
  distinct(year, source, obs) 
  
  
ggplot()+
  geom_line(data = est_catches, aes(year, estimates, color = run))+
  geom_line(data = obs_catches, aes(year, obs, color = source), 
            linetype = "dashed")
  
  filter(res == "catch_nonspat") |> 
  select(year, vals, source, obs) |> 
  left_join(new_run, by = "year") |> 
  pivot_longer(c(vals, new_effort_run), names_to = "run", 
               values_to = "catch") |> 
  ggplot()+
  geom_line(aes(year, catch, color = run))+
  geom_point(aes(year, obs, color = source))


effort_both <- eff |> 
  select(year, total_nom_active) |> 
  full_join(eff_ccamlr, by = "year") |> 
  arrange(year) |> 
  rename(ccmalr = Effort_kW_Days_yr, watson = total_nom_active) |> 
  pivot_longer(!year, names_to = "source", values_to = "effort")


fao48 <- effort_both |>
  ggplot(aes(year, effort, colour = source))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = c("#1f78b4", "#33a02c"),
                     labels = c("CCAMLR (regional)", "FishMIP (global)"))+
  guides(color = guide_legend(title = "Fishing effort sources"))+
  labs(y = ~paste("Total annual fishing effort (kW*days*", km^-2, ")"),
       subtitle = paste0(str_to_title(reg_name_plot), " (FAO ", fao_region,
                         ")"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(family = "sans", size = 14),
        axis.text = element_text(family = "sans", size = 12),
        legend.title = element_text(face = "bold"),
        legend.position = "right",
        legend.text = element_text(family = "sans", size = 12), 
        plot.margin = margin(5, 15, 5, 5, "pt"),
        plot.subtitle = element_text(family = "sans", size = 12, hjust = 0.5))


leg <- get_legend(fao58)

fout <- file.path(base_dir, "gridded_dbpm_outputs", 
                  paste0("effort_ts_watson-ccamlr_fao-", fao_region, 
                         "_1961-2010.png"))

ggsave(fout, device = "png")

fao48_noleg <- fao48+theme(legend.position = "none",
                           axis.title.y = element_blank())
fao58_noleg <- fao58+theme(legend.position = "none",
                           axis.title.y = element_blank())
fao88_noleg <- fao88+theme(legend.position = "none", 
                           axis.title.y = element_blank())+
  scale_y_continuous(labels = scales::scientific)

ytitle <- textGrob(~paste("Total annual fishing effort (kW*days*", km^-2, ")"),
                   rot = 90)

fig <- grid.arrange(arrangeGrob(plot_grid(plot_grid(fao48_noleg, fao58_noleg, 
                                                    ncol = 2), 
                                          plot_grid(fao88_noleg, leg), 
                                          nrow = 2), left = ytitle))

ggsave("outputs/composite_fig_effort_ccamlr-novaglio.png", fig)
