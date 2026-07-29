source("/Users/lakshmi_nccs/Desktop/NCCS/Projects/For_Publication/HNSCC_CosMx6k/helpers/design_elements.R")
source("/Users/lakshmi_nccs/Desktop/NCCS/Projects/For_Publication/HNSCC_CosMx6k/helpers/plotting_functions.R")

plot_results_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/Paper_Figures/Figure_3_µm"
library('anndata')
library('patchwork')
library(rstatix)

ad = reticulate::import("anndata")

pat_levels1 = c("2-011","2-020","2-017","2-028","2-013","2-026",
               "2-022","2-034","2-039","2-033","2-018","2-032",
               "2-036","2-040","20-010","20-002","2-041","2-012","2-029")
pat_levels2 =  c(
  "2-011", "2-012", "2-013", "2-017", "2-018", "2-020", "2-022", 
  "2-026", "2-028", "2-029", "2-032", "2-033", "2-034", "2-036", 
  "2-039", "2-040", "2-041", "20-002", "20-010"
)
old_tt_levels = c("Pre-NIVO Primary", "Post-NIVO Primary", "Post-NIVO Adv Front", "Post-NIVO Normal")
tt_levels = c("Pre-NIVO primary", "Post-NIVO primary", "Post-NIVO adv front", "Post-NIVO normal")
recode_tt_map <- setNames(tt_levels,old_tt_levels)

old_ct_levels = c("B Cells", "Dentritic Cells", "Endothelial Cells", "Epithelial Cells",
                  "Fibroblasts", "Macrophages", "Mast Cells", "Myogenic Cells",
                  "Neuronal Cells", "Neutrophils", "Plasma Cells", "T Cells", "Tumor")
ct_levels = c("B cells", "Dendritic cells", "Endothelial cells", "Epithelial cells",
              "Fibroblasts", "Macrophages", "Mast cells", "Myogenic cells",
              "Neuronal cells", "Neutrophils", "Plasma cells", "T cells", "Tumor")
recode_ct_map <- setNames(ct_levels,old_ct_levels)

old_resp_levels = c("Responder", "Non_Responder")
resp_levels = c("Response", "Non-response")
recode_resp_map <- setNames(resp_levels,old_resp_levels)

old_bin_levels = c("Q0: TN", "Q1: 0-250", "Q2: 250-500", "Q3: 500-800", "Q4: 800+")
bin_levels = c("B0: TN", "B1: 0–250", "B2: 250–500", "B3: 500–800", "B4: 800+")
bin_levels = c("B0: TN", "B1: 0–30", "B2: 30–60", "B3: 60–96", "B4: 96+")
recode_bin_map <- setNames(bin_levels,old_bin_levels)

old_niche = c("Epithelial_rich_niche","Fibroblasts_rich_niche","Neuronal_rich_niche","Plasma_rich_niche","T_rich_niche","Tumor_rich_niche")
new_niche = c("Epithelial-rich niche","Fibroblasts-rich niche","Neuronal-rich niche","Plasma-rich niche","T-rich niche","Tumor-rich niche")
recode_niche_map <- setNames(new_niche,old_niche)

old_con_niche <- c("B Cells_rich_niche", "Dentritic Cells_rich_niche","Endothelial Cells_rich_niche","Epithelial Cells_rich_niche",
                   "Fibroblasts_rich_niche","Macrophages_rich_niche","Mast Cells_rich_niche","Myogenic Cells_rich_niche","Neuronal Cells_rich_niche",
                   "Neutrophils_rich_niche","Plasma Cells_rich_niche","T Cells_rich_niche","Tumor_rich_niche", "No-contact")
new_con_niche <- c("B-rich niche","Dendritic-rich niche","Endothelial-rich niche","Epithelial-rich niche",
                   "Fibroblasts-rich niche","Macrophages-rich niche","Mast-rich niche","Myogenic-rich niche","Neuronal-rich niche",
                   "Neutrophils-rich niche","Plasma-rich niche","T-rich niche","Tumor-rich niche","No-contact niche")
recode_con_niche_map <- setNames(new_con_niche,old_con_niche)


###### CD8 Distance to Tum Histogram ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_tcells_dis_tn_tum"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary", "Post-NIVO adv front"), ]
adata = adata[adata$obs$Distance_to_Tumor <= 800,]

adata$obs$Distance_to_Tumor = adata$obs$Distance_to_Tumor * 0.120

p1 = plot_hist_lines_points(adata$obs, "Distance_to_Tumor", "TumorType", 50, tt_colors, 
                         c(" ", "Distance to tumor (µm)", "Number of CD8+ CD4- T cells"))+
  theme(plot.title = element_blank()) + labs(color = "Tumor type") +
  theme(legend.position = c(0.7, 0.7), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "CD8_hist_spat", 6,4)


###### CD4 Distance to Tum Histogram ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd4_tcells_dis_tn_tum"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary", "Post-NIVO adv front"), ]
adata = adata[adata$obs$Distance_to_Tumor <= 2000,]
adata$obs$Distance_to_Tumor = adata$obs$Distance_to_Tumor * 0.120
p1 = plot_hist_lines_points(adata$obs, "Distance_to_Tumor", "TumorType", 50, tt_colors, 
                              c(" ", "Distance to tumor  (µm)", "Number of CD8- CD4+ T cells"))+
  theme(plot.title = element_blank()) + labs(color = "Tumor type") +
  theme(legend.position = c(0.7, 0.7), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "CD4_hist_spat", 5,6)

###### Proportion per bin plot ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_tcells_dis_tn_tum"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

resp_df <- readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/response.xlsx", sheet = "Sheet1")
adata$obs <- adata$obs %>%
  left_join(resp_df, by = "Patient")
adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"), ]
#adata = adata[adata$obs$Distance_to_Tumor <= 1000,]
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)
adata$obs$spatial_bin = recode(as.character(adata$obs$spatial_bin), !!!recode_bin_map)
adata$obs$spatial_bin = factor(adata$obs$spatial_bin, levels = bin_levels)

for (i in unique(adata$obs$T_subtypes_l1)){
  print(i)
  p_final =  plot_prop_per_bin_overall(
    adata = adata,
    y_obs = "T_subtypes_l1",
    y_value = i,
    x_obs = "spatial_bin",
    fill_obs = "TumorType",
    split_obs = "Response",
    t_subtype_col = "T_subtypes_l1"
  )

  p_final = p_final + theme(panel.background = element_rect(fill = "white", color = NA),
                            plot.background = element_rect(fill = "white", color = NA))
    
  plot_save_tiff(p_final, paste0("PP_propplot_", i), 5,4)
}

###### Enrichment Bars: Responder ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_tcells_dis_tn_tum"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

resp_df <- readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/response.xlsx", sheet = "Sheet1")
adata$obs <- adata$obs %>%
  left_join(resp_df, by = "Patient")
adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"), ]
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

df = adata$obs

df <- df %>%
  mutate(spatial_bin = cut(
    Distance_to_Tumor,
    breaks = c(0, 500, 4354.142),
    labels = c("Near", "Far"),
    right = FALSE,   # intervals like [0,50), [50,100)
    include.lowest = TRUE
  ))

df = df[df$Response == "Response",]
agg_fov <- df %>%
  group_by(tma_sid, spatial_bin, Response, TumorType, T_subtypes_l1) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(tma_sid, spatial_bin, Response, TumorType) %>%
  mutate(
    Total = sum(Count),
    Proportion = Count / Total
  ) %>%
  ungroup()
plot_df <- agg_fov %>%
  group_by(TumorType, spatial_bin, T_subtypes_l1) %>%
  summarise(
    mean_prop = median(Proportion),
    se_prop   = sd(Proportion) / sqrt(n()),
    .groups = "drop"
  )
plot_df <- plot_df %>%
  mutate(
    plot_value = ifelse(spatial_bin == "Near", mean_prop, -mean_prop),
    ymin = plot_value - se_prop,
    ymax = plot_value + se_prop
  )
plot_df$spatial_bin = factor(plot_df$spatial_bin, levels = c("Near", "Far"))
plot_df$group4 <- interaction(plot_df$spatial_bin,plot_df$TumorType,sep = "::")
plot_df$group4 <- factor(plot_df$group4,levels = c("Near::Pre-NIVO primary","Near::Post-NIVO primary",
                                                   "Far::Pre-NIVO primary","Far::Post-NIVO primary"))

stat_df <- agg_fov %>%
  group_by(T_subtypes_l1, spatial_bin) %>%
  filter(n_distinct(TumorType) == 2) %>% 
  rstatix::wilcox_test(Proportion ~ TumorType) %>%   # TumorType = Pre/Post
  add_significance("p")
stat_df <- stat_df %>%
  left_join(
    plot_df %>%
      group_by(T_subtypes_l1, spatial_bin) %>%
      summarise(
        ymax_bar = max(ymax),
        ymin_bar = min(ymin),
        .groups = "drop"
      ),
    by = c("T_subtypes_l1", "spatial_bin")
  ) %>%
  mutate(
    x = T_subtypes_l1,
    y.position = ifelse(spatial_bin == "Near", ymax_bar + 0.1, ymin_bar - 0.1),
    group = NA, 
    label = ifelse(p.signif == "ns", "", p)
  )

p1 = ggplot(plot_df,aes(x = T_subtypes_l1,y = plot_value,fill = group4,group = TumorType)) +
  geom_col(position = position_dodge(width = 0.7),width = 0.6,colour = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                position = position_dodge(width = 0.7),width = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  theme_Publication(base_size = 12) +
  scale_fill_manual(
    values = c(
      "Near::Pre-NIVO primary" = "#4363d8",
      "Near::Post-NIVO primary" = "#3cb44b",
      "Far::Pre-NIVO primary" = "#142B7A",
      "Far::Post-NIVO primary" = "#145A1F")) +
  scale_y_continuous(labels = abs) +
  labs(y = "Median proportion of subtype per FOV \n Response",x = "CD4- CD8+ subtype",fill = "Status") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1))+
  theme( axis.text.x = element_text(angle = 15, hjust = 1), legend.title = element_text(hjust = 0.5))+
  #theme(legend.position = "bottom")+
  theme( axis.text.x = element_text(size = 10))+
  ggpubr::stat_pvalue_manual(
    stat_df,
    label = "label",
    x = "T_subtypes_l1",
    y.position = "y.position",
    group = "group",   # NA → line spans both bars
    tip.length = 0,
    size = 3
  )

plot_save_tiff(p1, "CD8_enrich_Responders", 6,9)


###### Enrichment Bars: Non_Responder ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_tcells_dis_tn_tum"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

resp_df <- readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/response.xlsx", sheet = "Sheet1")
adata$obs <- adata$obs %>%
  left_join(resp_df, by = "Patient")
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"), ]
df = adata$obs

df <- df %>%
  mutate(spatial_bin = cut(
    Distance_to_Tumor,
    breaks = c(0, 500, 4354.142),
    labels = c("Near", "Far"),
    right = FALSE,   # intervals like [0,50), [50,100)
    include.lowest = TRUE
  ))

df = df[df$Response == "Non-response",]
agg_fov <- df %>%
  group_by(tma_sid, spatial_bin, Response, TumorType, T_subtypes_l1) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(tma_sid, spatial_bin, Response, TumorType) %>%
  mutate(
    Total = sum(Count),
    Proportion = Count / Total
  ) %>%
  ungroup()

plot_df <- agg_fov %>%
  group_by(TumorType, spatial_bin, T_subtypes_l1) %>%
  summarise(
    mean_prop = median(Proportion),
    se_prop   = sd(Proportion) / sqrt(n()),
    .groups = "drop"
  )
plot_df <- plot_df %>%
  mutate(
    plot_value = ifelse(spatial_bin == "Near", mean_prop, -mean_prop),
    ymin = plot_value - se_prop,
    ymax = plot_value + se_prop
  )
plot_df$spatial_bin = factor(plot_df$spatial_bin, levels = c("Near", "Far"))
plot_df$group4 <- interaction(plot_df$spatial_bin,plot_df$TumorType,sep = "::")
plot_df$group4 <- factor(plot_df$group4,levels = c("Near::Pre-NIVO primary","Near::Post-NIVO primary",
                                                   "Far::Pre-NIVO primary","Far::Post-NIVO primary"))

stat_df <- agg_fov %>%
  group_by(T_subtypes_l1, spatial_bin) %>%
  filter(n_distinct(TumorType) == 2) %>% 
  rstatix::wilcox_test(Proportion ~ TumorType) %>%   # TumorType = Pre/Post
  add_significance("p")
stat_df <- stat_df %>%
  left_join(
    plot_df %>%
      group_by(T_subtypes_l1, spatial_bin) %>%
      summarise(
        ymax_bar = max(ymax),
        ymin_bar = min(ymin),
        .groups = "drop"
      ),
    by = c("T_subtypes_l1", "spatial_bin")
  ) %>%
  mutate(
    x = T_subtypes_l1,
    y.position = ifelse(spatial_bin == "Near", ymax_bar + 0.1, ymin_bar - 0.1),
    group = NA, 
    label = ifelse(p.signif == "ns", "", p)
  )

p1 = ggplot(plot_df,aes(x = T_subtypes_l1,y = plot_value,fill = group4,group = TumorType)) +
  geom_col(position = position_dodge(width = 0.7),width = 0.6,colour = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                position = position_dodge(width = 0.7),width = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  theme_Publication(base_size = 12) +
  scale_fill_manual(
    values = c(
      "Near::Pre-NIVO primary" = "#4363d8",
      "Near::Post-NIVO primary" = "#3cb44b",
      "Far::Pre-NIVO primary" = "#142B7A",
      "Far::Post-NIVO primary" = "#145A1F")) +
  scale_y_continuous(labels = abs) +
  labs(y = "Median proportion of subtype per FOV \n Non-response",x = "CD4- CD8+ subtype",fill = "Status") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1))+
  theme( axis.text.x = element_text(angle = 15, hjust = 1))+
  theme( axis.text.x = element_text(size = 10), legend.title = element_text(hjust = 0.5))+
  ggpubr::stat_pvalue_manual(
    stat_df,
    label = "label",
    x = "T_subtypes_l1",
    y.position = "y.position",
    group = "group",   # NA → line spans both bars
    tip.length = 0,
    size = 3
  )

plot_save_tiff(p1, "CD8_enrich_Non_Responders", 6,9)


####### Landscape of differences: Responder ######
library(ggridges)
library(ggstream)
library(dplyr)
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_tcells_dis_tn_tum"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

resp_df <- readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/response.xlsx", sheet = "Sheet1")
adata$obs <- adata$obs %>%
  left_join(resp_df, by = "Patient")
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"), ]
ridge_df <- adata$obs %>%
  #filter(TumorType == "Post-NIVO Primary") %>%
  dplyr::select(Distance_to_Tumor, T_subtypes_l1, Response, TumorType)
ridge_df <- ridge_df %>%
  filter(!is.na(T_subtypes_l1) & !is.na(Response) & !is.na(Distance_to_Tumor))


resp_df <- ridge_df[ridge_df$Response == "Response",]
resp_df = resp_df[resp_df$Distance_to_Tumor <= 2000,]
bin_width <- 75
bins <- seq(0, max(resp_df$Distance_to_Tumor, na.rm = TRUE) + bin_width, by = bin_width)
df_binned <- resp_df %>%
  mutate(distance_bin = cut(Distance_to_Tumor, breaks = bins, include.lowest = TRUE, right = TRUE)) %>%
  group_by(distance_bin, T_subtypes_l1, TumorType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(distance_bin) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  # compute numeric midpoints robustly
  mutate(distance_bin_char = as.character(distance_bin)) %>%
  rowwise() %>%
  mutate(
    lower = as.numeric(strsplit(gsub("\\(|\\[|\\]|\\)", "", distance_bin_char), ",")[[1]][1]),
    upper = as.numeric(strsplit(gsub("\\(|\\[|\\]|\\)", "", distance_bin_char), ",")[[1]][2]),
    Distance = (lower + upper)/2
  ) %>%
  ungroup()

df_diff <- df_binned %>%
  dplyr::select(Distance, T_subtypes_l1, TumorType, proportion) %>%
  pivot_wider(names_from = TumorType, values_from = proportion, values_fill = 0) %>%
  mutate(diff = `Post-NIVO primary` - `Pre-NIVO primary`)

df_stream_diff <- df_diff %>%
  dplyr::select(Distance, T_subtypes_l1, diff) %>%
  dplyr::rename(scaled_expr = diff, category = T_subtypes_l1)

df_stream_scaled <- df_stream_diff %>%
  group_by(Distance) %>%   # scale within each distance bin
  mutate(
    scaled_expr = (scaled_expr - min(scaled_expr)) / (max(scaled_expr) - min(scaled_expr))
  ) %>%
  ungroup()
df_stream_scaled$Distance = df_stream_scaled$Distance * 0.120
p1 = ggplot(df_stream_scaled, aes(x = Distance, y = scaled_expr, fill = category)) +
  geom_stream(type = "proportional", bw = 0.5) +
  theme_Publication(base_size = 14) +
  labs(
    x = "Distance to tumor (µm)",
    y = "CD8 Spatial enrichment score \n Response", #Scaled (Post − Pre) \n (Non_Responders)
    fill = "CD4- CD8+ subtype",
    title = NULL
  )+
  ggplot2::scale_fill_manual(values = cd8_l1)+
  scale_x_continuous(breaks = seq(0, 240, by = 30),expand = c(0, 0)) +   # 👈 removes padding
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.line = element_blank(),           # ❌ remove axis lines
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2, "pt"),     # smaller tick spacing
    plot.margin = margin(10, 10, 10, 10)       # tighter margins
  ) +
  theme(legend.title = element_text(hjust = 0.5))

plot_save_tiff(p1, "CD8_landscape_Resp", 5,11)


####### Landscape of differences: Non_Responder ######
library(ggridges)
library(ggstream)
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_tcells_dis_tn_tum"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

resp_df <- readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/response.xlsx", sheet = "Sheet1")
adata$obs <- adata$obs %>%
  left_join(resp_df, by = "Patient")
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"), ]
ridge_df <- adata$obs %>%
  #filter(TumorType == "Post-NIVO Primary") %>%
  dplyr::select(Distance_to_Tumor, T_subtypes_l1, Response, TumorType)
ridge_df <- ridge_df %>%
  filter(!is.na(T_subtypes_l1) & !is.na(Response) & !is.na(Distance_to_Tumor))


resp_df <- ridge_df[ridge_df$Response == "Non-response",]
resp_df = resp_df[resp_df$Distance_to_Tumor <= 2000,]
bin_width <- 75
bins <- seq(0, max(resp_df$Distance_to_Tumor, na.rm = TRUE) + bin_width, by = bin_width)
df_binned <- resp_df %>%
  mutate(distance_bin = cut(Distance_to_Tumor, breaks = bins, include.lowest = TRUE, right = TRUE)) %>%
  group_by(distance_bin, T_subtypes_l1, TumorType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(distance_bin) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  # compute numeric midpoints robustly
  mutate(distance_bin_char = as.character(distance_bin)) %>%
  rowwise() %>%
  mutate(
    lower = as.numeric(strsplit(gsub("\\(|\\[|\\]|\\)", "", distance_bin_char), ",")[[1]][1]),
    upper = as.numeric(strsplit(gsub("\\(|\\[|\\]|\\)", "", distance_bin_char), ",")[[1]][2]),
    Distance = (lower + upper)/2
  ) %>%
  ungroup()

df_diff <- df_binned %>%
  dplyr::select(Distance, T_subtypes_l1, TumorType, proportion) %>%
  pivot_wider(names_from = TumorType, values_from = proportion, values_fill = 0) %>%
  mutate(diff = `Post-NIVO primary` - `Pre-NIVO primary`)

df_stream_diff <- df_diff %>%
  dplyr::select(Distance, T_subtypes_l1, diff) %>%
  dplyr::rename(scaled_expr = diff, category = T_subtypes_l1)

df_stream_scaled <- df_stream_diff %>%
  group_by(Distance) %>%   # scale within each distance bin
  mutate(
    scaled_expr = (scaled_expr - min(scaled_expr)) / (max(scaled_expr) - min(scaled_expr))
  ) %>%
  ungroup()

df_stream_scaled$Distance = df_stream_scaled$Distance * 0.120
p1 = ggplot(df_stream_scaled, aes(x = Distance, y = scaled_expr, fill = category)) +
  geom_stream(type = "proportional", bw = 0.5) +
  theme_Publication(base_size = 14) +
  labs(
    x = "Distance to tumor (µm)",
    y = "CD8 spatial enrichment score \n Non-response", #Scaled (Post − Pre) \n (Non_Responders)
    fill = "CD4- CD8+ subtype",
    title = NULL
  )+
  ggplot2::scale_fill_manual(values = cd8_l1)+
  scale_x_continuous(breaks = seq(0, 240, by = 30),expand = c(0, 0)) +   # 👈 removes padding
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.line = element_blank(),           # ❌ remove axis lines
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2, "pt"),     # smaller tick spacing
    plot.margin = margin(10, 10, 10, 10)       # tighter margins
  ) +
  theme(legend.title = element_text(hjust = 0.5))

plot_save_tiff(p1, "CD8_landscape_Non_Resp", 5,11)

###### Heatmap of Distances of selected genes ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_tcells_dis_tn_tum"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

resp_df <- readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/response.xlsx", sheet = "Sheet1")
adata$obs <- adata$obs %>%
  left_join(resp_df, by = "Patient")
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)
adata$obs$Distance_to_Tumor = adata$obs$Distance_to_Tumor*0.120
adata_post = adata[adata$obs$TumorType %in% c("Post-NIVO primary"), ]
adata_pre = adata[adata$obs$TumorType %in% c("Pre-NIVO primary"), ]

adata_post_R = adata_post[adata_post$obs$Response %in% c("Response"), ]
adata_post_NR = adata_post[adata_post$obs$Response %in% c("Non-response"), ]
adata_pre_R = adata_pre[adata_pre$obs$Response %in% c("Response"), ]
adata_pre_NR = adata_pre[adata_pre$obs$Response %in% c("Non-response"), ]

gns = c("KRT5", "LAG3", "LGALS3", "HLA-DRA", "HLA-DRB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2",
        "HLA-DQB1", "HLA-DQB2")
gns = gns[gns %in% adata$var_names]
gns = c('KRT5', 'IGKC', 'COL1A1', 'PECAM1', 'CD68',
        "LGALS3","LGALS1","CXCL8","GNLY","PRF1",
        "IL7R",  "CTLA4", "LAG3", "CXCL13", "TOX", "TIGIT", 
         "GZMH", "CCL5", "CXCR4","KLF2","BMI1","LEF1" )
p1 = plot_gene_distance_heatmap(adata_post_R, 
                                gns, 
                                cbreaks = c(c(0, 50, 100, 150),seq(175, 620, by = 20), 
                                            seq(700, 1000, by = 100),
                                            c(2000, 3000))*0.120, 
                                max_distance = 3000*0.120, tit = "Post-R",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Heatmap_Dis_PostR.tiff"), 
       plot =grid::grid.draw(p1),
       width = 4, height = 4, bg = "white")


p2 = plot_gene_distance_heatmap(adata_post_NR, 
                                gns, 
                                cbreaks = c(c(0, 50, 100, 150),seq(175, 620, by = 20), 
                                            seq(700, 1000, by = 100),
                                            c(2000, 3000))*0.120, 
                                max_distance = 3000*0.120,tit = "Post-NR",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Heatmap_Dis_PostNR.tiff"), 
       plot =grid::grid.draw(p2),
       width = 4, height = 4, bg = "white")

p3 = plot_gene_distance_heatmap(adata_pre_R, 
                                gns, 
                                cbreaks = c(c(0, 50, 100, 150),seq(175, 620, by = 20), 
                                            seq(700, 1000, by = 100),
                                            c(2000, 3000))*0.120, 
                                max_distance = 3000*0.120,tit = "Pre-R",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Heatmap_Dis_PreR.tiff"), 
       plot =grid::grid.draw(p3),
       width = 4, height = 4, bg = "white")

p4 = plot_gene_distance_heatmap(adata_pre_NR, 
                                gns, 
                                cbreaks = c(c(0, 50, 100, 150),seq(175, 620, by = 20), 
                                            seq(700, 1000, by = 100),
                                            c(2000, 3000))*0.120, 
                                max_distance = 3000*0.120,tit = "Pre-NR",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Heatmap_Dis_PreNR.tiff"), 
       plot =grid::grid.draw(p4),
       width = 4, height = 4, bg = "white")



###### Niche Figures ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_tcells_spfp_500_niche"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$spatial_niche = recode(as.character(adata$obs$spatial_niche), !!!recode_niche_map)
adata$obs$spatial_niche = factor(adata$obs$spatial_niche, levels = new_niche)

p1 = plot_umap(dat = adata, 1, 2, FALSE, 1, niche_colors_spat, adata$obs$spatial_niche, rep(1, length(adata$obs$spatial_niche)),
               c(paste0("PN-Tum"), "Spatial niche", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(ncol = 1, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "right", plot.title = element_blank(), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "SN_UMAP", 5,7)


data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_spat_metrics"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

resp_df <- readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/response.xlsx", sheet = "Sheet1")
adata$obs <- adata$obs %>%
  left_join(resp_df, by = "Patient")
adata$obs$spatial_niche = recode(as.character(adata$obs$spatial_niche), !!!recode_niche_map)
adata$obs$spatial_niche = factor(adata$obs$spatial_niche, levels = new_niche)
adata$obs$contact_niche = recode(as.character(adata$obs$contact_niche), !!!recode_con_niche_map)
adata$obs$contact_niche = factor(adata$obs$contact_niche, levels = new_con_niche)
adata$obs$spatial_bin = recode(as.character(adata$obs$spatial_bin), !!!recode_bin_map)
adata$obs$spatial_bin = factor(adata$obs$spatial_bin, levels = bin_levels)

p1 = plot_umap(dat = adata, 1, 2, FALSE, 1.2, niche_colors_con, adata$obs$contact_niche, rep(1, length(adata$obs$contact_niche)),
               c(paste0("PN-Tum"), "Contact niche", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(ncol = 1, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "right", plot.title = element_blank(), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "CN_UMAP", 6,9)


p1 = plot_umap(dat = adata, 1, 2, FALSE, 1, niche_colors_spat, adata$obs$spatial_niche, rep(1, length(adata$obs$spatial_niche)),
               c(paste0("PN-Tum"), "Spatial niche", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(ncol = 1, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "right", plot.title = element_blank(), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "SN_all_UMAP", 5,7)


p1 = plot_umap(dat = adata, 1, 2, FALSE, 1, distance_colors_orange, adata$obs$spatial_bin, rep(1, length(adata$obs$spatial_niche)),
               c(paste0("PN-Tum"), "Spatial bin", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(ncol = 1, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "right", plot.title = element_blank(), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "SpatBins_all_UMAP", 5,7)


###### Fibroblasts Niche Post-Nivo ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_spat_metrics"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

resp_df <- readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/response.xlsx", sheet = "Sheet1")
adata$obs <- adata$obs %>%
  left_join(resp_df, by = "Patient")
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

adata$obs$spatial_bin = factor(adata$obs$spatial_bin, levels = names(distance_colors_orange))
adata$obs$spatial_niche = recode(as.character(adata$obs$spatial_niche), !!!recode_niche_map)
adata$obs$spatial_niche = factor(adata$obs$spatial_niche, levels = new_niche)
adata$obs$contact_niche = recode(as.character(adata$obs$contact_niche), !!!recode_con_niche_map)
adata$obs$contact_niche = factor(adata$obs$contact_niche, levels = new_con_niche)
adata$obs$spatial_bin = recode(as.character(adata$obs$spatial_bin), !!!recode_bin_map)
adata$obs$spatial_bin = factor(adata$obs$spatial_bin, levels = bin_levels)

adata = adata[adata$obs$TumorType %in% c("Post-NIVO primary"), ]$copy()
adata$obs$Response <- recode(
  adata$obs$Response,
  "Response"  = "Post-R",
  "Non-response" = "Post-NR"
)
pl = list()
for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_stacked_bar_percentage(ad_here$obs, "Response", "spatial_niche", c(niche_colors_spat), TRUE, FALSE, c(i, " "))+
    ggplot2::theme(legend.position = "bottom") + 
    theme(plot.title = element_blank(),
          axis.text.x = element_text(margin = margin(t = -8)))+
    labs(fill = "Spatial niche")+
    guides(fill = guide_legend(nrow = 1, keywidth = 1,keyheight = 1))+
    theme(plot.margin = margin(5.5, 4, 15, 10), legend.title = element_text(hjust = 0.5)) 
  if (i != "B0: TN") {p <- p +theme(axis.title.y = element_blank())}
  pl[[i]] = p
}
p1 = ggpubr::ggarrange(plotlist = pl, common.legend = TRUE, nrow = 1, legend = "bottom") +
  theme(panel.background = element_rect(fill = "white", color = NA),
                plot.background = element_rect(fill = "white", color = NA))
plot_save_tiff(p1, "SN_Prop_Plot_Post", 5,10)
plot_save_tiff(p1, "SN_Prop_Plot_Post_legend", 5,13)

pl = list()
for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  ad_here$obs$contact_niche = factor(ad_here$obs$contact_niche, levels = new_con_niche)
  p = plot_stacked_bar_percentage(ad_here$obs, "Response", "contact_niche", c(niche_colors_con), TRUE, FALSE, c(i, " "))+
    ggplot2::theme(legend.position = "none") + 
    theme(plot.title = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 1, margin = margin(t = -8)))+
    labs(fill = "Contact niche")+
    #guides(fill = guide_legend(ncol = 1, keywidth = 1,keyheight = 1))+
    theme(plot.margin = margin(5.5, 4, 15, 10)) 
  if (i != "B0: TN") {p <- p +theme(axis.title.y = element_blank())}
  pl[[i]] = p
}
p1 = ggpubr::ggarrange(plotlist = pl, common.legend = FALSE, ncol = 5) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
plot_save_tiff(p1, "CN_Prop_Plot_Post", 5,8)

p = plot_stacked_bar_percentage(adata$obs, "Response", "contact_niche", c(niche_colors_con), TRUE, FALSE, c(i, " "))+
  ggplot2::theme(legend.position = "bottom") + 
  theme(plot.title = element_blank(),
        axis.text.x = element_text(margin = margin(t = -8)))+
  labs(fill = "Contact niche")+
  guides(fill = guide_legend(nrow = 2, keywidth = 1,keyheight = 1))+
  theme(plot.margin = margin(5.5, 4, 15, 10)) +
  theme(legend.title = element_text(hjust = 0.5))
plot_save_tiff(p, "CN_Prop_Plot_Post_legend", 4,16.5)


for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_response_boxplot(ad_here$obs, "tma_sid", "Fibroblasts-rich niche","spatial_niche")+
    theme(legend.position = "none")+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  plot_save_tiff(p, paste0("SN_Post_Fibro_", i), 5,3)
}

for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_response_boxplot(ad_here$obs, "tma_sid", "Plasma-rich niche","spatial_niche")+
    theme(legend.position = "none")+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  plot_save_tiff(p, paste0("SN_Post_Plasma_", i), 5,3)
}


for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_response_boxplot(ad_here$obs, "tma_sid", "Fibroblasts-rich niche","contact_niche")
  plot_save_tiff(p, paste0("CN_Post_Fibro_", i), 6,6)
}


###### Fibroblasts Niche Pre-Nivo ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P3/")
nm = "cd8_spat_metrics"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

resp_df <- readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/response.xlsx", sheet = "Sheet1")
adata$obs <- adata$obs %>%
  left_join(resp_df, by = "Patient")
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

adata$obs$spatial_bin = factor(adata$obs$spatial_bin, levels = names(distance_colors_orange))
adata$obs$spatial_niche = recode(as.character(adata$obs$spatial_niche), !!!recode_niche_map)
adata$obs$spatial_niche = factor(adata$obs$spatial_niche, levels = new_niche)
adata$obs$contact_niche = recode(as.character(adata$obs$contact_niche), !!!recode_con_niche_map)
adata$obs$contact_niche = factor(adata$obs$contact_niche, levels = new_con_niche)
adata$obs$spatial_bin = recode(as.character(adata$obs$spatial_bin), !!!recode_bin_map)
adata$obs$spatial_bin = factor(adata$obs$spatial_bin, levels = bin_levels)

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary"), ]$copy()
adata$obs$Response <- recode(
  adata$obs$Response,
  "Response"  = "Pre-R",
  "Non-response" = "Pre-NR"
)
pl = list()
for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_stacked_bar_percentage(ad_here$obs, "Response", "spatial_niche", c(niche_colors_spat), TRUE, FALSE, c(i, " "))+
    ggplot2::theme(legend.position = "bottom") + 
    theme(plot.title = element_blank(),
          axis.text.x = element_text(margin = margin(t = -8)))+
    labs(fill = "Spatial niche")+
    guides(fill = guide_legend(nrow = 2, keywidth = 1,keyheight = 1))+
    theme(plot.margin = margin(5.5, 4, 15, 10)) 
  if (i != "B0: TN") {p <- p +theme(axis.title.y = element_blank())}
  pl[[i]] = p
}
p1 = ggpubr::ggarrange(plotlist = pl, common.legend = TRUE, nrow = 1, legend = "bottom") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
plot_save_tiff(p1, "SN_Prop_Plot_Pre", 6,9)

pl = list()
for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_stacked_bar_percentage(ad_here$obs, "Response", "contact_niche", c(niche_colors_con), TRUE, FALSE, c(i, " "))+
    ggplot2::theme(legend.position = "none") + 
    theme(plot.title = element_blank(),
          axis.text.x = element_text(margin = margin(t = -8)))+
    labs(fill = "Contact niche")+
    guides(fill = guide_legend(nrow = 4, keywidth = 1,keyheight = 1))+
    theme(plot.margin = margin(5.5, 4, 15, 10)) 
  if (i != "B0: TN") {p <- p +theme(axis.title.y = element_blank())}
  pl[[i]] = p
}
p1 = ggpubr::ggarrange(plotlist = pl, common.legend = FALSE, nrow = 1, legend = "none") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
plot_save_tiff(p1, "CN_Prop_Plot_Pre", 6,10)


for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_response_boxplot(ad_here$obs, "tma_sid", "Fibroblasts-rich niche","spatial_niche")+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  plot_save_tiff(p, paste0("SN_Pre_Fibro_", i), 6,6)
}

for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_response_boxplot(ad_here$obs, "tma_sid", "Fibroblasts-rich niche","contact_niche")
  plot_save_tiff(p, paste0("CN_Pre_Fibro_", i), 6,6)
}


###### Plasma enrichment ######
enrichment_df <- data.frame(
  cell_type = c(
    "Plasma cells",
    "B cells",
    "T cells",
    "Endothelial cells",
    "Myogenic cells",
    "Mast cells",
    "Macrophages",
    "Fibroblasts",
    "Neuronal cells",
    "Dendritic cells",
    "Neutrophils",
    "Tumor",
    "Epithelial cells"
  ),
  fold_enrichment = c(
    11.959760,
    1.688326,
    1.510111,
    1.161711,
    1.124090,
    1.097226,
    1.059553,
    0.983491,
    0.751844,
    0.577286,
    0.206498,
    0.100740,
    0.000000
  )
)

# Order from highest to lowest enrichment
enrichment_df <- enrichment_df[order(enrichment_df$fold_enrichment, decreasing = TRUE), ]
enrichment_df$cell_type <- factor(
  enrichment_df$cell_type,
  levels = enrichment_df$cell_type[order(enrichment_df$fold_enrichment)]
)

enrichment_df_t <- data.frame(
  cell_type = c(
    "CD8_T_N",
    "CD8_T_Q",
    "CD4- CD8-",
    "CD8_T_EM",
    "CD8_T_TR",
    "CD4+ CD8-",
    "CD8_T_DYS_1",
    "CD8_T_CYT_DYS",
    "CD4+ CD8+",
    "CD8_T_SL",
    "CD8_T_DYS_2"
  ),
  fold_enrichment = c(
    10.208720,
    2.416279,
    1.923556,
    1.586605,
    1.516254,
    1.510231,
    1.310910,
    1.044141,
    1.019091,
    0.979617,
    0.030727
  )
)
enrichment_df_t <- enrichment_df_t[order(enrichment_df_t$fold_enrichment, decreasing = TRUE), ]
enrichment_df_t$cell_type <- factor(
  enrichment_df_t$cell_type,
  levels = enrichment_df_t$cell_type[order(enrichment_df_t$fold_enrichment)]
)

library(ggplot2)


p1 = ggplot(enrichment_df,
       aes(x = fold_enrichment,
           y = cell_type,
           fill = fold_enrichment > 1)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#4C78A8",
                               "FALSE" = "#E45756"),
                    guide = "none") +
  theme_Publication(base_size = 14)+
  labs(
    x = "Enrichment ratio",
    y = NULL
  )


p2 = ggplot(enrichment_df_t,
            aes(x = fold_enrichment,
                y = cell_type,
                fill = fold_enrichment > 1)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             colour = "black",
             linewidth = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#4C78A8",
                               "FALSE" = "#E45756"),
                    guide = "none") +
  theme_Publication(base_size = 14)+
  labs(
    x = "Enrichment ratio",
    y = NULL
  )

plot_save_tiff(p1, "ER_All", 6,5)
plot_save_tiff(p2, "ER_T", 5,5)





###### LR Interactions ######
inter_cols = c("NR-specific" = response_colors[["Non-response"]],
               "R-specific" = response_colors[["Response"]],
               "Shared" = "grey70")
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P4/")
lr_files <- list.files(path = data_dir, pattern = "^LR", full.names = TRUE)
df_all <- do.call(rbind, lapply(list.files(path = data_dir, pattern = "^LR", full.names = TRUE), read.csv))
df_all[df_all$interaction_type == "Non_Responder-specific", 'interaction_type'] = "NR-specific"
df_all[df_all$interaction_type == "Responder-specific", 'interaction_type'] = "R-specific"

shortlisted = c('Fibroblasts-Tumor | COL4A1-CD44',
                'Fibroblasts-Tumor | COL4A1-SDC1',
                'Fibroblasts-T Cells (Q2: 250-500) | COL1A1-SDC1',
                'Fibroblasts-T Cells (Q2: 250-500) | COL1A2-SDC1',
                'T Cells (Q2: 250-500)-Fibroblasts | CCL5-SDC1',
                'T Cells (Q3: 500-800)-Fibroblasts | CCL5-SDC1',
                'Fibroblasts-T Cells (Q3: 500-800) | MIF-CD74',
                'Fibroblasts-T Cells (Q3: 500-800) | THBS2-ITGA4',
                'Fibroblasts-T Cells (Q3: 500-800) | MIF-CXCR4',
                'Fibroblasts-T Cells (Q3: 500-800) | VIM-CD44',
                'Fibroblasts-T Cells (Q3: 500-800) | MMP1-CD44',
                'Fibroblasts-T Cells (Q2: 250-500) | TIMP2-CD44',
                'Fibroblasts-T Cells (Q1: 0-250) | LUM-ITGB1',
                'Fibroblasts-T Cells (Q1: 0-250) | PLAU-LRP1',
                'Fibroblasts-T Cells (Q1: 0-250) | LAMA4-CD44',
                'Fibroblasts-Tumor | COL6A2-ITGA3'
                )

df_short = df_all[df_all$pair_ct_genes %in% shortlisted,]

df_short = df_short%>%
  mutate(pair_ct_genes = factor(pair_ct_genes, levels = unique(pair_ct_genes)))
library(stringr)
df_short$pair_ct_genes <- str_replace_all(
  df_short$pair_ct_genes,
  c(
    "Q0: TN"       = "B0: TN",
    "Q1: 0-250"   = "B1:0–30",
    "Q2: 250-500" = "B2:30–60",
    "Q3: 500-800" = "B3:60–96",
    "Q4: 800+"  = "96+ µm"
  )
)
df_short = df_short %>% arrange(score_diff)%>%
  mutate(pair_ct_genes = factor(pair_ct_genes, levels = unique(pair_ct_genes)))

p1 = plot_bubble(df_short, "score_diff", "pair_ct_genes", "lr_means_NonResp", "interaction_type",
                 inter_cols, c("Bubble Plot", expression(Delta~"LR Score (NR - R)"), "Ligand-Receptor Pairs", "Interaction Strength in NR", "Interaction Type"))
plot_save_tiff(p1, "LR_Plot", 6,10)

###### Spatial Reconstruction Figure ######

data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P2/")
nm = "all_adata_tcells_subtypes_cd8s"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$T_subtypes_l1 = recode(as.character(adata$obs$T_subtypes_l1), !!!recode_ct_map)
adata$obs$T_subtypes_l1 = factor(adata$obs$T_subtypes_l1)
adata$obs$ct = recode(as.character(adata$obs$ct), !!!recode_ct_map)
adata$obs$ct = factor(adata$obs$ct)

p1 = spatial_plot_6k_ii(adata, 1, "tma_sid", "T_subtypes_l1", 
                        unlist(c(celltype_cols[c("Tumor", "B cells", "Plasma cells")], cd8_l1)), "TMA1_32")+
  labs(fill = "Cell type")+theme(legend.title = element_text(hjust = 0.5))
plot_save(p1, "TMA1_32", 7,9)


p1 = spatial_plot_6k_ii(adata, 1, "tma_sid", "ct", 
                        unlist(c(celltype_cols[c("Tumor", "T cells", "Fibroblasts")])), "TMA1_23")+
  theme(legend.position = "bottom")+
  labs(fill = "Cell Type")+theme(legend.title = element_text(hjust = 0.5))
plot_save(p1, "TMA1_23", 7,7)



###### end ######
