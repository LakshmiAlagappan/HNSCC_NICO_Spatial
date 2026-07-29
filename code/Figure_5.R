source("/Users/lakshmi_nccs/Desktop/NCCS/Projects/For_Publication/HNSCC_CosMx6k/helpers/design_elements.R")
source("/Users/lakshmi_nccs/Desktop/NCCS/Projects/For_Publication/HNSCC_CosMx6k/helpers/plotting_functions.R")

plot_results_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/Paper_Figures/Figure_4_mm"
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

old_niche = c("Epithelial_rich_niche","Fibroblasts_rich_niche","Macrophages_rich_niche", "Neuronal_rich_niche","Neutrophils_rich_niche", "Plasma_rich_niche","T_rich_niche","Tumor_rich_niche")
new_niche = c("Epithelial-rich niche","Fibroblasts-rich niche","Macrophages-rich niche", "Neuronal-rich niche", "Neutrophils-rich niche", "Plasma-rich niche","T-rich niche","Tumor-rich niche")
recode_niche_map <- setNames(new_niche,old_niche)

old_con_niche <- c("B Cells_rich_niche", "Dentritic Cells_rich_niche","Endothelial Cells_rich_niche","Epithelial Cells_rich_niche",
                   "Fibroblasts_rich_niche","Macrophages_rich_niche","Mast Cells_rich_niche","Myogenic Cells_rich_niche","Neuronal Cells_rich_niche",
                   "Neutrophils_rich_niche","Plasma Cells_rich_niche","T Cells_rich_niche","Tumor_rich_niche", "No-contact")
new_con_niche <- c("B-rich niche","Dendritic-rich niche","Endothelial-rich niche","Epithelial-rich niche",
                   "Fibroblasts-rich niche","Macrophages-rich niche","Mast-rich niche","Myogenic-rich niche","Neuronal-rich niche",
                   "Neutrophils-rich niche","Plasma-rich niche","T-rich niche","Tumor-rich niche","No-contact niche")
recode_con_niche_map <- setNames(new_con_niche,old_con_niche)

old_scc_bin_levels = c("Q1: 0-150", "Q2: 150-300", "Q3: 300-500", "Q4: 500+")
scc_bin_levels = c("B1: 0-150", "B2: 150-300", "B3: 300-500", "B4: 500+")
scc_bin_levels = c("B1: 0-18", "B2: 18-36", "B3: 36-60", "B4: 60+")
recode_scc_bin_map <- setNames(scc_bin_levels,old_scc_bin_levels)

comet_old_ct_levels =c("B Cells","CD4+","CD8+","Dendritic Cells","Endothelial Cells","Fibroblasts",
                       "Macrophages","Mast Cells","Myeloid Cells","Tumor")
comet_ct_levels =c("B cells","CD4+ T cells","CD8+ T cells","Dendritic cells","Endothelial cells","Fibroblasts",
                   "Macrophages","Mast cells","Myeloid cells","Tumor")
recode_comet_ct_map = setNames(comet_ct_levels,comet_old_ct_levels)


###### Landscape of differences:  Advancing Front Responder ######
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

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO adv front"), ]
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
  mutate(diff = `Post-NIVO adv front` - `Pre-NIVO primary`)

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
  scale_x_continuous(breaks = seq(0, 240, by = 30), expand = c(0, 0)) +   # 👈 removes padding
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.line = element_blank(),           # ❌ remove axis lines
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2, "pt"),     # smaller tick spacing
    plot.margin = margin(10, 10, 10, 10)       # tighter margins
  ) +
  theme(legend.title = element_text(hjust = 0.5))

plot_save_tiff(p1, "CD8_landscape_Resp_AdvF", 5,11)


###### Landscape of differences:  Advancing Front Non_Responder ######
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

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO adv front"), ]
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
  mutate(diff = `Post-NIVO adv front` - `Pre-NIVO primary`)

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
    y = "CD8 Spatial enrichment score \n Non-response", #Scaled (Post − Pre) \n (Non_Responders)
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

plot_save_tiff(p1, "CD8_landscape_Non_Resp_AdvF", 5,11)








###### Enrichment Bars:  Advancing Front Responder######
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

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO adv front"), ]
df = adata$obs

df <- df %>%
  mutate(spatial_bin = cut(
    Distance_to_Tumor,
    breaks = c(0, 500, 5329.787),
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
plot_df$group4 <- factor(plot_df$group4,levels = c("Near::Pre-NIVO primary","Near::Post-NIVO adv front",
                                                   "Far::Pre-NIVO primary","Far::Post-NIVO adv front"))

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
      "Near::Post-NIVO adv front" = "#E91E63",
      "Far::Pre-NIVO primary" = "#142B7A",
      "Far::Post-NIVO adv front" = "#9e1e49")) +
  scale_y_continuous(labels = abs) +
  labs(y = "Median proportion of subtype \n per FOV - Response",x = "CD4- CD8+ subtype",fill = "Status") +
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

plot_save_tiff(p1 + theme(legend.position = "None"), "CD8_enrich_Responders", 4,6)
plot_save_tiff(p1 + theme(legend.position = "bottom"), "CD8_enrich_Responders_legend", 6,12)


###### Enrichment Bars:  Advancing Front Non_Responder######
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

adata = adata[adata$obs$TumorType %in% c("Pre-NIVO primary", "Post-NIVO adv front"), ]
df = adata$obs

df <- df %>%
  mutate(spatial_bin = cut(
    Distance_to_Tumor,
    breaks = c(0, 500, 5329.787),
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
plot_df$group4 <- factor(plot_df$group4,levels = c("Near::Pre-NIVO primary","Near::Post-NIVO adv front",
                                                   "Far::Pre-NIVO primary","Far::Post-NIVO adv front"))

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
      "Near::Post-NIVO adv front" = "#E91E63",
      "Far::Pre-NIVO primary" = "#142B7A",
      "Far::Post-NIVO adv front" = "#9e1e49")) +
  scale_y_continuous(labels = abs) +
  labs(y = "Median proportion of subtype \n per FOV - Non-response",x = "CD4- CD8+ subtype",fill = "Status") +
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

plot_save_tiff(p1+ theme(legend.position = "None"), "CD8_enrich_Non_Responders", 4,6)





###### Fibroblasts Niche Advancing Front ######
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

adata$obs$spatial_niche = recode(as.character(adata$obs$spatial_niche), !!!recode_niche_map)
adata$obs$spatial_niche = factor(adata$obs$spatial_niche, levels = new_niche)
adata$obs$contact_niche = recode(as.character(adata$obs$contact_niche), !!!recode_con_niche_map)
adata$obs$contact_niche = factor(adata$obs$contact_niche, levels = new_con_niche)
adata$obs$spatial_bin = recode(as.character(adata$obs$spatial_bin), !!!recode_bin_map)
adata$obs$spatial_bin = factor(adata$obs$spatial_bin, levels = bin_levels)

adata = adata[adata$obs$TumorType %in% c("Post-NIVO adv front"), ]$copy()
adata$obs$Response <- recode(
  adata$obs$Response,
  "Response"  = "Post-Adv-R",
  "Non-response" = "Post-Adv-NR"
)
pl = list()
for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_stacked_bar_percentage(ad_here$obs, "Response", "spatial_niche", c(niche_colors_spat), TRUE, FALSE, c(i, " "))+
    ggplot2::theme(legend.position = "bottom") + 
    theme(plot.title = element_blank(),
          plot.text = element_text(size = 16),
          axis.text.x = element_text(size = 12, margin = margin(t = -10), angle = 20, hjust = 1),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))+
    labs(fill = "Spatial niche")+
    guides(fill = guide_legend(nrow = 2, keywidth = 1,keyheight = 1))+
    theme(plot.margin = margin(5.5, 4, 15, 10)) 
  if (i != "B0: TN") {p <- p +theme(axis.title.y = element_blank())}
  pl[[i]] = p
}
p1 = ggpubr::ggarrange(plotlist = pl, common.legend = TRUE, nrow = 1, legend = "bottom") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
plot_save_tiff(p1, "SN_Prop_Plot_Adv", 6,10)

pl = list()
for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_stacked_bar_percentage(ad_here$obs, "Response", "contact_niche", c(niche_colors_con), TRUE, FALSE, c(i, " "))+
    ggplot2::theme(legend.position = "bottom") + 
    theme(plot.title = element_blank(),
          axis.text.x = element_text(size = 10, margin = margin(t = -8)),
          axis.title = element_text(size = 10))+
    labs(fill = "Contact niche")+
    guides(fill = guide_legend(nrow = 4, keywidth = 1,keyheight = 1))+
    theme(plot.margin = margin(5.5, 4, 15, 10)) 
  if (i != "B0: TN") {p <- p +theme(axis.title.y = element_blank())}
  pl[[i]] = p
}
p1 = ggpubr::ggarrange(plotlist = pl, common.legend = TRUE, nrow = 1, legend = "bottom") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
plot_save_tiff(p1, "CN_Prop_Plot_Adv", 8,10)


for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_response_boxplot(ad_here$obs, "tma_sid", "Fibroblasts-rich niche","spatial_niche")
  plot_save_tiff(p, paste0("SN_Post-Adv_Fibro_", i), 6,6)
}

for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_response_boxplot(ad_here$obs, "tma_sid", "Fibroblasts-rich niche","contact_niche")
  plot_save_tiff(p, paste0("CN_Post-Adv_Fibro_", i), 6,6)
}


###### Advancing Front Heatmap of selected Genes ######
library(viridis)
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
adata$obs$Distance_to_Tumor = adata$obs$Distance_to_Tumor *0.120
adata_pre = adata[adata$obs$TumorType %in% c("Pre-NIVO primary"), ]
adata_pre_R = adata_pre[adata_pre$obs$Response %in% c("Response"), ]
adata_pre_NR = adata_pre[adata_pre$obs$Response %in% c("Non-response"), ]

adata_post = adata[adata$obs$TumorType %in% c("Post-NIVO adv front"), ]
adata_post_R = adata_post[adata_post$obs$Response %in% c("Response"), ]
adata_post_NR = adata_post[adata_post$obs$Response %in% c("Non-response"), ]


gns = c("KRT5", "LAG3", "LGALS3", "HLA-DRA", "HLA-DRB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2",
        "HLA-DQB1", "HLA-DQB2")

gns = c('KRT5', 'IGKC', 'COL1A1', 'PECAM1', 'CD68',
        "LGALS3","LGALS1","CXCL8","GNLY","PRF1",
        "IL7R",  "CTLA4", "LAG3", "CXCL13", "TOX", "TIGIT", 
        "GZMH", "CCL5", "CXCR4","KLF2","BMI1","LEF1" )
p1 = plot_gene_distance_heatmap(adata_post_R, 
                                gns, 
                                cbreaks = c(c(0, 50, 100, 150),seq(200, 2500, by = 100), 
                                            c(3000, 4000, 5000, 6000))*0.120, 
                                max_distance = 6000*0.120, tit = "Post-Adv-R",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Heatmap_Dis_PostAdvR.tiff"), 
       plot =grid::grid.draw(p1),
       width = 4, height = 4, bg = "white")


p2 = plot_gene_distance_heatmap(adata_post_NR, 
                                gns, 
                                cbreaks = c(c(0, 50, 100, 150),seq(200, 2500, by = 100), 
                                            c(3000, 4000, 5000, 6000))*0.120, 
                                max_distance = 6000*0.120,tit = "Post-Adv-NR",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Heatmap_Dis_PostAdvNR.tiff"), 
       plot =grid::grid.draw(p2),
       width = 4, height = 4, bg = "white")






###### COMET Analysis UMAP, Stack, Markers ######
adata = ad$read_h5ad("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P5/comet_s1_adata_res_0.5_annot_2.h5ad") 
adata$obs$ct = adata$obs$ct_2
colnames(adata$obs)[2] = "CellId"
adata$obs$ct = recode(as.character(adata$obs$ct), !!!recode_comet_ct_map)
adata$obs$ct[adata$obs$ct == "Tumor" & adata$obs$TumorType == "Post-NIVO Normal"] = "Epithelial cells"
adata$obs$ct = factor(adata$obs$ct, levels = c(comet_ct_levels, "Epithelial cells"))

vn <- adata$var_names$to_list()
vn[vn == "LAG.3"] <- "LAG3"
vn[vn == "E.Cad"] <- "E-Cadherin"
adata$var_names <- vn

#adata$var_names[adata$var_names == "LAG.3"] <- "LAG3"
#adata$var_names[adata$var_names == "E.Cad"] <- "E-Cadherin"

p1 = plot_umap(dat = adata, 1, 2, FALSE, 0.5, celltype_cols, adata$obs$ct, rep(1, length(adata$obs$ct)),
               c(paste0("PN-Tum"), "Cell type", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(nrow = 4, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "none", plot.title = element_blank())
plot_save_tiff(p1, "Comet_umap_ct", 6,7)

p2 = plot_stacked_bar_percentage(adata$obs, "Patient", "ct", c(celltype_cols), FALSE, FALSE, c("Patient", " "))+
  labs(fill = "Cell type")+
  ggplot2::theme(legend.position = "right", legend.title.align = 0.5) + 
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14,margin = margin(t = -10)))+
  guides(fill = guide_legend(ncol = 1, keywidth = 1,keyheight = 1))
plot_save_tiff(p2, "Comet_stack_pat_ct", 6,12)
plot_save_tiff(p2+ggplot2::theme(legend.position = "bottom")+
                 guides(fill = guide_legend(nrow = 1, keywidth = 1,keyheight = 1)), "Comet_stack_pat_ct_legend", 6,20)


markers = c("CD3", "CD4", "CD8", "LAG3", 
            "CK", "E-Cadherin")
p1 = plot_marker_umap(adata, markers[1:6], 0.2)
plot_save_tiff(p1, paste0('comet_markers_m1'), 6, 8)

###### COMET Heatmap and Lollipop of Selected Markers ######
adata = ad$read_h5ad("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P5/comet_s1_cd8_dist_resp.h5ad") 
adata$obs$ct = adata$obs$ct_2
colnames(adata$obs)[2] = "CellId"
adata$obs$ct = recode(as.character(adata$obs$ct), !!!recode_comet_ct_map)
adata$obs$ct = factor(adata$obs$ct, levels = comet_ct_levels)

adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)
vn <- adata$var_names$to_list()
vn[vn == "LAG.3"] <- "LAG3"
vn[vn == "E.Cad"] <- "E-Cadherin"
adata$var_names <- vn

#adata$var_names[adata$var_names == "LAG.3"] <- "LAG3"
#adata$var_names[adata$var_names == "E.Cad"] <- "E-Cadherin"
adata$obs$Distance_to_Tumor = adata$obs$Distance_to_Tumor * 0.280
adata_pre = adata[adata$obs$TumorType %in% c("Pre-NIVO primary"), ]$copy()
adata_pre$obs$Response <- recode(
  adata_pre$obs$Response,
  "Response"  = "Pre-R",
  "Non-response" = "Pre-NR"
)
adata_pre_R = adata_pre[adata_pre$obs$Response %in% c("Pre-R"), ]$copy()
adata_pre_NR = adata_pre[adata_pre$obs$Response %in% c("Pre-NR"), ]$copy()


adata_post = adata[adata$obs$TumorType %in% c("Post-NIVO primary"), ]$copy()
adata_post$obs$Response <- recode(
  adata_post$obs$Response,
  "Response"  = "Post-R",
  "Non-response" = "Post-NR"
)
adata_post_R = adata_post[adata_post$obs$Response %in% c("Post-R"), ]$copy()
adata_post_NR = adata_post[adata_post$obs$Response %in% c("Post-NR"), ]$copy()

gns = c('CK', 'CD19', 'Podoplanin', "VISTA", "CD163", "CD45RO",  'CD3', 'CD4',
        "CD8","LAG3","GZMB",  "EOMES", "Ki67", "CD45RA","FoxP3","PD1","ICOS")

#gns = c("CK", "LAG3", "HLA.DR")
p1 = plot_gene_distance_heatmap(adata_pre_R, 
                                gns, 
                                cbreaks = c(seq(0, 500, by = 20))*0.280,
                                max_distance = 1000*0.280, tit = "Pre-R",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Comet_Heatmap_Dis_PreR.tiff"), 
       plot =grid::grid.draw(p1),
       width = 4, height = 4, bg = "white")

p1 = plot_gene_distance_heatmap(adata_post_R, 
                                gns, 
                                cbreaks = c(seq(0, 500, by = 20))*0.280,
                                max_distance = 1000*0.280, tit = "Post-R",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Comet_Heatmap_Dis_PostR.tiff"), 
       plot =grid::grid.draw(p1),
       width = 4, height = 4, bg = "white")

p1 = plot_gene_distance_heatmap(adata_pre_NR, 
                                gns, 
                                cbreaks = c(seq(0, 500, by = 20))*0.280,
                                max_distance = 1000*0.280, tit = "Pre-NR",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Comet_Heatmap_Dis_PreNR.tiff"), 
       plot =grid::grid.draw(p1),
       width = 4, height = 4, bg = "white")

p1 = plot_gene_distance_heatmap(adata_post_NR, 
                                gns, 
                                cbreaks = c(seq(0, 500, by = 20))*0.280,
                                max_distance = 1000*0.280, tit = "Post-NR",
                                distance_obs = "Distance_to_Tumor")
ggsave(paste0(plot_results_dir, "/Comet_Heatmap_Dis_PostNR.tiff"), 
       plot =grid::grid.draw(p1),
       width = 4, height = 4, bg = "white")


adata_post$obs$Response = factor(adata_post$obs$Response, levels = c("Post-NR", "Post-R"))
p1 = plot_mirrored_lollipop(adata_post, gns, "Response", 0.8, cbreaks = seq(0, 500, by = 20)*0.280,
                                  distance_obs = "Distance_to_Tumor", pal = g4_colors,
                            cc = c("Distance to tumor (µm)","", "Status", "Expression","COMET: Post-R vs Post-NR"))
  
plot_save_tiff(p1, "COMET_post", 6,10)

adata_pre$obs$Response = factor(adata_pre$obs$Response, levels = c("Pre-NR", "Pre-R"))
p1 = plot_mirrored_lollipop(adata_pre, gns, "Response", 0.8, cbreaks = seq(0, 500, by = 20)*0.280,
                            distance_obs = "Distance_to_Tumor", pal = g4_colors,
                            cc = c("Distance to tumor  (µm)","", "Status", "Expression","COMET: Pre-R vs Pre-NR"))

plot_save_tiff(p1, "COMET_pre", 6,10)

###### COMET Spatial Reconstruction ######
adata = ad$read_h5ad("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P5/comet_s1_adata_res_0.5_annot_2.h5ad") 
adata$obs$ct = adata$obs$ct_2
colnames(adata$obs)[2] = "CellId"
adata$obs$ct = recode(as.character(adata$obs$ct), !!!recode_comet_ct_map)
adata$obs$ct = factor(adata$obs$ct, levels = comet_ct_levels)


unique(adata$obs$coreID)
p = spatial_plot_comet_ii(adata, 1, "coreID", "ct", celltype_cols, "B4")
plot_save(p, paste0("Marker_TMA1_coreID_", i, "_celltype_spatialplot"), 6,8)



###### SCC UMAP ######
adata = ad$read_h5ad("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P6/all_adata_tcell_type_centroid.h5ad")
adata$obs$ct = recode(as.character(adata$obs$ct), !!!recode_ct_map)
adata$obs$ct = factor(adata$obs$ct)
patients <- table(
  adata$obs$Patient,
  adata$obs$TumorType
)
patients_keep <- rownames(patients)[
  rowSums(patients[, c("Baseline","Post-dose")] > 0) == 2
]

resp_table <- unique(adata$obs[, c("Patient", "Response")])
resp_table <- resp_table[resp_table$Patient %in% patients_keep, ]
drug_table <- unique(adata$obs[, c("Patient", "Treatment Drug")])
drug_table <- drug_table[drug_table$Patient %in% patients_keep, ]

merge(resp_table, drug_table, by = "Patient")

p1 = plot_umap(dat = adata, 1, 2, FALSE, 0.5, celltype_cols, adata$obs$ct, rep(1, length(adata$obs$ct)),
               c(paste0("PN-Tum"), "Cell type", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(nrow = 2, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "bottom", plot.title = element_blank(), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "scc_all_umap_ct", 6,6)

p2 = plot_stacked_bar_percentage(adata$obs, "Patient", "Response", c(c20), FALSE, FALSE, c("Patient", " "))+
  labs(fill = "Cell type")+
  ggplot2::theme(legend.position = "right", legend.title.align = 0.5) + 
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14,margin = margin(t = -10)))+
  guides(fill = guide_legend(ncol = 1, keywidth = 1,keyheight = 1))
plot_save_tiff(p2, "Comet_stack_pat_ct", 6,12)

adata = ad$read_h5ad("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P6/tcells_harm_resolved.h5ad")
p1 = plot_umap(dat = adata, 1, 2, FALSE, 0.5, tsub_l1, adata$obs$tcell_type_resolved, rep(1, length(adata$obs$ct)),
               c(paste0("PN-Tum"), "T lineage", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(nrow = 2, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "bottom", plot.title = element_blank(), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "scc_T_umap", 6,6)


adata = ad$read_h5ad("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P6/all_adata_tcell_type_centroid.h5ad")
adata$obs$ct = recode(as.character(adata$obs$ct), !!!recode_ct_map)
adata$obs$ct = factor(adata$obs$ct)
p1 = ggplot(adata[adata$obs$tma_sid == '1_7']$obs, aes(x = centroid_x, y = centroid_y, color = ct)) +
  geom_point(size = 1) +
  scale_color_manual(values = celltype_cols)+
  ggplot2::guides(colour=ggplot2::guide_legend(ncol = 1, override.aes = list(shape = 15, size = 4)))+
  theme_Publication(base_size = 14) +
  theme(legend.title = element_text(hjust = 0.5))+
  labs(
    x = "Centroid X",
    y = "Centroid Y",
    color = "Cell type"
  )+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
plot_save_tiff(p1, "scc_Recon", 5,6)

###### SCC Histogram #######
adata = ad$read_h5ad("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P6/cd8_tcells_dis_tn_tum.h5ad")
adata = adata[adata$obs$TumorType %in% c("Baseline", "Post-dose"),]
adata = adata[adata$obs$Response %in% c("Responder", "Non_Responder"),]

tab = table(adata$obs$Patient,adata$obs$TumorType)
valid_patients <- rownames(tab)[tab[, "Baseline"] > 0 & tab[, "Post-dose"] > 0]
adata <- adata[adata$obs$Patient %in% valid_patients, ]

p1 = plot_hist_lines_points(adata$obs, 
                            "Distance_to_Tumor", "TumorType", 50, c("Baseline" = "#4363d8","Post-dose" = "#3cb44b"), 
                            c(" ", "Distance to tumor", "Number of CD8+ CD4- T cells"))+
  theme(plot.title = element_blank()) + labs(color = "Biopsy timepoint") +
  theme(legend.position = c(0.7, 0.7), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "SCC_CD8_hist_spat", 5,5)

####### SCC Shells ######
adata = ad$read_h5ad("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P6/cd8_tcells_dis_tn_tum.h5ad")
adata = adata[adata$obs$TumorType %in% c("Baseline", "Post-dose"),]
adata = adata[adata$obs$Response %in% c("Responder", "Non_Responder"),]

tab = table(adata$obs$Patient,adata$obs$TumorType)
valid_patients <- rownames(tab)[tab[, "Baseline"] > 0 & tab[, "Post-dose"] > 0]
adata <- adata[adata$obs$Patient %in% valid_patients, ]
krt5_idx <- which(adata$var_names$to_list() == "KRT5")
remove_cells1 <- (
  (adata$obs$Distance_to_Tumor > 250) & 
    (adata$X[,krt5_idx] > 0)
)
adata <- adata[!remove_cells1, ]

bb = c(0,100,150,250, seq(275, 450, 25))
gns = c('KRT5', 'COL1A1', 'PECAM1', 'CD68',
        "LGALS3","LGALS1","CXCL8","GNLY","PRF1",
        "IL7R",  "CTLA4", "LAG3", "CXCL13", "TOX", "TIGIT", 
        "GZMH", "CCL5", "CXCR4","KLF2","BMI1","LEF1" )
gns = gns[gns %in% adata$var_names$to_list()]
adata$obs$TumorType = recode(adata$obs$TumorType, "Baseline" = "Pre-NIVO Primary", "Post-dose" = "Post-NIVO Primary")

adata_post = adata[adata$obs$TumorType %in% c("Post-NIVO Primary"),]
adata_post$obs$Distance_to_Tumor = adata_post$obs$Distance_to_Tumor*0.120
p1 = plot_spatialshells_wrapper(adata_post, c("KRT5", "LAG3", "CXCL13", "PRF1", "TOX"),1000,bb=bb*0.120, scaling = "none")
plot_save_tiff(p1[[1]], "SCC_multi", 9,9)

###### SCC Niche UMAP ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P6/")
nm = "cd8_tcells_spfp_200_niche"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$spatial_niche = recode(as.character(adata$obs$spatial_niche), !!!recode_niche_map)
adata$obs$spatial_niche = factor(adata$obs$spatial_niche, levels = new_niche)

p1 = plot_umap(dat = adata, 1, 2, FALSE, 0.5, niche_colors_spat, 
               adata$obs$spatial_niche, rep(1, length(adata$obs$spatial_niche)),
               c(paste0("PN-Tum"), "Spatial Niche", " "))+
  theme_Publication(base_size = 10)+
  ggplot2::theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(ncol = 1, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "right", plot.title = element_blank(), legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1+ggplot2::theme(legend.position = "None"), "scc_spatnich_umap", 3,3)

###### SCC Niche Proportion ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P6/")
nm = "cd8_spat_metrics"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$spatial_niche = recode(as.character(adata$obs$spatial_niche), !!!recode_niche_map)
adata$obs$spatial_niche = factor(adata$obs$spatial_niche, levels = new_niche)
adata$obs$spatial_bin = recode(as.character(adata$obs$spatial_bin), !!!recode_scc_bin_map)
adata$obs$spatial_bin = factor(adata$obs$spatial_bin, levels = scc_bin_levels)

adata = adata[!adata$obs$Response %in% c("Unknown"),]

adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)


adata = adata[adata$obs$TumorType %in% c("Post-dose"),]
adata$obs$Response <- recode(
  adata$obs$Response,
  "Response"  = "Post-dose-R",
  "Non-response" = "Post-dose-NR"
)
pl = list()
for (i in levels(adata$obs$spatial_bin)[1:3]){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_stacked_bar_percentage(ad_here$obs, "Response", "spatial_niche", niche_colors_spat, TRUE, FALSE, c(i, " "))+
    ggplot2::theme(legend.position = "bottom") + 
    theme(plot.title = element_blank(),
          plot.text = element_text(size = 14),
          #axis.text.x = element_text(size = 12, margin = margin(t = -8)),
          axis.text.x = element_text(size = 12, margin = margin(t = -10), angle = 10, hjust = 0.8),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))+
    labs(fill = "Spatial niche")+
    guides(fill = guide_legend(nrow = 2, keywidth = 1,keyheight = 1))+
    theme(plot.margin = margin(5.5, 4, 15, 10)) 
  if (i != "B1: 0-18") {p <- p +theme(axis.title.y = element_blank())}
  pl[[i]] = p
}
p1 = ggpubr::ggarrange(plotlist = pl, common.legend = TRUE, nrow = 1, legend = "bottom") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
plot_save_tiff(p1, "SCC_SN_Prop_Plot", 5,6)

 plot_save_tiff(p1, "SCC_SN_Prop_Plot_legend", 6,10)

for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_response_boxplot(ad_here$obs, "Patient", "Fibroblasts-rich niche","spatial_niche")
  plot_save_tiff(p, paste0("SCC_SN_Post_Fibro_", i), 6,6)
}








###### SCC SN #######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P6/")
nm = "cd8_spat_metrics"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$spatial_niche = recode(as.character(adata$obs$spatial_niche), !!!recode_niche_map)
adata$obs$spatial_niche = factor(adata$obs$spatial_niche, levels = new_niche)

adata$obs$spatial_bin = recode(as.character(adata$obs$spatial_bin), !!!recode_scc_bin_map)
adata$obs$spatial_bin = factor(adata$obs$spatial_bin, levels = scc_bin_levels)


adata = adata[adata$obs$TumorType %in% c("Post-dose"), ]$copy()

adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)
adata$obs$Response <- recode(
  adata$obs$Response,
  "Responder"  = "Post-R",
  "Non_Responder" = "Post-NR"
)
pl = list()
for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_stacked_bar_percentage(ad_here$obs, "Response", "spatial_niche", c20, TRUE, FALSE, c(i, " "))+
    ggplot2::theme(legend.position = "bottom") + 
    theme(plot.title = element_blank(),
          axis.text.x = element_text(margin = margin(t = -8)))+
    labs(fill = "Spatial niche")+
    guides(fill = guide_legend(nrow = 2, keywidth = 1,keyheight = 1))+
    theme(plot.margin = margin(5.5, 4, 15, 10), legend.title = element_text(hjust = 0.5)) 
  if (i != "B1") {p <- p +theme(axis.title.y = element_blank())}
  pl[[i]] = p
}
p1 = ggpubr::ggarrange(plotlist = pl, common.legend = TRUE, nrow = 1, legend = "bottom") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
plot_save_tiff(p1, "SN_Prop_Plot_Post", 5,10)

for (i in levels(adata$obs$spatial_bin)){
  ad_here = adata[adata$obs$spatial_bin == i]$copy()
  p = plot_response_boxplot(ad_here$obs, "tma_sid", "Fibroblasts_rich_niche","spatial_niche")+
    theme(legend.position = "none")
  plot_save_tiff(p, paste0("SCC_SN_Post_Fibro_", i), 5,3)
}










###### Explants Shells ######
data_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Explants/results/"
plot_results_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Explants/results/"
plot_results_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/Paper_Figures/Figure_4"
gns = c("CCL5", "CXCL8", "CXCL13", "GAL1", "GAL3", "KLF2", "LAG3", "PDPN", "PRF1")

#Note HN440, HN421, HN369 and HN389 are clinically validated
old_pat_levels = c("HN440", "HN431", "P13", "P15", "HN421", "P8")
pat_levels = c("P1", "P2", "P3", "P4", "P5", "P6")
recode_pat_map <- setNames(pat_levels,old_pat_levels)

i = "CXCL8"
adata_PRF1 = ad$read_h5ad(paste0(data_dir, "Explants_", i, ".h5ad"))
adata_PRF1$obs$Patient = recode(as.character(adata_PRF1$obs$Patient), !!!recode_pat_map)
p1 = plot_spatialshells_wrapper_explants(adata_PRF1, i,dist = "Distance_norm", bb=seq(0, 1, by = 0.1), scaling = "min-max")
plot_save_tiff(p1[[1]], paste0("Explants_", i), 6,6)

i = "CXCL13"
adata_PRF1 = ad$read_h5ad(paste0(data_dir, "Explants_", i, ".h5ad"))
adata_PRF1$obs$Patient = recode(as.character(adata_PRF1$obs$Patient), !!!recode_pat_map)
p1 = plot_spatialshells_wrapper_explants(adata_PRF1, i,dist = "Distance_norm", bb=seq(0, 1, by = 0.1), scaling = "min-max")
plot_save_tiff(p1[[1]], paste0("Explants_", i), 6,6)

i = "CCL5"
adata_PRF1 = ad$read_h5ad(paste0(data_dir, "Explants_", i, ".h5ad"))
adata_PRF1$obs$Patient = recode(as.character(adata_PRF1$obs$Patient), !!!recode_pat_map)
p1 = plot_spatialshells_wrapper_explants(adata_PRF1, i,dist = "Distance_norm", bb=seq(0, 1, by = 0.1), scaling = "min-max")
plot_save_tiff(p1[[1]], paste0("Explants_", i), 6,6)


i = "GAL1"
adata_PRF1 = ad$read_h5ad(paste0(data_dir, "Explants_", i, ".h5ad"))
p1 = plot_spatialshells_wrapper_explants(adata_PRF1, i,dist = "Distance_norm", bb=seq(0, 1, by = 0.1), scaling = "min-max")
plot_save_tiff(p1[[1]], paste0("Explants_", i), 8,8)

i = "GAL3"
adata_PRF1 = ad$read_h5ad(paste0(data_dir, "Explants_", i, ".h5ad"))
p1 = plot_spatialshells_wrapper_explants(adata_PRF1, i,dist = "Distance_norm", bb=seq(0, 1, by = 0.1), scaling = "min-max")
plot_save_tiff(p1[[1]], paste0("Explants_", i), 8,8)

i = "KLF2"
adata_PRF1 = ad$read_h5ad(paste0(data_dir, "Explants_", i, ".h5ad"))
adata_PRF1$obs$Patient = recode(as.character(adata_PRF1$obs$Patient), !!!recode_pat_map)
p1 = plot_spatialshells_wrapper_explants(adata_PRF1, i,dist = "Distance_norm", bb=seq(0, 1, by = 0.1), scaling = "min-max")
plot_save_tiff(p1[[1]], paste0("Explants_", i), 6,6)

i = "LAG3"
adata_PRF1 = ad$read_h5ad(paste0(data_dir, "Explants_", i, ".h5ad"))
adata_PRF1$obs$Patient = recode(as.character(adata_PRF1$obs$Patient), !!!recode_pat_map)
p1 = plot_spatialshells_wrapper_explants(adata_PRF1, i,dist = "Distance_norm", bb=seq(0, 1, by = 0.1), scaling = "min-max")
plot_save_tiff(p1[[1]], paste0("Explants_", i), 6,6)

i = "PRF1"
adata_PRF1 = ad$read_h5ad(paste0(data_dir, "Explants_", i, ".h5ad"))
adata_PRF1$obs$Patient = recode(as.character(adata_PRF1$obs$Patient), !!!recode_pat_map)
p1 = plot_spatialshells_wrapper_explants(adata_PRF1, i,dist = "Distance_norm", bb=seq(0, 1, by = 0.1), scaling = "min-max")
plot_save_tiff(p1[[1]], paste0("Explants_", i), 6,6)

