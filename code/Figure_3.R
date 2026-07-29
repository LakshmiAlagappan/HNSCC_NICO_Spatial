source("/Users/lakshmi_nccs/Desktop/NCCS/Projects/For_Publication/HNSCC_CosMx6k/helpers/design_elements.R")
source("/Users/lakshmi_nccs/Desktop/NCCS/Projects/For_Publication/HNSCC_CosMx6k/helpers/plotting_functions.R")

plot_results_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/Paper_Figures/Figure_2"
library('anndata')
library('patchwork')

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

###### CD4 vs CD8 Resolving of DN ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P1/")
nm = "tcells_harm_resolved"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
#adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

p1 = plot_umap_list(adata, c1=1, c2=2, quant = FALSE, ptsize = 1, pal = c(tsub_l1),
                    color_group_list=list(adata$obs$tcell_type, adata$obs$tcell_type_resolved),
                    c("tcell_type", "tcell_type_resolved"), " ")
plot_save_tiff(p1[[1]] + theme(plot.title = element_blank())+
            ggplot2::guides(colour=ggplot2::guide_legend(nrow = 2, override.aes = list(shape = 15, size = 4)))+
            labs(color = "T lineage")+theme(legend.position = "none"),
          "T_umap_Tcellsubtype", 6,6)
plot_save_tiff(p1[[2]] + theme(plot.title = element_blank())+
            ggplot2::guides(colour=ggplot2::guide_legend(nrow = 2, override.aes = list(shape = 15, size = 4)))+
            labs(color = "T lineage")+theme(legend.position = "none"),
          "T_umap_Tcellsubtype_resolved", 6,6)

p2 = plot_stacked_bar_percentage(adata$obs, "Patient", "tcell_type", c(tsub_l1), FALSE, FALSE, c("Patient", ""))+
  ggplot2::theme(legend.position = "bottom", plot.title = element_blank(), 
                 axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = -2))) + 
  guides(fill = guide_legend(nrow = 2, keywidth = 0.6,keyheight = 0.6))+
  labs(fill = "T lineage")+theme(legend.position = "none")
  
plot_save_tiff(p2, "T_stack_pat_lin_tcelltype", 5,7)

p2 = plot_stacked_bar_percentage(adata$obs, "TumorType", "tcell_type", c(tsub_l1), FALSE, FALSE, c("TumorType", ""))+
  ggplot2::theme(legend.position = "bottom", plot.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(nrow = 2, keywidth = 0.6,keyheight = 0.6))+
  labs(fill = "T lineage")
plot_save_tiff(p2, "T_stack_tt_lin_tcelltype", 6,5)

p2 = plot_stacked_bar_percentage(adata$obs, "Patient", "tcell_type_resolved", c(tsub_l1), FALSE, FALSE, c("Patient", ""))+
  ggplot2::theme(legend.position = "right", plot.title = element_blank(), 
                 axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = -2))) + 
  guides(fill = guide_legend(ncol = 1, keywidth = 0.6,keyheight = 0.6))+theme(legend.title = element_text(hjust = 0.5))+
  labs(fill = "T lineage")
plot_save_tiff(p2, "T_stack_pat_lin_tcelltype_res", 5,8)
p2 = plot_stacked_bar_percentage(adata$obs, "TumorType", "tcell_type_resolved", c(tsub_l1), FALSE, FALSE, c("TumorType", ""))+
  ggplot2::theme(legend.position = "bottom", plot.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(nrow = 2, keywidth = 0.6,keyheight = 0.6))
plot_save_tiff(p2, "T_stack_tt_lin_tcelltype_res", 6,5)+
  labs(fill = "T lineage")


df = adata$obs %>% select (tcell_type, tcell_type_resolved)
colnames(df) = c("Original", "Resolved")
df = df %>% pivot_longer(cols = c("Original", "Resolved"), names_to = "Type", values_to = "T lineage")
p1 = plot_stacked_bar_percentage(df, "Type", "T lineage", c(tsub_l1), TRUE, FALSE, c(" ", ""))+
  ggplot2::theme(legend.position = "bottom", plot.title = element_blank(),
                 axis.title.x = element_blank())+
              #axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(nrow = 2, keywidth = 0.6,keyheight = 0.6))  
plot_save_tiff(p1, "T_stack_type_lin", 6,4)





###### CD8 Subtypes ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P2/")
nm = "cd8_tcells_adata_res_1.5_annot"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
#adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

p1 = plot_umap(dat = adata, 1, 2, FALSE, 1, cd8_l1, adata$obs$T_subtypes_l1, rep(1, length(adata$obs$ct)),
               c(paste0("PN-Tum"), "CD4- CD8+\nsubtypes", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(nrow = 3, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "none", plot.title = element_blank())+
  theme(legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "CD8_umap_l1", 5,6)

 p1 = plot_umap(dat = adata, 1, 2, FALSE, 1, c20, adata$obs$T_subtypes_l2, rep(1, length(adata$obs$ct)),
               c(paste0("PN-Tum"), "CD4- CD8+\nsubtypes", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(ncol = 2, override.aes = list(shape = 15, size = 5)))+
  ggplot2::theme(legend.position = "right", plot.title = element_blank())+
  theme(legend.text = element_text(size = 18), axis.text =element_text(size = 18))+
  theme(legend.title = element_text(hjust = 0.5, size = 18))
plot_save_tiff(p1, "CD8_umap_l2", 8,23)

p2 = plot_stacked_bar_percentage(adata$obs, "Patient", "T_subtypes_l1", c(cd8_l1), FALSE, FALSE, c("Patient", " "))+
  ggplot2::theme(legend.position = "right", plot.title = element_blank(), 
                 axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14,margin = margin(t = -8)))+
  guides(fill = guide_legend(ncol = 1, keywidth = 1,keyheight = 1)) + 
  labs(fill = "CD4- CD8+\nsubtypes") + theme(legend.title = element_text(hjust = 0.5))
plot_save_tiff(p2 + theme(legend.position = "none"), "CD8_stack_l1", 5,6)

plot_save_tiff(p2 + 
                 guides(fill = guide_legend(nrow = 1, keywidth = 1,keyheight = 1)) + 
                 theme(legend.position = "bottom"), "CD8_stack_l1_legend", 6,14)

library(scales)
gns = c("CD3D","CD3E","CD3G","CD247",#'TRAC', "TRBC1", "TRBC2","CD8A", "CD8B", #T Cells
        "G0S2", "PTEN", #CD_T_SL
        "BTG1", #CD_T_Q
        "BTG2", "LEF1",#CD_T_N
        "IL7R", "CCR7", "SELL", "CD27", "CD28", "CD69",  "GZMK", "GZMM",  #CD_T_EM
        "BMI1", "CCL5", "GZMA", "GNLY", "GZMH", "TIGIT",  "CTLA4", "LAG3", "PRF1", "GZMB","ENTPD1", "HAVCR2","TOX", 
        "LGALS1", "FOS", "STAT2", 
        "JUN","LGALS3", "CD44", 
        "BCL6","ITGAE", "ITGA7", "NCAM1"
        #,"BMI1", "FOXP1", "PRDM1" Add this for l2 if needed
)

data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P2/")
expr_summary <- read.csv(paste0(data_dir, "cd8_tcells_deg_dotplot_summary.csv"))
expr_summary$cluster <- factor(expr_summary$cluster, levels = unique(expr_summary$cluster))
expr_summary$gene <- factor(expr_summary$gene, levels = unique(expr_summary$gene))
expr_summary <- expr_summary %>% filter(gene %in% gns)
expr_summary <- expr_summary %>%
  group_by(gene) %>%
  mutate(logFC_scaled = rescale(logFC, to = c(-5, 5))) %>%
  ungroup()

expr_summary$gene <- factor(expr_summary$gene, levels = gns)
expr_summary$cluster <- factor(expr_summary$cluster, levels = c("CD8_T_SL", "CD8_T_Q", "CD8_T_N", 
                                                                "CD8_T_EM", 
                                                                "CD8_T_CYT_DYS", "CD8_T_DYS_1",
                                                                "CD8_T_DYS_2", "CD8_T_TR"))


p1 = plot_dotplot(expr_summary, "gene", "cluster", "logFC_scaled", "pct_expr", c("logFC", "% expressed"))
plot_save_tiff(p1, "dotplot_cd8", 7,12)

###### CD4 Subtypes ######

data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P2/")
nm = "cd4_tcells_adata_res_1.5_annot"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
#adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

p1 = plot_umap(dat = adata, 1, 2, FALSE, 1.2, cd4_l1, adata$obs$T_subtypes_l1, rep(1, length(adata$obs$ct)),
               c(paste0("PN-Tum"), "CD4+ CD8-\nsubtypes", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(nrow = 3, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "none", plot.title = element_blank()) + theme(legend.title = element_text(hjust = 0.5))
plot_save_tiff(p1, "CD4_umap_l1", 6,7)

p1 = plot_umap(dat = adata, 1, 2, FALSE, 1.2, c20, adata$obs$T_subtypes_l2, rep(1, length(adata$obs$ct)),
               c(paste0("PN-Tum"), "CD4+ CD8-\nsubtypes", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(ncol = 2, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "right", legend.text = element_text(size = 18), axis.text =element_text(size = 18), 
                 plot.title = element_blank()) + theme(legend.title = element_text(hjust = 0.5, size = 18))
plot_save_tiff(p1, "CD4_umap_l2", 8,20)

p2 = plot_stacked_bar_percentage(adata$obs, "Patient", "T_subtypes_l1", c(cd4_l1), FALSE, FALSE, c("Patient", " "))+
  ggplot2::theme(legend.position = "right", plot.title = element_blank(), 
                 axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14,margin = margin(t = -8)))+
  guides(fill = guide_legend(ncol = 1, keywidth = 1,keyheight = 1)) + 
  labs(fill = "CD4+ CD8-\nsubtypes") + theme(legend.title = element_text(hjust = 0.5))
plot_save_tiff(p2, "CD4_stack_l1", 5,9)

###### Response based broad alluvial plots for CD8s Resp vs non Resp ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P2/")
nm = "cd8_tcells_adata_res_1.5_annot_resp"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

resp_pp  = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
resp_pp = resp_pp%>% filter(Response %in% c("Response"))
p1 = plot_alluvium(resp_pp, "TumorType", "T_subtypes_l1", cd8_l1, c("Response", "Cell type", "Cell type"))
p1 = p1 + coord_cartesian(clip = "off") + annotate("text",x = 0.55,y = 0.5,label = "Subtype proportion (%)",
                                                   angle = 90,hjust = 0.5,vjust = 0.5,size = 4)
plot_save_tiff(p1, "CD8_alluvial_resp", 5,6)

nresp_pp  = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
nresp_pp = nresp_pp%>% filter(Response %in% c("Non-response"))
p2 = plot_alluvium(nresp_pp, "TumorType", "T_subtypes_l1", cd8_l1, c("Non-response", "Cell type", "Cell type"))
p2 = p2 + coord_cartesian(clip = "off") + annotate("text",x = 0.55,y = 0.5,label = "Subtype proportion (%)",
                                                   angle = 90,hjust = 0.5,vjust = 0.5,size = 4)
plot_save_tiff(p2, "CD8_alluvial_non_resp", 5,6)


data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P2/")
nm = "cd4_tcells_adata_res_1.5_annot_resp"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

resp_pp  = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
resp_pp = resp_pp%>% filter(Response %in% c("Response"))
p1 = plot_alluvium(resp_pp, "TumorType", "T_subtypes_l1", cd4_l1, c("Response", "Cell type", "Cell type"))
p1 = p1 + coord_cartesian(clip = "off") + annotate("text",x = 0.55,y = 0.5,label = "Subtype proportion (%)",
                                                   angle = 90,hjust = 0.5,vjust = 0.5,size = 4)
plot_save_tiff(p1, "CD4_alluvial_resp", 6,7)

nresp_pp  = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
nresp_pp = nresp_pp%>% filter(Response %in% c("Non-response"))
p2 = plot_alluvium(nresp_pp, "TumorType", "T_subtypes_l1", cd4_l1, c("Non-response", "Cell type", "Cell type"))
p2 = p2 + coord_cartesian(clip = "off") + annotate("text",x = 0.55,y = 0.5,label = "Subtype proportion (%)",
                                                   angle = 90,hjust = 0.5,vjust = 0.5,size = 4)
plot_save_tiff(p2, "CD4_alluvial_non_resp", 6,7)


###### Categ 4 Boxplots ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P2/")
nm = "cd8_tcells_adata_res_1.5_annot_resp"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

df = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
for (celltype_of_interest in unique(adata$obs$T_subtypes_l1)) {
  p1 = plot_g4_paired_boxplot(df, celltype_of_interest, "T_subtypes_l1")
  p1 = p1 + ggplot2::theme(legend.position = "none")
  plot_save_tiff(p1, paste0("CD8_boxplot_g4_", celltype_of_interest), 5,4)
}

###### Categ 4 Boxplots ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P2/")
nm = "cd4_tcells_adata_res_1.5_annot_resp"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

df = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
for (celltype_of_interest in unique(adata$obs$T_subtypes_l1)) {
  p1 = plot_g4_paired_boxplot(df, celltype_of_interest, "T_subtypes_l1")
  p1 = p1 + ggplot2::theme(legend.position = "none")
  plot_save_tiff(p1, paste0("CD4_boxplot_g4_", celltype_of_interest), 5,4)
}



###### Delta Tile Plot ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P2/")
nm = "cd8_tcells_adata_res_1.5_annot_resp"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

df = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
ct_col = "T_subtypes_l1"

frac_df <- df %>%
  group_by(Patient, TumorType) %>%
  mutate(total_cells = n()) %>%
  ungroup() %>%
  group_by(Patient, TumorType, Response, !!sym(ct_col)) %>%
  summarise(
    celltype_cells = n(),
    total_cells = first(total_cells),
    fraction = celltype_cells / total_cells,
    .groups = "drop"
  ) %>%
  rename(celltype = !!sym(ct_col))

paired_patients <- frac_df %>%
  distinct(Patient, TumorType) %>%
  count(Patient) %>%
  filter(n == 2) %>%        # has both Pre and Post
  pull(Patient)

delta_df <- frac_df %>%
  filter(Patient %in% paired_patients) %>%     # 👈 key line
  select(Patient, Response, celltype, TumorType, fraction) %>%
  pivot_wider(
    names_from = TumorType,
    values_from = fraction,
    values_fill = list(fraction = 0)            # safe now
  ) %>%
  mutate(
    delta = `Post-NIVO primary` - `Pre-NIVO primary`
  )

heat_df <- delta_df %>%
  select(Patient, celltype, delta) %>%
  pivot_wider(
    names_from = Patient,
    values_from = delta
  )

heat_mat <- heat_df %>%
  column_to_rownames("celltype") %>%
  as.matrix()

library(ComplexHeatmap)

patient_annot <- delta_df %>%
  select(Patient, Response) %>%
  distinct() %>%
  arrange(match(Patient, colnames(heat_mat)))

ha_col <- HeatmapAnnotation(
  Patient  = patient_annot$Patient,
  Response = patient_annot$Response,
  col = list(
    Patient  = patient_colors,   # named vector: names = patient IDs
    Response = response_colors   # named vector: names = response levels
  )
)

ht <- Heatmap(
  heat_mat,
  name = "Delta",
  top_annotation = ha_col,
  col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
  rect_gp = gpar(col = "black", lwd = 0.3),  # <-- black outline
  na_col = "grey90",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_gap = unit(0, "mm"),  
  column_split = patient_annot$Response,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90
)

tiff(paste0(plot_results_dir, "/heatmap.tiff"), width = 2200, height = 1200, res = 300)
draw(ht)
dev.off()




