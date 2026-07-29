source("/Users/lakshmi_nccs/Desktop/NCCS/Projects/For_Publication/HNSCC_CosMx6k/helpers/design_elements.R")
source("/Users/lakshmi_nccs/Desktop/NCCS/Projects/For_Publication/HNSCC_CosMx6k/helpers/plotting_functions.R")

plot_results_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/Paper_Figures/Figure_1"
library('anndata')
library(dplyr)
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


###### Bulk ######
adata = ad$read_h5ad("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/1_nico_nextseq/bulk_adata.h5ad")
tab = readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/Tables.xlsx", sheet = "Metadata")
colnames(tab) = c("Patient", "Age", "Sex", "Necrosis", "Delta TIL", "BiopsyBin", "Biopsy PDL1 CPS score", "Swimmers Rank",
                                  "Survivers", "Targetted Bulk-RNASeq PC3 rank", "Notes", "Cycles", "Site", "Status")

adata$obs$DeltaTil = tab$`Delta TIL`[match(adata$obs$Patient, tab$Patient)]
adata$obs$Necrosis = tab$Necrosis[match(adata$obs$Patient, tab$Patient)]
del_cols = c(
  "1" = "#4CAF50",
  "-1" = "#F44336",
  "0" = "#CCCCCC")

p1 = plot_pca(adata,1,2,FALSE, 8, del_cols, adata$obs$DeltaTil, rep(1,length(adata$obs$DeltaTil)), c("Tit", "Delta TIL", ""))+
  theme(plot.title = element_blank(), legend.position = "none")
p2 = plot_pca(adata,3,4,FALSE, 8, del_cols, adata$obs$DeltaTil, rep(1,length(adata$obs$DeltaTil)), c("Tit", "Delta TIL", ""))+
  theme(plot.title = element_blank(), legend.title.align = 0.5)+
  guides(shape = "none", color = guide_legend(ncol = 1, override.aes=list(size =3)))
plot_save_tiff(p1+p2, "bulk_pca_delta_til", 5,10)


p1 = plot_pca(adata,1,2,FALSE, 8, patient_colors_til, adata$obs$Patient, rep(1, length(adata$obs$Patient)), c("Tit", "Patient", ""))+
  theme(plot.title = element_blank(), legend.position = "none")
p2 = plot_pca(adata,3,4,FALSE, 8, patient_colors_til, adata$obs$Patient, rep(1, length(adata$obs$Patient)), c("Tit", "Patient", ""))+
  theme(plot.title = element_blank(), legend.title.align = 0.5)+
  guides(shape = "none", color = guide_legend(ncol = 2, override.aes=list(size =3)))
plot_save_tiff(p1+p2, "bulk_pca", 5,10)


PC3 = data.frame("PC3" = adata$obsm['X_pca']$X_pca[,3],
                 "Patient" = adata$obs$Patient)
PC3 = PC3[order(PC3$PC3),]
p1 = ggplot(PC3, aes(x = reorder(Patient, PC3), y = PC3)) +
  geom_point(shape = 95, size = 10) +
  geom_hline(yintercept = 3.5, linetype = "dotted") +
  theme_Publication(14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Patient", y = "PC3 Value")


p1 = ggplot(PC3, aes(x = PC3, y = reorder(Patient, PC3))) +
  geom_point(shape = 95, size = 10) +
  geom_vline(xintercept = 3.5, linetype = "dotted") +
  scale_x_continuous(breaks = seq(floor(min(PC3$PC3)), ceiling(max(PC3$PC3)), by = 5)) +
  theme_Publication(14) +
  theme(axis.text.y = element_text(hjust = 1)) +
  labs(x = "PC3 value", y = "Patient")
plot_save_tiff(p1, "bulk_pc3_values", 6,5)


loadings_df <- data.frame(Gene = adata$var_names, Loading = adata$varm$PCs[,3])
# Select top 10 positive and negative
top_pos <- loadings_df %>% arrange(desc(Loading)) %>% slice(1:10)
top_neg <- loadings_df %>% arrange(Loading) %>% slice(1:10)
plot_df <- bind_rows(top_pos, top_neg)

p1 = ggplot(plot_df, aes(x = reorder(Gene, Loading), y = Loading,
                         fill = Loading)) +
  geom_col() +
  scale_fill_gradient2(low = "#3d5497",   # negative: dark blue
                       mid = "white",     # zero
                       high = "#9e3e3e",  # positive: dark red
                       midpoint = 0)+
  #scale_fill_manual(values = c("TRUE" = "#9e3e3e", "FALSE" = "#3d5497")) +
  coord_flip() +  # horizontal bars
  labs(x = "Gene", y = "PC3 loading",
       title = "Top Positive and Negative Genes for PC3") +
  theme_Publication(14)+
  theme(legend.position = "none", plot.title = element_blank())
plot_save_tiff(p1, "bulk_pc3_loadings", 6,4)


###### Metadata Tile ######
tab = readxl::read_excel("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/Tables.xlsx", sheet = "Metadata")
colnames(tab) = c("Patient", "Age", "Sex", "Necrosis", "Delta TIL", "BiopsyBin", "Biopsy PDL1 CPS score", "Swimmers Rank",
                  "Survivers", "Targetted Bulk-RNASeq PC3 rank", "Notes", "Cycles", "Site", "TP53 status", "TMB", "Status")

p1 = plot_tile_metadata(tab)+
  theme(legend.justification = "left")

plot_save_tiff(p1, paste0("metadata"), 7, 12)

###### Pre-Processing: Before ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P1/")
nm = "all_adata_raw_qc"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))

p1 = plot_adata_qcmetrics(adata, "Patient", patient_colors)
p2 = plot_tc_tg(adata)
p3 = plot_histogram(adata$obs$total_counts, ll = c(" ", "Total counts", "Frequency"))
right_column <- cowplot::plot_grid(p2, p3, ncol = 1, #labels = c("B", "C"),
                          label_fontfamily = "sans")
combined <- cowplot::plot_grid(p1, right_column, ncol = 2, 
                      rel_widths = c(0.9, 1),  # adjust widths if needed
                      #labels = "A",          # only label for p1
                      label_fontfamily = "sans")
plot_save_tiff(combined, "all_qc_before", 10,12)

###### Pre-Processing: After ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/tma12_1111/P1/")
nm = "all_adata_pp"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))

p1 = plot_adata_qcmetrics(adata, "Patient", patient_colors)
p2 = plot_tc_tg(adata)
p3 = plot_histogram(adata$obs$total_counts, ll = c(" ", "Total counts", "Frequency"))
right_column <- cowplot::plot_grid(p2, p3, ncol = 1,# labels = c("B", "C"),
                                   label_fontfamily = "sans")
combined <- cowplot::plot_grid(p1, right_column, ncol = 2, 
                               rel_widths = c(0.9, 1),  # adjust widths if needed
                               #labels = "A",          # only label for p1
                               label_fontfamily = "sans")
plot_save_tiff(combined, "all_qc_after", 10,12)


###### Batch Correction: After ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P1/")
nm = "all_adata_res_2_hier_annot"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

p1 = plot_umap_list(adata, c1=1, c2=2, quant = FALSE, ptsize = 0.5, pal = c(patient_colors, tt_colors),
                    color_group_list=list(adata$obs$Patient, adata$obs$TumorType), 
                    c("Patient", "Tumor type"), " ")
pp3 <- cowplot::plot_grid(p1[[1]]+
                            ggplot2::guides(colour=ggplot2::guide_legend(nrow = 4, 
                                                                         override.aes = list(shape = 15, size = 4))),
                          p1[[2]]+
                            ggplot2::guides(colour=ggplot2::guide_legend(nrow = 4, 
                                                                         override.aes = list(shape = 15, size = 4))),
                          ncol = 2, label_fontfamily = "sans")

pp1 = p1[[1]]+ ggplot2::guides(colour=ggplot2::guide_legend(nrow = 4,override.aes = list(shape = 15, size = 4)))
pp2 = p1[[2]]+ ggplot2::guides(colour=ggplot2::guide_legend(nrow = 4,override.aes = list(shape = 15, size = 4)))
plot_save_tiff(pp1, "all_umap_patient", 8,7)
plot_save_tiff(pp2+theme(legend.position = "none"), "all_umap_tt", 6,6)


p2 = plot_stacked_bar_percentage(adata$obs, "Patient", "TumorType", tt_colors, FALSE, FALSE, c("Patient", " "))+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14,margin = margin(t = -6)))+
  guides(fill = guide_legend(nrow = 2, keywidth = 0.6, keyheight = 0.6))+
  theme(legend.position = "none")
plot_save_tiff(p2, "all_stack_pat_tt", 6,8)


pp3 = ggpubr::ggarrange(pp2+ ggplot2::guides(colour=ggplot2::guide_legend(nrow = 1,override.aes = list(shape = 15, size = 4))),
                        p2, ncol = 2, common.legend = TRUE, legend = "bottom")

plot_save_tiff(pp3, "all_umap_stack_tt", 6,12)

###### Broad Annotation: UMAP, Stack and Markers ###### 
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P1/")
nm = "all_adata_res_2_hier_annot"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)

adata$obs$ct = recode(as.character(adata$obs$ct), !!!recode_ct_map)
adata$obs$ct = factor(adata$obs$ct, levels = ct_levels)

markers = c("IGHM", "HLA-DRA", "PECAM1", "KRT13", 
            "COL1A1", "CD68", "CD163", "CPA3", "TTN", "CRHR2",
            "CSF3R", "IGHG1/2", "IGKC", "CD3D", "CD3E", "KRT5")

p1 = plot_marker_umap(adata, markers[1:16], 0.2)
plot_save_tiff(p1, paste0('all_markers_m1'), 16, 16)

p1 = plot_umap(dat = adata, 1, 2, FALSE, 0.5, celltype_cols, adata$obs$ct, rep(1, length(adata$obs$ct)),
               c(paste0("PN-Tum"), "Cell type", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(nrow = 4, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "none", plot.title = element_blank())
plot_save_tiff(p1, "all_umap_ct", 6,7)

p2 = plot_tsne(dat = adata, 1, 2, FALSE, 0.5, celltype_cols, adata$obs$ct, rep(1, length(adata$obs$ct)),
               c(paste0("PN-Tum"), "Cell type", " "))+
  ggplot2::guides(shape = "none")+
  ggplot2::guides(colour=ggplot2::guide_legend(nrow = 4, override.aes = list(shape = 15, size = 4)))+
  ggplot2::theme(legend.position = "bottom", plot.title = element_blank())
plot_save_tiff(p2, "all_tsne_ct", 9,8)

p2 = plot_stacked_bar_percentage(adata$obs, "Patient", "ct", c(celltype_cols), FALSE, FALSE, c("Patient", " "))+
  labs(fill = "Cell type")+
  ggplot2::theme(legend.position = "right", legend.title.align = 0.5) + 
  theme(plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14,margin = margin(t = -10)))+
  guides(fill = guide_legend(ncol = 1, keywidth = 1,keyheight = 1))
plot_save_tiff(p2, "all_stack_pat_ct", 6,12)


###### Dotplot for the markers used for the annotation ###### 
library(ggplot2)
library(viridis)
library(dplyr)
library(scales)
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P1/")
expr_summary <- read.csv(paste0(data_dir, "alldeg_dotplot_summary.csv"))
expr_summary$cluster <- factor(expr_summary$cluster, levels = unique(expr_summary$cluster))
expr_summary$gene <- factor(expr_summary$gene, levels = unique(expr_summary$gene))
expr_summary <- expr_summary %>%
  group_by(gene) %>%
  mutate(logFC_scaled = rescale(logFC, to = c(-2, 2))) %>%
  ungroup()

expr_summary$cluster = recode(as.character(expr_summary$cluster), !!!recode_ct_map)
expr_summary$cluster = factor(expr_summary$cluster, levels = ct_levels)
expr_summary$source_cluster = recode(as.character(expr_summary$source_cluster), !!!recode_ct_map)
expr_summary$source_cluster = factor(expr_summary$source_cluster, levels = ct_levels)

p1 = plot_dotplot(expr_summary, "gene", "cluster", "logFC_scaled", "pct_expr", c("logFC", "% expressed"))
plot_save_tiff(p1, "all_dotplot_ct", 6,15)

###### Response based broad alluvial plots for Resp vs non Resp ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P1/")
nm = "all_adata_res_2_hier_annot_resp"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)
adata$obs$ct = recode(as.character(adata$obs$ct), !!!recode_ct_map)
adata$obs$ct = factor(adata$obs$ct, levels = ct_levels)
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)

resp_pp  = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
resp_pp = resp_pp%>% filter(Response %in% c("Response"))
p1 = plot_alluvium(resp_pp, "TumorType", "ct", celltype_cols, c("Response", "Cell type proportion (%)", "Cell type"))
p1 = p1 + coord_cartesian(clip = "off") + annotate("text",x = 0.55,y = 0.5,label = "Cell type proportion (%)",
    angle = 90,hjust = 0.5,vjust = 0.5,size = 4)
plot_save_tiff(p1, "all_alluvial_resp", 6,7)

nresp_pp  = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
nresp_pp = nresp_pp%>% filter(Response %in% c("Non-response"))
p2 = plot_alluvium(nresp_pp, "TumorType", "ct", celltype_cols, c("Non-response", "Cell type proportion (%)", "Cell type"))
p2 = p2 + coord_cartesian(clip = "off") + annotate("text",x = 0.55,y = 0.5,label = "Cell type proportion (%)",
                                                   angle = 90,hjust = 0.5,vjust = 0.5,size = 4)
plot_save_tiff(p2, "all_alluvial_non_resp", 6,7)

###### Status Boxplots with pval ######
data_dir = paste0("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/results/final/P1/")
nm = "all_adata_res_2_hier_annot_resp"
adata = ad$read_h5ad(paste0(data_dir, nm, ".h5ad"))
adata$obs$Patient = factor(adata$obs$Patient, levels = pat_levels2)
adata$obs$TumorType = recode(as.character(adata$obs$TumorType), !!!recode_tt_map)
adata$obs$TumorType = factor(adata$obs$TumorType, levels = tt_levels)
adata$obs$ct = recode(as.character(adata$obs$ct), !!!recode_ct_map)
adata$obs$ct = factor(adata$obs$ct, levels = ct_levels)
adata$obs$Response = recode(as.character(adata$obs$Response), !!!recode_resp_map)
adata$obs$Response = factor(adata$obs$Response, levels = resp_levels)


df = adata$obs %>% filter(TumorType %in% c("Pre-NIVO primary", "Post-NIVO primary"))
for (celltype_of_interest in unique(adata$obs$ct)) {
  p1 = plot_g4_paired_boxplot(df, celltype_of_interest, 'ct')
  p1 = p1 + ggplot2::theme(legend.position = "none")
  plot_save_tiff(p1, paste0("all_boxplot_resp_", celltype_of_interest), 6,6)
}

p1 = plot_g4_paired_boxplot(df, "Plasma cells", 'ct')
p1 = p1 + ggplot2::theme(legend.position = "none")
plot_save_tiff(p1, paste0("all_boxplot_resp_", "Plasma cells"), 6,4)

p1 = plot_g4_paired_boxplot(df, "Tumor", 'ct')
p1 = p1 + ggplot2::theme(legend.position = "none")
plot_save_tiff(p1, paste0("all_boxplot_resp_", "Tumor"), 6,4)

p1 = plot_g4_paired_boxplot(df, "Tumor", 'ct')
p1 = p1 + ggplot2::theme(legend.position = "right")
plot_save_tiff(p1, paste0("all_boxplot_resp_", "Legend_Ref"), 4,4)
