#library(randomcoloR)

library(RColorBrewer)

tt_colors <- c(
  "Pre-NIVO primary" = '#4363d8',  # Red
  "Post-NIVO primary" = '#3cb44b',  # Blue
  "Post-NIVO normal" = '#ffe119',  # Green
  "Post-NIVO adv front" = "#E91E63",   # Purple
  "Post-NIVO uninvolved node distant (all)" = '#3cb44b',  # Blue
  "Post-NIVO uninvolved node local (all)" = '#ffe119',  # Green
  "Post-NIVO involved node (+ve only)" = "#E91E63"   # Purple
)

response_colors <- c(
  "Response" = "#4CAF50",
  "Non-response" = "#F44336")

g4_colors <- c(
  "Pre-R" = "#65b968",
  "Post-R" = "#054407",
  "Pre-NR" = "#f3776e",
  "Post-NR" = "#830c03",
  "Post-Adv-R" = "#67bf24",
  "Post-Adv-NR" = "#bf4624")

distance_colors_orange <- c(
  "Q0: TN"      = "#FDD0A2",
  "Q1: 0-250"   = "#FDAE6B",
  "Q2: 250-500" = "#F16913",
  "Q3: 500-800" = "#D94801",
  "Q4: 800+"    = "#8C2D04"
)
patient_colors <- c(
  "2-011" = "#1f77b4",
  "2-020" = "#ff7f0e",
  "2-017" = "#2ca02c",
  "2-028" = "#d62728",
  "2-013" = "#9467bd",
  "2-026" = "#8c564b",
  "2-022" = "#e377c2",
  "2-034" = "#7f7f7f",
  "2-039" = "#bcbd22",
  "2-033" = "#17becf",
  "2-018" = "#aec7e8",
  "2-032" = "#ffbb78",
  "2-036" = "#98df8a",
  "2-040" = "#ff9896",
  "20-010" = "#c5b0d5",
  "20-002" = "#c49c94",
  "2-041" = "#f7b6d2",
  "2-012" = "#c7c7c7",
  "2-029" = "#dbdb8d"
)

celltype_cols = c('Tumor' = "#FF007F",
                  'Epithelial cells' = '#f032e6','T cells' =  "#2f53ed",
                  "B cells" = "#80b736", "NK cells" ='#000075',
                  'Macrophages' = '#ffe119', 'Fibroblasts' = '#f58231',
                  'Endothelial cells' = "#a60401", "Plasma cells" = '#175616',
                  "Dendritic cells" = "#4bc3b7", "Myogenic cells" =  '#911eb4',
                  "Mast cells" = "#cf4f33", "Unknown" = "#4B4B4B", "Neuronal cells" = "#d79a36", "Neutrophils" =  "#4d2a22",
                  "Adipocytes" = "#34495E","Myeloid cells" = '#911eb4',
                  'CD8+ T cells' =  "#2f53ed","CD4+ T cells" = '#7b82e0')

tsub_l1 = c( "CD4+ CD8+" = "#3cb44b",
             "CD4+ CD8-" = "#4f95cc",
             "CD4- CD8+" = "#F77F00",
             "CD4- CD8-" = "#dcbeff")


cd8_l1 <- c(
  "CD8_T_CYT_DYS" = "#c76b0a",    # Red-Orange
  "CD8_T_EM" = "#f2cb0a",        # Blue
  "CD8_T_SL" = "#8E44AD",       # Purple
  "CD8_T_DYS_2" = "#1F4BFF",        # Yellow
  "CD8_T_DYS_1" = "#00BFFF",
  "CD8_T_CM" = "#e4859d",        # Light Green
  "CD8_T_Q" = "#5cd90f",        # Teal
  "CD8_T_N" = "#f032e6",     # Light Blue
  "CD8_T_TR" = "#a60401"        # Dark Red
)

cd4_l1 <- c(
  "CD4_T_TH17" = "#c76b0a",    # Red-Orange
  "CD4_T_EM" = "#f2cb0a",        # Blue
  "CD4_T_SL" = "#8E44AD",       # Purple
  "CD4_T_DYS" = "#1F4BFF",  
  "CD4_T_TFH" = "#00BFFF",# Yello
  "CD4_T_CM" = "#e4859d",        # Light Green
  "CD4_T_REG" = "#5cd90f",        # Teal
  "CD4_T_N" = "#f032e6",     # Light Blue
  "CD4_T_TR" = "#a60401",        # Dark Red,
  "CD4_T_TH9" = "#dcbeff"
)

niche_colors_con <- c(
  "Tumor-rich niche"             = "#E41A1C",  # red
  "Neuronal-rich niche"    = "#FF7F00",  # purple
  "T-rich niche"           = "#377EB8",  # blue
  "Epithelial-rich niche"  = "#f032e6",  # orange
  "Fibroblasts-rich niche"       = "#4DAF4A",  # green
  "Mast-rich niche"        = "#A65628",  # brown
  "Macrophages-rich niche"       = "#FFD92F",  # pink
  "Dendritic-rich niche"   = "#66C2A5",  # grey
  "Endothelial-rich niche" = "#00BFFF",  # teal
  "B-rich niche"           = "#E78AC3",  # yellow
  "Myogenic-rich niche"    = "#c4794d",  # light pink
  "Plasma-rich niche"      = "#984EA3",  # aqua
  "Neutrophils-rich niche"       = "#4d2a22",   # lavender
  "No-contact niche" = "#34495E"
)

niche_colors_spat <- c(
  "Tumor-rich niche"        = "#E41A1C",  # strong red
  "Fibroblasts-rich niche"  = "#4DAF4A",  # green
  "Plasma-rich niche"       = "#984EA3",  # purple
  "T-rich niche"            = "#377EB8",  # blue
  "Neuronal-rich niche"     = "#FF7F00",   # orange
  "Epithelial-rich niche" = "#f032e6",
  "Neutrophils-rich niche" = "#4d2a22",
  "Macrophages-rich niche" = "#FFD92F"
)

spfp_colors <- c(
  "Epit"                 = "#E41A1C",  # strong red
  "Tum"                  = "#377EB8",  # blue
  "Neutro"               = "#4DAF4A",  # green
  "Tum_Fib_Mac"          = "#984EA3",  # purple
  "Mac_Den_Frib_Mast_NK" = "#FF7F00",  # orange
  "Myo_Endo"             = "#A65628",  # brown
  "Plasma_Fib"           = "#F781BF",  # pink
  "Neuron"               = "#999999",  # grey
  "T_NK_B"               = "#66C2A5",  # teal
  "Tum_T"                = "#FFD92F"   # yellow
)

c6 = c("#2f53ed","#537d32","#db9b48","#d94a1f","#7e0077","#ff005e")

col_rb = c("#a60401", "#00418f")

pre_post_col = c("#a74959", "#4bc3b7")
c20 = c("#FF007F", "#1D3557", "#3B8A57", "#F77F00", 
        "#6A0572", "#00BFFF", "#FFCC00", "#008C8C",
        "#4B4B4B", "#E63946", "#c4794d","#e4859d",
        "#9a56ca", "#a60401", "#4bc3b7", "#826e2d",
        "#2f53ed", "#80b736", "#d77dd5", "#cf4f33",
        "#4d2a22", "#00418f", "#a293da","#c4794d")

c10 = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
        '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', 
        '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', 
        '#000075', '#a9a9a9', '#ffffff', '#000000',
        "#FF007F", "#1D3557", "#3B8A57", "#F77F00", 
        "#6A0572", "#00BFFF", "#FFCC00", "#008C8C",
        "#4B4B4B", "#E63946", "#c4794d","#e4859d",
        "#9a56ca", "#a60401", "#4bc3b7", "#826e2d",
        "#2f53ed", "#80b736", "#d77dd5", "#cf4f33",
        "#4d2a22", "#00418f", "#a293da","#c4794d")


col6 = c("#a60401", "#00418f", "#4d2a22", "#0d8d01", "#0981b4", "#ad6105")

col_wr = c("#D3D3D3","#00418f")
col_wr = c("#808080", "#F0F921", "#a60401")

col_wr = c("#D3D3D3","#0981b4","#2f53ed")

cb = c('#000000')

celltype_cols1 = c('Tumor' = "#FF007F", 'Tumour' = "#FF007F",
                  'Epithelial Cells' = '#f032e6','T Cells' =  "#2f53ed",
                  "B Cells" = "#80b736", "NK Cells" ='#000075',
                  'Macrophages' = '#ffe119', 'Fibroblasts' = '#f58231',
                  'Endothelial Cells' = "#a60401", "Plasma Cells" = '#175616',
                  "Dentritic Cells" = "#4bc3b7", "Myogenic Cells" =  '#911eb4',
                  "Mast Cells" = "#cf4f33", "Unknown" = "#4B4B4B", "Neuronal Cells" = "#d79a36", "Neutrophils" =  "#4d2a22",
                  "Adipocytes" = "#34495E")

theme_Publication  =  function(base_size=14, base_family="sans") {
  ggthemes::theme_foundation(base_size=base_size, base_family=base_family) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      size = ggplot2::rel(1.2), hjust = 0.5),
                   text = ggplot2::element_text(),
                   panel.background = ggplot2::element_rect(fill = "white",colour = NA),
                   plot.background = ggplot2::element_rect(fill = "white",colour = NA),
                   panel.border = ggplot2::element_rect(colour = NA),
                   axis.title = ggplot2::element_text(size = ggplot2::rel(1)),
                   axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                   axis.title.x = ggplot2::element_text(vjust = -0.2),
                   axis.text = element_text(size = rel(1), color = "black"),
                   axis.line = ggplot2::element_line(colour="black"),
                   axis.ticks = ggplot2::element_line(),
                   panel.grid.major = ggplot2::element_blank(), #element_line(colour="#f0f0f0"),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.key = ggplot2::element_rect(colour = NA),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(1), color = "black"),      
                   legend.position = "right",
                   legend.direction = NULL,         
                   legend.justification = "center",
                   legend.key.size= ggplot2::unit(1, "cm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
                   legend.title = ggplot2::element_text(),
                   plot.title.position = "plot" ,
                   plot.margin = unit(c(12, 6, 8, 6), "mm"),
                   strip.background=ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                   strip.text = ggplot2::element_text(face="bold"))}


theme_th1  =  function(base_size=18, base_family="Times") {
  ggthemes::theme_foundation(base_size=base_size, base_family=base_family) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      size = ggplot2::rel(1.2), hjust = 0.5),
                   text = ggplot2::element_text(),
                   panel.background = element_rect(fill = "white", colour = NA),
                   axis.title = ggplot2::element_text(size = ggplot2::rel(1)),
                   axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                   axis.title.x = ggplot2::element_text(vjust = -0.2),
                   axis.text = element_text(size = rel(1.1), color = "black"),
                   axis.line = ggplot2::element_line(colour="black"),
                   axis.ticks = ggplot2::element_line(),
                   panel.grid.major = ggplot2::element_blank(), #element_line(colour="#f0f0f0"),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.key = ggplot2::element_rect(colour = NA),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(1), color = "black"),      
                   legend.position = "right",
                   legend.direction = NULL,         
                   legend.justification = "center",
                   legend.key.size= ggplot2::unit(1, "cm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
                   legend.title = ggplot2::element_text(),
                   plot.title.position = "plot" ,
                   plot.margin = unit(c(12, 6, 8, 6), "mm"),
                   strip.background=ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                   strip.text = ggplot2::element_text(face="bold")
    )
}

plot_save = function(p, name, bh = 3, bw = 4){
  fname = paste0(plot_results_dir, "/", name, ".png")
  cowplot::save_plot(fname, p+  theme(
    panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white")),  base_height = bh, base_width = bw)
} 


plot_save_tiff = function(p, name, bh = 3, bw = 4, dpi = 300) {
  fname = file.path(plot_results_dir, paste0(name, ".tiff"))
  cowplot::save_plot(
    filename = fname,
    plot = p,
    base_height = bh,
    base_width = bw,
    dpi = dpi,
    device = "tiff"
  )
  message("Saved: ", fname)
}