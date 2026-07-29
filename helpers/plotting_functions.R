###### Scatters ######
plot_scatter = function(df,quant=FALSE, ptsize = 0.5, pal = palette("default"),
                        ll = c("Simple Scatter Plot","Comp1","Comp2",
                               "Color Factor", "Shape Factor")){
  colnames(df) = c("Comp1", "Comp2", "Col_gr", "Sh_gr")
  if (quant) {
    color_group = df[, 3]
  }
  else {
    color_group = factor(df[, 3])
  }
  shape_group = factor(df[, 4])
  sh_chars = c(20,15, seq(1:25), 0,33,35,36,37,38,42,43,63,64)
  p = ggplot2::ggplot(df, ggplot2::aes(x = Comp1, y = Comp2)) +
    ggplot2::geom_point(ggplot2::aes(shape = shape_group, color = color_group),
                        size = ptsize, alpha = 0.9, stroke = 0.2, position = "jitter") +
    ggplot2::scale_shape_manual(values = sh_chars[1:length(unique(shape_group))])
  if (quant) {
    p = ggplot2::ggplot(df, ggplot2::aes(x = Comp1, y = Comp2)) +
      ggplot2::geom_point(ggplot2::aes(shape = shape_group, color = color_group),
                          size = ptsize, alpha = 0.7, stroke = 0.2, position = "jitter") +
      ggplot2::scale_shape_manual(values = sh_chars[1:length(unique(shape_group))])
    p1 = p + ggplot2::scale_color_gradientn(colors = pal)
  }
  else {
    p1 = p + ggplot2::scale_color_manual(values = pal)
  }
  p2 = p1 + theme_Publication(base_size = 14) +
    ggplot2::labs(title = ll[1],x = ll[2], y = ll[3],color = ll[4], shape = ll[5])
  return(p2)
}

plot_tc_tg = function(adata){
  #tg = rowSums(as.matrix(adata$X) > 0)
  #tc =  rowSums(as.matrix(adata$X))
  tg <- Matrix::rowSums(adata$X > 0)
  tc <- Matrix::rowSums(adata$X)
  df = data.frame(Comp1 = tc, Comp2 = tg,
                  color_group = 1, shape_group = 0)
  p =plot_scatter(df, FALSE, 0.5, cb, c("cvsg", "Total counts per cell", "Total genes per cell", "Patients", "Sh"))+
    ggplot2::guides(shape = "none")+
    ggplot2::guides(colour="none")+
    theme_Publication()+
    ggplot2::theme(plot.title = ggplot2::element_blank(), legend.position = "none")
  return(p)
}

plot_pca = function(dat, c1=1, c2=2, quant=FALSE, ptsize = 0.5, pal = palette("default"),
                    color_group, shape_group,
                    ll = c("PCA Plot", "Colour Factor", "Shape Factor")) {
  var_pc = dat$uns['pca']$pca$variance_ratio
  pc_labels = sprintf("PC%d (%.2f%%)", 1:10, var_pc * 100)
  df = data.frame("PC1" = dat$obsm$X_pca[, c1],
                  "PC2" = dat$obsm$X_pca[, c2],
                  "Col_gr" = color_group,
                  "Sh_gr" = shape_group)
  p = plot_scatter(df, quant, ptsize, pal,
                   c(ll[1], pc_labels[c1], pc_labels[c2], ll[2], ll[3])) +
    ggplot2::theme(axis.line = ggplot2::element_blank())#+
  return(p)
}

plot_umap = function(dat, c1=1, c2=2,quant=FALSE, ptsize = 0.5, pal = palette("default"),
                     color_group, shape_group,
                     ll = c("PCA Plot", "Colour Factor", "Shape Factor")) {
  umap_labels = sprintf("UMAP%d", 1:10)
  df = data.frame("Comp1" = dat$obsm$X_umap[, c1],
                  "Comp2" = dat$obsm$X_umap[, c2],
                  "Col_gr" = color_group,
                  "Sh_gr" = shape_group)
  p = plot_scatter(df, quant,ptsize, pal,
                   c(ll[1], umap_labels[c1], umap_labels[c2], ll[2], ll[3])) +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank())
  return(p)
}

plot_tsne = function(dat, c1=1, c2=2,quant=FALSE, ptsize = 0.5, pal = palette("default"),
                     color_group, shape_group,
                     ll = c("PCA Plot", "Colour Factor", "Shape Factor")) {
  umap_labels = sprintf("TSNE%d", 1:10)
  df = data.frame("Comp1" = dat$obsm$X_tsne[, c1],
                  "Comp2" = dat$obsm$X_tsne[, c2],
                  "Col_gr" = color_group,
                  "Sh_gr" = shape_group)
  p = plot_scatter(df, quant,ptsize, pal,
                   c(ll[1], umap_labels[c1], umap_labels[c2], ll[2], ll[3])) +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank())
  return(p)
}


plot_pca_list = function(dat, c1=1, c2=2, quant = FALSE, ptsize = 0.5, pal = palette("default"),
                         color_group_list, color_group_labels){
  plot_list = list()
  for (i in 1:length(color_group_list)){
    plot_list[[i]] = plot_pca(dat, c1, c2, quant, ptsize, c20, color_group_list[[i]], rep(1, length(dat$obs$fov)),
                              c("PCA Plot", color_group_labels[[i]], " "))+
      ggplot2::guides(shape = "none")+
      ggplot2::guides(colour=ggplot2::guide_legend(nrow = 4, override.aes = list(shape = 15,size = 4)))+
      ggplot2::theme(plot.title = ggplot2::element_blank(), legend.position = "bottom")
  }
  return(plot_list)
} 

plot_umap_list = function(dat, c1=1, c2=2, quant = FALSE, ptsize = 0.5, pal = palette("default"),
                          color_group_list, color_group_labels, tit){
  plot_list = list()
  for (i in 1:length(color_group_list)){
    plot_list[[i]] = plot_umap(dat, c1, c2, quant, ptsize, pal, color_group_list[[i]], rep(1, length(color_group_list[[i]])),
                               c(tit, color_group_labels[[i]], " "))+
      ggplot2::guides(shape = "none")+
      ggplot2::guides(colour=ggplot2::guide_legend(nrow = 3, override.aes = list(shape = 15, size = 4)))+
      ggplot2::theme(legend.position = "bottom")
  }
  return(plot_list)
} 

plot_tsne_list = function(dat, c1=1, c2=2, quant = FALSE, ptsize = 0.5, pal = palette("default"),
                          color_group_list, color_group_labels, tit){
  plot_list = list()
  for (i in 1:length(color_group_list)){
    plot_list[[i]] = plot_tsne(dat, c1, c2, quant, ptsize, pal, color_group_list[[i]], rep(1, length(color_group_list[[i]])),
                               c(tit, color_group_labels[[i]], " "))+
      ggplot2::guides(shape = "none")+
      ggplot2::guides(colour=ggplot2::guide_legend(nrow = 4, override.aes = list(shape = 15,size = 4)))+
      ggplot2::theme(legend.position = "bottom")
  }
  return(plot_list)
} 

plot_marker_umap = function(dat, markers, sz){
  plot_list = list()
  for (i in seq(1,length(markers))){
    print(i)
    print(markers[i])
    #scaled_ex = scale(dat$X[,which(dat$var_names == markers[i])])[,1]
    #scaled_ex = min_max_normalize(dat$X[,which(dat$var_names == markers[i])])
    #scaled_ex = z_norm(dat$X[,which(dat$var_names == markers[i])])
    #colorspace::sequential_hcl(5, palette = "Heat 2"),
    scaled_ex = dat$X[,which(dat$var_names == markers[i])]
    plot_list[[i]] = plot_umap(dat,1,2,TRUE, sz, viridis::viridis(5),  scaled_ex , 1, 
                               c(markers[i], " ", ""))+
      ggplot2::theme(legend.position = "bottom",
                     axis.title = element_blank())+
      ggplot2::guides(shape = "none")+
      #ggplot2::guides(colour=ggplot2::guide_legend(nrow = 2, override.aes = list(size = 4)))+
      ggplot2::theme(legend.position = "right")
  }
  pf = ggpubr::ggarrange(plotlist = plot_list, common.legend = TRUE, legend = "right")+
    theme_Publication()
  return(pf)
}

plot_bubble = function(dat, xx, yy, size_factor, color_factor, pal, cc = c("Bubble Plot", "Xlab", "Ylab", "sizeFac", "colFac")){
  p1 =  ggplot(dat, aes(x = !!sym(xx), y = !!sym(yy))) +
    geom_point(aes(size = !!sym(size_factor), color = !!sym(color_factor)), alpha = 0.8) +
    scale_size_continuous(range = c(2, 10)) +  
    scale_color_manual(values = pal) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
    theme_Publication(base_size = 14) +
    theme(axis.text.y = element_text(size = 10)) +
    labs(x = cc[[2]],
         y =cc[[3]],
         size = cc[[4]],
         color = cc[[5]],
         title = cc[[1]])
  
  return(p1)
}
###### Metadata Plot Tile ######
plot_tile_metadata <- function(df) {
  df <- df %>%
    mutate(across(all_of(c("Sex", "Necrosis", "Delta TIL", "Site", "Status")), as.factor))
  df$Patient <- factor(df$Patient, levels = df$Patient[order(df$`Targetted Bulk-RNASeq PC3 rank`)])
  df <- df %>%
    mutate(`Targetted Bulk-RNASeq PC3 rank` =
             na_if(`Targetted Bulk-RNASeq PC3 rank`, 6.5)) 
  
  df$Status = factor(df$Status, levels = resp_levels)
  p1 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "Status", fill = Status), color = "white") +
    scale_fill_manual(values = response_colors, name = "Status")+
    theme_Publication()+
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_legend(keywidth = 0.6,keyheight = 0.6, nrow = 1))
  
  p2 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "Age", fill = Age), color = "white") +
    scale_fill_gradient(
      low  = "#CCCCCC",
      high = "#70410a"
    )+
    
    #scale_fill_viridis(option = "C", direction = 1, na.value = "grey90", name = "Age") +
    theme_Publication()+
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_colorbar(direction = "horizontal", 
                                 ticks = FALSE,
                                 label = TRUE, #keep it False, make it true to see and annotate manually in ppt
                                 barwidth = 10,
                                 barheight = 0.35))
  
  df$Sex = factor(df$Sex, levels = c("M", "F"))
  sex_cols = c("M" = "#2e5b9b" ,
               "F" = "#ffaaf0")
  p22 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "Sex", fill = Sex), color = "white") +
    scale_fill_manual(values = sex_cols, name = "Sex")+
    theme_Publication()+
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_legend(keywidth = 0.6,keyheight = 0.6, nrow = 1))
  
  
  
  
  
  df$Necrosis = factor(df$Necrosis, levels = c("Yes", "No"))
  nec_cols = c("Yes" = "#2e5b9b" ,
               "No" = "#CCCCCC")
  p3 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "Necrosis", fill = Necrosis), color = "white") +
    scale_fill_manual(values = nec_cols, name = "Necrosis")+
    theme_Publication()+
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_legend(keywidth = 0.6,keyheight = 0.6, nrow = 1))
  
  df$`Delta TIL` = factor(df$`Delta TIL`, levels = c(1,0,-1))
  del_cols = c(
    "1" = "#4CAF50",
    "-1" = "#F44336",
    "0" = "#CCCCCC")
  p4 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "Delta TIL", fill = `Delta TIL`), color = "white") +
    scale_fill_manual(values = del_cols, name = "Delta TIL")+
    theme_Publication()+
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_legend(keywidth = 0.6,keyheight = 0.6, nrow = 1))
  
  
  cps_cols = c("<1%" = "#CCCCCC",
               "2-3%" = "#a8b6cb",
               "5-10%" = "#6e84a4",
               "10-20%" =  "#9ECAE1",
               "20-30%" = "#4A90C2",
               "60-70%" = "#2e5b9b")
  df$BiopsyBin = factor(df$BiopsyBin, levels = c("<1%","2-3%","5-10%","10-20%","20-30%","60-70%"))
  p5 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "Biopsy PDL1 \n CPS score", fill = BiopsyBin), color = "white") +
    scale_fill_manual(values = cps_cols, name = "Biopsy PDL1\nCPS score")+
    #scale_fill_viridis(option = "C", direction = 1, na.value = "grey90", name = "Biopsy PDL1 \n CPS Score") +
    theme_Publication()+
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_legend(keywidth = 0.6,keyheight = 0.6, nrow = 2))
  #guides(fill = guide_colorbar(direction = "horizontal",
  #                             ticks = FALSE,
  #                             label = FALSE,
  #                             barwidth = 10,
  #                             barheight = 0.35))
  
  p6 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "Targetted\nBulk-RNASeq\nPC3 rank", fill = `Targetted Bulk-RNASeq PC3 rank`), color = "white") +
    scale_fill_viridis(option = "C", direction = 1, na.value = "white", name = "Targetted Bulk-RNASeq \n PC3 rank") +
    theme_Publication()+
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_colorbar(direction = "horizontal", 
                                 ticks = FALSE,
                                 label = TRUE,
                                 barwidth = 10,
                                 barheight = 0.35))
  
  
  p7 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "Site", fill = Site), color = "white") +
    scale_fill_manual(values = c20, name = "Site")+
    theme_Publication()+
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_legend(keywidth = 0.6,keyheight = 0.6, nrow = 3))
  
  tight <- theme(plot.margin = margin(0, 0, 0, 0))+
    theme(legend.justification = "left")
  
  p8 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "Swimmers Rank", fill = `Swimmers Rank`), color = "white") +
    scale_fill_viridis(option = "C", direction = 1, na.value = "white", name = "Swimmers\nRank") +
    theme_Publication()+
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_colorbar(direction = "horizontal", 
                                 ticks = FALSE,
                                 label = FALSE,
                                 barwidth = 10,
                                 barheight = 0.35))
  
  df$`TP53 status` = factor(df$`TP53 status`, levels = c("WT", "Mutant"))
  TP53_cols = c(
    "Mutant" = "#ffaaf0",
    "WT" = "#2e5b9b")
  p9 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "TP53 status", fill = `TP53 status`), color = "white") +
    scale_fill_manual(values = TP53_cols, name = "TP53 status",na.value = "white",
                      breaks = c("WT", "Mutant"),   # 🔥 key line
                      drop = TRUE)+
    theme_Publication()+
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_legend(keywidth = 0.6,keyheight = 0.6, nrow = 1))
  
  p10 <- ggplot() +
    geom_tile(data = df, aes(x = Patient, y = "TMB", fill = `TMB`), color = "white") +
    scale_fill_viridis(option = "C", direction = 1, na.value = "white", name = "TMB") +
    theme_Publication()+
    theme(axis.text.x = element_blank(), # show names
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(), 
          axis.title.x = element_blank(),  # remove x-axis title
          axis.title.y = element_blank(),  # remove y-axis title
          panel.grid = element_blank(),
          legend.title = element_blank())+
    guides(fill = guide_colorbar(direction = "horizontal", 
                                 ticks = FALSE,
                                 label = TRUE,
                                 barwidth = 10,
                                 barheight = 0.35))
  
  
  
  
  p_all <- (p1 + tight) /
    (p2 + tight) /
    (p22 + tight) /
    (p3 + tight) /
    (p4 + tight) /
    (p5 + tight) /
    (p6 + tight) /
    (p7 + tight) / #/ (p8 + tight)
    (p9 + tight) /
    (p10 + tight)
  
  return(p_all)
}


###### Violins ######
plot_violin = function(data, value_col, condition_col, pal = palette("default"), violin = TRUE,
                       ll = c("Violin Plot", "Condition", "Frequency", "Condition")) {
  if (violin == TRUE){
    p = ggplot2::ggplot(data, ggplot2::aes(x = !!ggplot2::sym(condition_col), 
                                           y = !!ggplot2::sym(value_col), 
                                           fill = !!ggplot2::sym(condition_col))) +
      ggplot2::geom_violin(draw_quantiles = c(0.5), trim = FALSE)
  }
  else{
    p = ggplot2::ggplot(data, ggplot2::aes(x = !!ggplot2::sym(condition_col), 
                                           y = !!ggplot2::sym(value_col), 
                                           fill = !!ggplot2::sym(condition_col))) +
      #ggplot2::geom_violin(draw_quantiles = c(0.5)) +
      ggplot2::geom_boxplot(outlier.shape = NA)
  }
  pf = p + ggplot2::scale_fill_manual(values = pal, name = condition_col) +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 3))+
    ggplot2::labs(title = ll[1], x = ll[2], y = ll[3])+
    theme_Publication()
  return(pf)
}

plot_adata_qcmetrics = function(adata, categ, c6){
  print(adata$obs$n_genes_by_counts[1:10])
  p1 = plot_violin(adata$obs, "n_genes_by_counts", categ, c6, TRUE,
                   c("Total genes per cell", categ, "Total genes", "Patients"))+
    ggplot2::theme(legend.position = "None", plot.title = ggplot2::element_blank())
  
  p2 = plot_violin(adata$obs, "total_counts", categ, c6, TRUE,
                   c("Total counts per cell", categ, "Total counts", "Slide"))+
    ggplot2::theme(legend.position = "None", plot.title = ggplot2::element_blank())
  
  p3 = plot_violin(adata$obs, "pct_counts_Negative", categ, c6, TRUE,
                   c("% of Negative Probes per cell", categ, "% of negative probes", "Slide"))+
    ggplot2::theme(legend.position = "None", plot.title = ggplot2::element_blank())
  
  p4 = plot_violin(adata$obs, "pct_counts_SystemControl", categ, c6, TRUE,
                   c("% of SysControl Probes per cell", categ, "% of control probes", "Slide"))+
    ggplot2::theme(legend.position = "None", plot.title = ggplot2::element_blank())
  
  p11 = p1+ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                          axis.ticks.x = ggplot2::element_blank(),
                          axis.title.x = ggplot2::element_blank(),
                          plot.margin = ggplot2::unit(c(0.5,0,0,0.5), "cm"))
  p22 = p2+ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                          axis.ticks.x = ggplot2::element_blank(),
                          axis.title.x = ggplot2::element_blank(),
                          plot.margin = ggplot2::unit(c(0.5,0,0,0.5), "cm"))
  p33 = p3+ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                          axis.ticks.x = ggplot2::element_blank(),
                          axis.title.x = ggplot2::element_blank(),
                          plot.margin = ggplot2::unit(c(0.5,0,0,0.5), "cm"))
  p44 = p4+ggplot2::theme(plot.margin = ggplot2::unit(c(0.8,0,0.5,0.5), "cm"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p1234 = p11/p22/p33/p44
  return (p1234)
}



###### Histograms ######
plot_histogram <- function(values, bins = 30, pal = palette("default"),
                           ll = c("Histogram of Column", "Column", "Frequency")) {
  data <- data.frame(values)
  colnames(data) = c("counts")
  
  ggplot2::ggplot(data, ggplot2::aes(x = counts)) +
    ggplot2::geom_histogram(bins = bins, fill = "#c4794d", color = 'black', position = "identity") +
    ggplot2::labs(title = ll[[1]],
                  x = ll[[2]],
                  y = ll[[3]]) +
    theme_Publication()+
    ggplot2::theme(plot.title = ggplot2::element_blank())
}

plot_overlayed_histogram <- function(df, x_col, cat_col, bins = 30, pal, ll = c("Overlayed Histograms", "Distance", "Frequency")) {
  # Convert column names to symbols for tidy evaluation
  x_sym <- rlang::sym(x_col)
  cat_sym <- rlang::sym(cat_col)
  
  ggplot(df, aes(x = !!x_sym, fill = factor(!!cat_sym))) +
    #geom_histogram(position = "identity",bins = bins,alpha = 0.8,color = "black") +
    geom_histogram(position = position_identity(), bins = bins, alpha = 0.4, color = "black", linewidth = 0.2)+
    scale_fill_manual(values = pal) +
    labs(title = ll[1], x = ll[2], y = ll[3]) +
    theme_Publication(base_size = 12)+
    theme(legend.position = "bottom")+
    guides(fill = guide_legend(keywidth = 0.6,keyheight = 0.6))}


plot_hist_lines_points <- function(df, x_col, cat_col, bins = 30, pal,
                                   ll = c("Histogram Lines", "Distance", "Cell count")) {
 
  x_sym <- rlang::sym(x_col)
  cat_sym <- rlang::sym(cat_col)
  # Bin + count
  df_binned <- df %>%
    mutate(bin = cut(!!x_sym, breaks = bins)) %>%
    group_by(bin, !!cat_sym) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      bin_center = (as.numeric(sub("\\((.+),.*", "\\1", bin)) +
                      as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", bin))) / 2
    )
  
  # Plot
  p1 = ggplot(df_binned, aes(x = bin_center, y = count, color = factor(!!cat_sym))) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 0.8) +
    scale_color_manual(values = pal) +
    theme_Publication(base_size = 12)+
    labs(title = ll[1], x = ll[2], y = ll[3]) +
    theme(legend.position = "bottom")
    
  return(p1)
}


###### Barplots ######
plot_stacked_bar_percentage <- function(data, x_obs, y_obs, pal, labell = FALSE, horizontal = TRUE, cc = "") {
  print(cc)
  summary_data <- data %>%
    dplyr::group_by(!!rlang::sym(x_obs), !!rlang::sym(y_obs)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
    dplyr::group_by(!!rlang::sym(x_obs)) %>%
    dplyr::mutate(percentage = count / sum(count) * 100)
  
  # Base ggplot
  p <- ggplot(summary_data, aes(x = !!sym(x_obs), y = percentage, fill = !!sym(y_obs))) +
    geom_bar(stat = "identity", position = "stack", width = 0.9, color = "black", linewidth = 0.2) +
    #scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    scale_fill_manual(values = pal) +
    theme_Publication(base_size = 14) + theme(axis.line = element_blank())+
    guides(fill = guide_legend(nrow = 1, keywidth = 0.6,keyheight = 0.6))+
    labs(x = cc[[1]], y = "Percentage (%)", fill = y_obs, title = cc[[2]]) 
  
  if (horizontal) {
    p_hor <- p + coord_flip()+
      theme(legend.position = "bottom",
            panel.grid.major.y = element_blank(),
            axis.ticks = element_blank(),
            #axis.text.y = element_text(angle = 0, hjust = 2, vjust = 0.5),
            axis.text.x = element_blank())+
      labs(x = cc[[1]], y = "Percentage (%)", fill = y_obs, title = cc[[2]]) 
    
    if (labell){
      p_hor = p_hor + geom_text(aes(label = ifelse(percentage >= 4, paste0(round(percentage, 1), "%"), "")),
                                position = position_stack(vjust = 0.5), size = 4)}
    return(p_hor)
  }
  else{
    p_vert = p + theme(legend.position = "bottom",
                       panel.grid.major.y = element_blank(),
                       axis.ticks = element_blank(),
                       #axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.7),
                       axis.text.y = element_blank())+
      labs(x = cc[[1]], y = "Percentage (%)", fill = y_obs, title = cc[[2]]) 
    
    if (labell){
      p_vert = p_vert + geom_text(aes(label = ifelse(percentage >= 4, paste0(round(percentage, 1), "%"), "")),
                                  position = position_stack(vjust = 0.5), size = 4)}
    return(p_vert)
  }
}

plot_prop_per_bin_overall <- function(
    adata,
    y_obs,
    y_value,
    x_obs = "spatial_bin",       # bins like Q0, Q1, etc.
    fill_obs = "TumorType",      # Pre vs Post
    split_obs = "Response",      # Responder vs Non-Responder
    t_subtype_col = "T_subtypes_l1",
    pal = tt_colors
) {
  
  df <- as.data.frame(adata$obs)
  overall_frac <- df %>%
    group_by(!!sym(split_obs), !!sym(fill_obs), !!sym(x_obs)) %>%
    summarise(
      frac = (sum(!!sym(t_subtype_col) == y_value) / n())*100,  # numerator / all T cells in bin
      .groups = "drop"
    )
  
  resp_levels <- sort(unique(overall_frac[[split_obs]]))
  plot_list <- list()
  global_limits <- range(overall_frac$frac, na.rm = TRUE)
  lower <- 0.00
  upper <- ceiling(global_limits[2] * 100) / 100
  rounded_limits <- c(lower, upper)
  print(rounded_limits)
  for(i in 1:length(resp_levels)){
    resp <- resp_levels[i]
    df_here <- overall_frac %>% filter(!!sym(split_obs) == resp)
    show_x_axis <- ifelse(i == length(resp_levels), TRUE, FALSE)
    print(show_x_axis)
    p <- ggplot(df_here, aes(x = !!sym(x_obs), y = frac, fill = !!sym(fill_obs))) +
      geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.8), linewidth = 0.3, color = "black") +
      scale_fill_manual(values = pal) +
      scale_x_discrete(labels = function(x) gsub(": ", ":\n", x)) + 
      #scale_color_manual(values = pal) +
      labs(
        x = ifelse(show_x_axis, "Spatial bin", ""), # only bottom plot shows x-axis
        y = paste0("% of ", y_value, "\n", resp),
        fill = "Tumor type"
      ) +
      theme_Publication(base_size = 12)+
      coord_cartesian(ylim = rounded_limits)+
      scale_y_continuous(expand = c(0, 0))+
      #scale_y_continuous(limits = global_limits)+
      theme(
        plot.margin = unit(c(0,5,0,5), "pt"),
        
        axis.text.x = element_text(
          angle = if (show_x_axis) 0 else 0,
          hjust = if (show_x_axis) 0.5 else 0,
          vjust = if (show_x_axis) 0.5 else 0
        ),
        
        axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
        
        axis.text.x.bottom =
          if (show_x_axis)
            element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10)
        else
          element_blank()
      )
    plot_list[[resp]] <- p
  }
  
  # Stack vertically like the reference function
  final_plot <- wrap_plots(plot_list, ncol = 1, heights = rep(1, length(plot_list)), guides = "collect") &
    theme(legend.position = "none")
  
  return(final_plot)
}
###### Dotplots ######
plot_dotplot = function(df, x_obs, y_obs, color_obs, size_obs, ll = c("col label", 'size_label')){
  p = ggplot(df, aes(x = !!sym(x_obs), y = !!sym(y_obs))) +
    geom_point(aes(size = !!sym(size_obs), color = !!sym(color_obs))) +
    #scale_color_viridis(option = "magma") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size(range = c(2, 8)) +
    theme_Publication(base_size = 14)+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_blank(),
      axis.title.x = element_blank(), 
      axis.title.y = element_blank()    
    ) + labs(color = ll[[1]], size = ll[[2]])
  return(p)
}



###### Alluvials ######
library(ggalluvial)
plot_alluvium = function(df, x, y, pal, ll = c("X", "Y", "Fill")){
  df_counts = df %>% group_by(!!sym(x), !!sym(y)) %>% summarise(count = n()) %>% 
    mutate(prop = count / sum(count)) %>%   # proportion within each TumorType
    ungroup()
  
  p1 = ggplot(data = df_counts,
              aes(x = !!sym(x), stratum = !!sym(y), alluvium = !!sym(y),
                  y = prop)) +
    geom_alluvium(aes(fill = !!sym(y)), alpha = 0.8, width = 0.7) +
    geom_stratum(aes(fill = !!sym(y)), width = 0.7, color = "black", alpha = 0) +
    #geom_stratum(width = 1/6, fill = "grey", color = "black") +
    geom_text(stat = "stratum",size = 4,
              aes(label = ifelse(prop >= 0.03, as.character(after_stat(stratum)), ""))) +
    scale_fill_manual(values = pal) +
    guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
    theme_Publication(base_size = 12)+
    labs(fill = ll[[3]], x = ll[[1]])+
    theme(legend.title = element_text(hjust = 0.5),
          #axis.title.x = element_blank(),    # remove x-axis label
          axis.title.y = element_blank(),    # remove y-axis label
          #axis.text.x = element_blank(),     # remove x-axis text
          
          axis.text.y = element_blank(),     # remove y-axis text
          axis.ticks = element_blank(),
          axis.line = element_blank(),       # remove axis lines
          panel.grid = element_blank()) +
    theme(
      plot.margin = unit(c(0.2, 0.2, 0.2, 0), "cm")  # top, right, bottom, left
    )
  return(p1)
}


###### Boxplots ######
plot_g4_paired_boxplot = function(df, celltype_of_interest, ct_col){
  plot_df <- df %>%
    group_by(Patient,TumorType) %>%
    mutate(total_cells = n()) %>%
    ungroup() %>%
    group_by(Patient, TumorType, Response) %>%
    summarise(
      total_cells = first(total_cells),
      celltype_cells = sum(!!sym(ct_col) == celltype_of_interest),
      fraction = (celltype_cells / total_cells)*100,
      .groups = "drop"
    )
  
  plot_df <- plot_df %>%
    mutate(
      group = case_when(
        TumorType == "Pre-NIVO primary" & Response == "Response"     ~ "Pre-R",
        TumorType == "Pre-NIVO primary" & Response == "Non-response"  ~ "Pre-NR",
        TumorType == "Post-NIVO primary" & Response == "Response"    ~ "Post-R",
        TumorType == "Post-NIVO primary" & Response == "Non-response" ~ "Post-NR",
        TRUE ~ "Other"  # catch any unexpected cases
      ),
      group = factor(group, levels = c("Pre-R", "Post-R", "Pre-NR","Post-NR"))
    )
  
  pre_post_R <- plot_df %>%
    filter(group %in% c("Pre-R", "Post-R")) %>%
    select(Patient, group, fraction) %>%
    tidyr::pivot_wider(names_from = group, values_from = fraction)
  pre_post_NR <- plot_df %>%
    filter(group %in% c("Pre-NR", "Post-NR")) %>%
    select(Patient, group, fraction) %>%
    tidyr::pivot_wider(names_from = group, values_from = fraction)
  wilcox_R <- wilcox.test(pre_post_R$`Pre-R`, pre_post_R$`Post-R`, paired = TRUE,exact = TRUE)
  wilcox_NR <- wilcox.test(pre_post_NR$`Pre-NR`, pre_post_NR$`Post-NR`, paired = TRUE,exact = TRUE)
  pvals_bracket <- data.frame(
    group1 = c("Pre-R", "Pre-NR"),
    group2 = c("Post-R", "Post-NR"),
    y.position = c(0.4, 0.4),  # adjust depending on your data
    label = c(
      paste0("p = ", signif(wilcox_R$p.value, 3)),
      paste0("p = ", signif(wilcox_NR$p.value, 3))
    )
  )
  
  pvals_bracket <- pvals_bracket %>%
    rowwise() %>%
    mutate(
      # Find the max fraction among the two groups in each row
      max_fraction = max(plot_df$fraction[plot_df$group %in% c(group1, group2)], na.rm = TRUE),
      y.position = max_fraction * 1.05  # 5% above the max for spacing
    ) %>%
    ungroup() %>%
    select(-max_fraction)
  
  print(pvals_bracket)
  p1 = ggplot(plot_df, aes(x = group, y = fraction) )+
    geom_boxplot(aes(fill = group), alpha = 1, outlier.shape = NA, width = 0.5) +
    scale_fill_manual(values = g4_colors)+
    geom_line(aes(group = Patient, colour = Patient), alpha = 0.8, linewidth = 0.3) +
    geom_point(aes(colour = Patient), size = 0.3) + 
    scale_color_manual(values = patient_colors)+
    labs(y = paste0("% of ", celltype_of_interest), x = NULL, fill = "Status") +
    theme_Publication(base_size = 14) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    guides(color = "none")+
    ggpubr::stat_pvalue_manual(
      pvals_bracket,
      label = "label",
      y.position = "y.position",
      tip.length = 0.01
    )
  return(p1)
}

plot_response_boxplot = function(df, groubby, obs, obs_col){
  plot_df <- df %>%
    group_by(!!sym(groubby),Response) %>%
    mutate(total_cells = n()) %>%
    ungroup() %>%
    group_by(!!sym(groubby), Response) %>%
    summarise(
      total_cells = dplyr::first(total_cells),
      total_obs = sum(!!sym(obs_col) == obs),
      fraction = (total_obs / total_cells)*100,
      .groups = "drop"
    )
  wil_pval = wilcox.test(fraction ~ Response, data = plot_df, exact = TRUE)
  print(levels(plot_df$Response))
  pvals_bracket <- data.frame(
    group1 = levels(plot_df$Response)[1],
    group2 = levels(plot_df$Response)[2],
    label = paste0("p = ", signif(wil_pval$p.value, 3)),
    y.position = 0  # placeholder, will compute next
  )
  
  pvals_bracket <- pvals_bracket %>%
    rowwise() %>%
    mutate(
      max_fraction = max(plot_df$fraction[plot_df$Response %in% c(group1, group2)], na.rm = TRUE),
      y.position = max_fraction * 1.05  # 5% above the max for spacing
    ) %>%
    ungroup() %>%
    dplyr::select(-max_fraction)
  
  p1 = ggplot(plot_df, aes(x = Response, y = fraction, fill = Response)) +
    geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.3) +
    scale_fill_manual(values = g4_colors)+
    geom_jitter(width = 0.2) +
    labs(y = paste0("% of CD8+ T cells in \n", obs), x = df$spatial_bin[[1]], fill = "Status") +
    theme_Publication(base_size = 14) +
    ggpubr::stat_pvalue_manual(
      pvals_bracket,
      label = "label",
      y.position = "y.position",
      tip.length = 0.01
    )

  return(p1)
}
###### Density Plots ######
plot_ridges <- function(df, x_col, cat_col, pal, ll = c("Overlayed Density", "Distance", "Density")) {
  # Convert column names to symbols for tidy evaluation
  x_sym <- rlang::sym(x_col)
  cat_sym <- rlang::sym(cat_col)
  
  ggplot(df, aes(x = !!x_sym, y = factor(!!cat_sym), fill = factor(!!cat_sym))) +
    #ggridges::geom_density_ridges(stat = "binline", bins = 30, scale = 1, alpha = 0.6) +
    ggridges::geom_density_ridges(alpha = 0.4, linewidth = 1) +
    scale_fill_manual(values = pal) +
    labs(title = ll[1], x = ll[2], y = ll[3]) +
    theme_Publication(base_size = 12)+
    theme(legend.position = "bottom")+
    guides(fill = guide_legend(keywidth = 0.6,keyheight = 0.6))}


plot_density_fill_split <- function(
    adata,
    y_obs = "tcell_type",
    y_value = NULL,
    x_obs = "Distance_to_Tumor",
    x_thres = 10000,
    fill_obs = "TumType",
    split_obs = "Response",
    pal,
    t_subtype_col = "T_subtypes_l1"
) {
  
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  df <- as.data.frame(adata$obs)
  
  df <- df %>%
    filter(!!sym(x_obs) < (x_thres)) %>%
    filter(!!sym(y_obs) == y_value)
  
  print(max(df$Distance_to_Tumor))
  
  # ----------------------------
  # Compute global max density
  # ----------------------------
  x_range <- range(df[[x_obs]], na.rm = TRUE)
  
  dens_list <- lapply(split(df, list(df[[split_obs]], df[[fill_obs]])), function(subdf) {
    density(subdf[[x_obs]], adjust = 0.5, from = x_range[1], to = x_range[2])
  })
  
  y_max <- max(sapply(dens_list, function(d) max(d$y))) # padding
  
  
  # ----------------------------
  # Loop over Response groups
  # ----------------------------
  plot_list <- list()
  unq <- sort(unique(df[[split_obs]]))
  
  for (i in seq_along(unq)) {
    
    grp <- unq[i]
    
    df_here <- df %>% filter(!!sym(split_obs) == grp)
    ks_grps <- unique(df_here[[fill_obs]])
    
    ks_res <- ks.test(
      df_here[[x_obs]][df_here[[fill_obs]] == ks_grps[[1]]],
      df_here[[x_obs]][df_here[[fill_obs]] == ks_grps[[2]]]
    )
    
    ks_label <- paste0("KS p = ", signif(ks_res$p.value, 3))
    grpl = ifelse(grp == "Response", "Response", "Non-response")
    y_label <- paste0("Density (", df_here[[t_subtype_col]][1], ")\n", grpl)
    
    median_df <- df_here %>%
      group_by(!!sym(fill_obs)) %>%
      do({
        x_vals <- .[[x_obs]]
        dens <- density(x_vals, adjust = 0.5)
        med_val <- median(x_vals, na.rm = TRUE)
        
        # find density value closest to median
        y_at_med <- dens$y[which.min(abs(dens$x - med_val))]
        
        data.frame(
          med = med_val,
          y_med = y_at_med
        )
      })
    
    
    p <- ggplot(df_here,
                aes(x = !!sym(x_obs),
                    fill = !!sym(fill_obs))) +
      
      geom_density(alpha = 0.8, adjust = 0.5, color = "black", linewidth = 0.5, trim = FALSE) +
      scale_fill_manual(values = pal) +
      
      coord_cartesian(
        xlim = c(0, x_thres-1000),
        ylim = c(0, y_max * 1.05)   # increased slightly to avoid cut-off
      ) +
      
      labs(
        x = ifelse(i == length(unq), "Distance to tumor", ""),  # only bottom plot shows x-axis
        y = y_label
      ) +
      
      annotate(
        "text",
        x = (x_thres-1000) * 0.85,
        y = y_max * 1.05,
        label = ks_label,
        hjust = 1,
        size = 4
      ) + 
      geom_segment(
        data = median_df,
        aes(x = med, xend = med, y = 0, yend = y_med),
        color = "black",
        linewidth = 0.5,
        linetype = "dashed",
        inherit.aes = FALSE
      ) +
      
      
      theme_Publication(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        legend.title = element_blank(),
        plot.margin = unit(c(0, 5, 0, 5), "pt")  # minimal spacing
      )
    
    plot_list[[grp]] <- p
  }
  
  # ----------------------------
  # Stack plots vertically with minimal gap
  # ----------------------------
  final_plot <- wrap_plots(plot_list, ncol = 1, heights = rep(1, length(plot_list)), guides = "collect") &
    theme(legend.position = "none")
  
  return(final_plot)
}




###### Heatmap ######
library(viridis)
plot_gene_distance_heatmap <- function(adata, genes, cbreaks, max_distance = 5000,tit = "Pheatmap",
                                       distance_obs = "Distance_to_Tumor" ) {
  
  adata = adata[adata$obs[[distance_obs]]<max_distance]
  df <- as.data.frame(adata$obs)
  expr <- as.data.frame(as.matrix(adata$X))
  colnames(expr) <- adata$var_names$to_list()
  print(head(expr))
  expr <- expr[, genes, drop = FALSE]
  df_expr <- cbind(df, expr)
  df_expr <- df_expr %>% arrange(!!sym(distance_obs))
  
  n_bins <- length(cbreaks) - 1
  first_label <- paste0("[", cbreaks[1], ", ", cbreaks[2], "]")
  other_labels <- paste0("(", cbreaks[2:n_bins], ", ", cbreaks[3:(n_bins+1)], "]")
  labels <- c(first_label, other_labels)
  
  df_expr <- df_expr %>%
    mutate(distance_bin = cut(
      !!sym(distance_obs),
      breaks = cbreaks,
      labels = labels,
      include.lowest = TRUE,
      right = TRUE 
    ))
  print(table(df_expr$distance_bin))
  bin_avg <- df_expr %>%
    group_by(distance_bin) %>%
    summarise(across(all_of(genes), mean, na.rm = TRUE))
  
  mat <- as.matrix(bin_avg[, genes])
  rownames(mat) <- bin_avg$distance_bin
  #mat_scaled <- t(scale(mat)) 
  my_colors <- viridis(1000)
  mat_scaled <- t(apply(t(mat), 1, function(x) (x - min(x)) / (max(x) - min(x))))
  #mat_scaled = t(apply(t(mat), 1, function(x) x))
  
  
  
  #mat_scaled = t(apply(t(mat), 1, function(x) (x/max(x))))
  p1 = pheatmap::pheatmap(
    mat_scaled,
    color = my_colors,       # same color scale
    #breaks = c(
    #  seq(0, 0.25, length.out = 100),
    #  seq(0.25, 1, length.out = 901)[-1]
    #),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    border_color = NA,           # removes borders like raster in ComplexHeatmap
    fontsize_row = 8,           # adjust as needed
    fontsize_col = 8,
    angle_col = 90,
    silent = TRUE,               # make it plot immediately
    main = tit
  )
  return(p1)
}


plot_mirrored_lollipop = function(adata, genes, mirror_cond, threshold, cbreaks = seq(0, 500, by = 20),
                                  distance_obs = "Distance_to_Tumor", pal = response_colors, cc = c("Distance to Tumor",
                                                                                                    "", "Response", "Expression",
                                                                                                    "Lollipop")){
  df <- as.data.frame(adata$obs)
  expr <- as.data.frame(as.matrix(adata$X))
  colnames(expr) <- adata$var_names$to_list()
  expr <- expr[, genes, drop = FALSE]
  df_expr <- cbind(df, expr)
  df_expr <- df_expr %>% arrange(!!sym(distance_obs))
  
  n_bins <- length(cbreaks) - 1
  first_label <- paste0("[", cbreaks[1], ", ", cbreaks[2], "]")
  other_labels <- paste0("(", cbreaks[2:n_bins], ", ", cbreaks[3:(n_bins+1)], "]")
  labels <- c(first_label, other_labels)
  
  df_expr <- df_expr %>%
    mutate(distance_bin = cut(
      !!sym(distance_obs),
      breaks = cbreaks,
      labels = labels,
      include.lowest = TRUE,
      right = TRUE 
    ))
  print(table(df_expr$distance_bin))
  bin_avg <- df_expr %>%
    group_by(distance_bin, !!sym(mirror_cond)) %>%
    summarise(across(all_of(genes), mean, na.rm = TRUE))
  
  conds = unique(bin_avg[[mirror_cond]])
  
  bin_avg_1 = bin_avg[bin_avg[[mirror_cond]] == conds[[1]],]
  mat_1 <- as.matrix(bin_avg_1[, genes])
  rownames(mat_1) <- paste0(bin_avg_1$distance_bin, ":", conds[[1]])
  mat_scaled_1 <- t(apply(t(mat_1), 1, function(x) (x - min(x)) / (max(x) - min(x))))
  #mat_scaled_1 <- t(apply(t(mat_1), 1, function(x) (x)))
  
  bin_avg_2 = bin_avg[bin_avg[[mirror_cond]] == conds[[2]],]
  mat_2 <- as.matrix(bin_avg_2[, genes])
  rownames(mat_2) <- paste0(bin_avg_2$distance_bin, ":", conds[[2]])
  mat_scaled_2 <- t(apply(t(mat_2), 1, function(x) (x - min(x)) / (max(x) - min(x))))
  #mat_scaled_2 <- t(apply(t(mat_2), 1, function(x) (x)))
  mat_scaled = cbind(mat_scaled_1, mat_scaled_2)
  
  df_long <- as.data.frame(mat_scaled) %>%
    tibble::rownames_to_column("marker") %>%
    pivot_longer(-marker, names_to = "distance_bin", values_to = "expression")
  df_long <- df_long %>%
    separate(distance_bin, into = c("distance_bin", "condition"), sep = ":")
  df_long <- df_long %>%
    mutate(distance_mid = sapply(distance_bin, function(x) {
      nums <- as.numeric(unlist(regmatches(x, gregexpr("\\d+", x))))
      mean(nums)
    }))
  
  lolli_df <- df_long %>%
    filter(expression >= threshold)
  
  lolli_df <- lolli_df %>%
    mutate(dist_mirror = ifelse(condition == conds[[2]], -distance_mid, distance_mid))
  lolli_df$marker = factor(lolli_df$marker, levels = rev(genes))
  max_dist <- max(df_long$distance_mid)
  p = ggplot(lolli_df, aes(x = dist_mirror, y = marker)) +
    geom_segment(aes(xend = 0, yend = marker, color = condition), size = 0.8, alpha = 0.7) +
    geom_point(aes(size = expression, color = condition), alpha = 0.9) +
    geom_vline(xintercept = 0, color = "black", size = 0.8) +
    scale_color_manual(values = pal)+
    scale_size_continuous(range = c(2,6)) +
    scale_x_continuous(limits = c(-max_dist, max_dist), labels = function(x) abs(x)) +
    theme_Publication(base_size = 14) +
    labs(x = cc[[1]], y = cc[[2]],
         color = cc[[3]], size = cc[[4]],
         title = cc[[5]])+
    theme(
      axis.ticks.y = element_blank(),       # remove y ticks
      axis.line.y = element_blank(),        # remove y labels
      panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"), # horizontal grid
      panel.grid.minor.y = element_blank(), # remove minor y grid
      panel.grid.major.x = element_blank(), # optional: remove vertical grid
      panel.grid.minor.x = element_blank()#,legend.title = element_text(hjust = 0.5)
    )
  
  return(p)
}


source("/Users/lakshmi_nccs/Desktop/NCCS/Projects/Helpers/Library/spatialshells.R")
plot_spatialshells_wrapper = function(adata, gn,thresh,bb=seq(0, 500, by = 30), scaling = "min-max"){
  
  adata = adata[adata$obs$TumorType %in% c("Pre-NIVO Primary", "Post-NIVO Primary"), ]
  df = adata$obs%>% dplyr::select(Distance_to_Tumor, TumorType, Response, tma_sid)
  df[['group']] = paste0(df$TumorType , "_", df$Response)
  df[gn] <- adata$X[, match(gn, adata$var_names$to_list())]
  df_long <- df %>%
    pivot_longer(
      cols = all_of(c(gn)),
      names_to = "genes",     # this will be your level-2 variable
      values_to = "value"      # the corresponding measurement
    )
  df_long <- df_long %>%
    mutate(
      group = case_when(
        group == "Pre-NIVO Primary_Responder"    ~ "Pre-R",
        group == "Pre-NIVO Primary_Non_Responder" ~ "Pre-NR",
        group == "Post-NIVO Primary_Responder"   ~ "Post-R",
        group == "Post-NIVO Primary_Non_Responder"~ "Post-NR",
        TRUE                                  ~ group  # keep others unchanged
      )
    )
  df_long$group = factor(df_long$group, levels = c("Pre-R", "Post-R", "Post-NR", "Pre-NR"))
  df_long$genes = factor(df_long$genes, levels = gn)
  p1 = make_spatial_shell_plot(
    df = df_long[df_long$Distance_to_Tumor < thresh,],
    dist_col = "Distance_to_Tumor",
    level1_col = "group",
    level2_col = "genes",
    value_col = "value",
    cut_type = "custom",
    n_bins = 50,  # numeric bin width, number of bins, or vector of breaks
    custom_breaks = bb,
    scaling = scaling,
    title = gn,
    fill_name = "Expression",
    viridis_palette = "viridis"
  ) 
  return(p1)
}


plot_spatialshells_wrapper_explants = function(adata, gn, dist = "Distance_norm", bb = seq(0,1, by= 0.05), scaling = "min-max"){
  
  df = adata$obs%>% select(!!sym(dist), Patient, Response)
  df[['group']] = paste0(df$Patient , "_", df$Response)
  df[gn] <- adata$X[, match(gn, adata$var_names$to_list())]
  df_long <- df %>%
    pivot_longer(
      cols = all_of(c(gn)),
      names_to = "genes",     # this will be your level-2 variable
      values_to = "value"      # the corresponding measurement
    )
  df_long$genes = factor(df_long$genes, levels = gn)
  p1 = make_spatial_shell_plot(
    df = df_long,
    dist_col = dist,
    level1_col = "group",
    level2_col = "genes",
    value_col = "value",
    cut_type = "custom",
    n_bins = 50,  # numeric bin width, number of bins, or vector of breaks
    custom_breaks = bb,
    scaling = scaling,
    title = gn,
    fill_name = "Expression",
    viridis_palette = "viridis"
  ) 
  return(p1)
}


###### Spatial Recon ######
spatial_plot_6k_ii = function(adata, tma, plotby, phencol, cc, ii){
  adata = adata[adata$obs['TMA'] == paste0("TMA",tma)]$copy()
  poly_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data/"
  pol1 = readRDS( paste0(poly_dir, "/tma", tma, "_polygon.rds"))
  adata_selec = adata$obs %>% dplyr::select(cell, Patient,tma_sid,TumorType, SID,!!sym(phencol))
  pol1 = pol1%>% right_join(adata_selec, by = "cell")
  pol11 = pol1[pol1[[plotby]] == ii,]
  print(phencol)
  p1 = ggplot(pol11, aes(x = x_global_px, y = y_global_px, group = cell, fill = !!sym(phencol))) +
    geom_polygon(color = "black", alpha = 0.8, size = 0.4) +
    geom_path(color = "black", alpha = 0.4, size = 0.1) +
    scale_fill_manual(values = cc, drop = TRUE) + 
    theme_Publication(base_size = 14) +
    labs(title = paste0(pol11$Patient[[1]], ": ", pol11$TumorType[[1]], ": ", plotby, ": ", ii),
         x = "Global X Position (px)",
         y = "Global Y Position (px)",
         fill = "Cell type") +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )+
    theme(legend.position = "right")
  #+guides(fill = guide_legend(ncol = 1, keywidth = 1,keyheight = 1))
  return(p1)
  #plot_save(p1, paste0(plotby, "_", ii), 10,16)
}

spatial_plot_comet_ii = function(adata, tma, plotby, phencol, cc, ii){
  adata = adata[adata$obs['TMA'] == paste0("TMA",tma)]$copy()
  poly_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx6k/processed_data"
  pol1 = readRDS( paste0(poly_dir, "/comet_tma", tma, "_polygon.rds"))
  adata_selec = adata$obs %>% select(CellId, Patient,TumorType,!!sym(phencol))
  print(dim(adata_selec))
  pol1 = pol1%>% right_join(adata_selec, by = "CellId")
  pol11 = pol1[pol1[[plotby]] == ii,]
  print(phencol)
  p1 = ggplot(pol11, aes(x = XCoordinate, y = YCoordinate, group = CellId, fill = !!sym(phencol))) +
    geom_polygon(color = "black", alpha = 0.7, size = 0.5) +
    geom_path(color = "black", alpha = 0.4, size = 0.5) +
    scale_fill_manual(values = cc, drop = FALSE) + 
    theme_Publication(base_size = 14) +
    labs(title = paste0(pol11$Patient[[1]], ": ", pol11$TumorType[[1]], ": ", plotby, ": ", ii),
         x = "Global X Position (px)",
         y = "Global Y Position (px)",
         fill = "Cell type") +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )+
    theme(legend.position = "right", legend.title = element_text(hjust = 0.5))
  #+guides(fill = guide_legend(ncol = 1, keywidth = 1,keyheight = 1))
  return(p1)
  #plot_save(p1, paste0(plotby, "_", ii), 10,16)
}


spatial_plot_ii_1k = function(adata, samidnum, plotby, phencol, cc, ii){
  adata = adata[adata$obs['samid'] == paste0("s",samidnum)]$copy()
  poly_dir = "/Users/lakshmi_nccs/Desktop/NCCS/Projects/Cosmx1k/processed_data/"
  pol1 = readRDS( paste0(poly_dir, "/s", samidnum, "_polygon.rds"))
  pol1$uniqid = paste0("s",samidnum, "_", pol1$fov, "_", pol1$cell_ID)
  adata_selec = adata$obs %>% select(uniqid, patient,tma_sid,phenotype, fov,!!sym(phencol), !!sym(plotby))
  pol1 = pol1%>% right_join(adata_selec, by = "uniqid")
  pol11 = pol1[pol1[[plotby]] == ii,]
  print(phencol)
  p1 = ggplot(pol11, aes(x = x_global_px, y = y_global_px, group = uniqid, fill = !!sym(phencol))) +
    #geom_point(size = 0.5)+
    geom_polygon(color = "black", alpha = 0.7, size = 0.5) +
    geom_path(color = "black", alpha = 0.4, size = 0.5) +
    #scale_fill_manual(values = cc, drop = FALSE) + #keeps all
    scale_fill_manual(values = cc) + 
    theme_minimal(base_size = 14) +
    labs(title = paste0(pol11$patient[[1]], ": ", pol11$phenotype[[1]], ": ", plotby, ": ", ii),
         x = "Global X Position (px)",
         y = "Global Y Position (px)",
         fill = "Cell type") +
    theme(text = element_text(family = "Arial"),
          legend.position = "bottom",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_blank(),
          legend.margin = margin(t = -15),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.box.background = element_rect(fill = "white", color = NA))
  #+guides(fill = guide_legend(ncol = 1, keywidth = 1,keyheight = 1))
  return(p1)
  #plot_save(p1, paste0(plotby, "_", ii), 10,16)
}