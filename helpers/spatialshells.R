library(dplyr)
library(ggplot2)

make_spatial_shell_plot <- function(
    df,
    dist_col = "Distance_to_Tumor",
    level1_col = "level1",
    level2_col = "level2",
    value_col = "value",
    cut_type = c("cut", "pcut", "custom"),
    n_bins = 50,  # numeric bin width, number of bins, or vector of breaks
    custom_breaks = NULL,
    scaling = c("min-max", "z-score","none"),
    title = "Spatial Shell Plot",
    fill_name = "Value",
    col_palette = NULL,
    viridis_palette = c("viridis", "magma", "plasma", "inferno", "cividis", "turbo")
) {
  
  library(dplyr)
  library(ggplot2)
  df_in <- df %>%
    rename(
      distance = !!sym(dist_col),
      level1   = !!sym(level1_col),
      level2   = !!sym(level2_col),
      value    = !!sym(value_col)
    )
  
  if (!is.numeric(df_in$distance)) {
    # Distance is categorical → use as-is
    df_in <- df_in %>%
      mutate(
        distance_bin = factor(distance, levels = unique(distance)))
    df_in <- df_in %>%
      mutate(distance_mid = sapply(distance_bin, function(x) {
        nums <- as.numeric(unlist(regmatches(x, gregexpr("\\d+", x))))
        mean(nums)
      }))
    
  }
  else{
    cut_type <- match.arg(cut_type)
    # --- BINNING ---
    if (cut_type == "cut") {
      # Equal-width bins based on number of bins
      min_dist <- min(df_in$distance, na.rm = TRUE)
      max_dist <- max(df_in$distance, na.rm = TRUE)
      breaks <- seq(min_dist, max_dist, length.out = n_bins + 1)
      labels <- paste0("[", round(head(breaks, -1),1), "-", round(tail(breaks, -1),1), ")")
      df_in <- df_in %>%
        mutate(distance_bin = cut(distance, breaks = breaks,labels = labels,
                                  include.lowest = TRUE, right = FALSE))
      bin_midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2
      df_in$distance_mid <- bin_midpoints[as.numeric(df_in$distance_bin)]
      
    } else if (cut_type == "pcut") { 
      # Equal-count bins (percentiles) 
      probs <- seq(0, 1, length.out = n_bins + 1) 
      breaks <- quantile(df_in$distance, probs = probs, na.rm = TRUE) 
      labels <- paste0("[", round(head(breaks, -1),1), "-", round(tail(breaks, -1),1), ")")
      df_in <- df_in %>% mutate(distance_bin = cut(distance, breaks = breaks, labels = labels,
                                                   include.lowest = TRUE, right = FALSE))
      bin_midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2
      df_in$distance_mid <- bin_midpoints[as.numeric(df_in$distance_bin)]
    } else if (cut_type == "custom") {
      # Custom breaks
      labels <- paste0("[", round(head(custom_breaks, -1),1), "-", round(tail(custom_breaks, -1),1), ")")
      df_in <- df_in %>%
        mutate(distance_bin = cut(distance, breaks = c(custom_breaks, Inf) ,labels = c(labels, paste0("[", round(tail(custom_breaks, 1),1), "+)")),
                                  include.lowest = TRUE, right = FALSE))
      print("Here")
      print(unique(df_in$distance_bin))
      bin_midpoints <- (head(custom_breaks, -1) + tail(custom_breaks, -1)) / 2
      df_in$distance_mid <- bin_midpoints[as.numeric(df_in$distance_bin)]
    }
  }
  
  print(levels(df_in$level1))
  print(levels(df_in$level2))
  x_levels <- df_in %>%
    select(level1, level2) %>%
    distinct() %>%
    arrange(level1, level2) %>%   # sort by level1 factor first, then level2 factor
    mutate(x_axis = paste(level1, level2, sep = "::")) %>%
    pull(x_axis)                # extract as vector
  
  print(x_levels)
  
  df_in$x_axis <- factor(paste(df_in$level1, df_in$level2, sep="::"), levels = x_levels)
  
  df_sum <- df_in %>%
    group_by(x_axis, distance_bin) %>%
    summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")



  scaling = match.arg(scaling)
  if (scaling == "min-max") {
    print("Scaling here check")
    df_sum <- df_sum %>%
      group_by(x_axis) %>%
      mutate(scaled_value = (mean_value - min(mean_value)) / (max(mean_value) - min(mean_value))) %>%
      ungroup()
  } else if (scaling == "z-score") {
    df_sum <- df_sum %>%
      group_by(x_axis) %>%
      mutate(scaled_value = (mean_value - mean(mean_value)) / sd(mean_value)) %>%
      ungroup()
  } else {
    df_sum <- df_sum %>%
      group_by(x_axis) %>%
      mutate(scaled_value = mean_value) %>%
      ungroup()
    
    #df_sum$scaled_value <- df_sum$mean_value
  }
  df_sum_complete <- df_sum %>%
    complete(
      x_axis,
      distance_bin
    )%>%
    mutate(
      mean_value = ifelse(is.na(mean_value), 0, mean_value),   # replace NA with 0
      scaled_value = ifelse(is.na(scaled_value), 0, scaled_value)  # if you already have scaled_value
    )
  
  viridis_palette <- match.arg(viridis_palette)
  
  # --- PLOT ---
  p <- ggplot(df_sum_complete, aes(
    x    = x_axis,
    y    = as.numeric(distance_bin),
    fill = scaled_value
  )) +
    geom_tile(color = "black", linewidth = 0.3) +
    coord_polar(theta = "x", start = 3*pi/2) +
    scale_fill_viridis_c(option = viridis_palette, name = fill_name, na.value = "grey") +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = title)
  
  
  x_numeric <- seq_along(levels(df_sum_complete$x_axis))  # numeric positions of segments
  n_segments <- length(x_numeric)
  label_df <- data.frame(
    x = x_numeric,
    y = max(as.numeric(df_sum_complete$distance_bin)) * 1.1,
    label = levels(df_sum_complete$x_axis)
  )
  
  label_df$angle <- (360 / n_segments) * (label_df$x - 0.5)+ 270
  label_df$tangent_angle <- 180 - label_df$angle
  label_df$tangent_wrapped <- ((label_df$tangent_angle + 180) %% 360) - 180
  flip_idx <- abs(label_df$tangent_wrapped) > 90
  label_df$angle_final <- label_df$tangent_wrapped
  label_df$angle_final[flip_idx] <- label_df$angle_final[flip_idx] + 180
  
  p1 = p +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label, angle = angle_final, hjust = 0.5),  # center text on midpoint
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4
    )
  return(list(p1,df_sum))
}