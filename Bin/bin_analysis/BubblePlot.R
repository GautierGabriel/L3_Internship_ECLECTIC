GOBubble_module <- function(data, colors, display = "single", title = "", labels = FALSE, bg.col = FALSE) {
  
  # Data preparation
  colnames(data) <- tolower(colnames(data))
  if (!"count" %in% colnames(data)) {
    data$count <- rep(1, nrow(data))
  }
  data$adj_pval <- -log10(data$adj_pval)
  sub <- data[!duplicated(data$term), ]
  
  # Base plot setup
  g <- ggplot(sub, aes(zscore, adj_pval, fill = module_names, size = count)) + 
    geom_point(shape = 21, col = "black", alpha = 0.5) + 
    geom_hline(yintercept = 1.3, col = "orange") + 
    scale_size(range = c(3, 30), guide = "none") +  # Set minimum size to 3
    labs(title = title, x = "z-score", y = "-log (adj p-value)") +
    scale_fill_manual(values = unique(combined_data$module_names)) +  # Manually set colors
    theme(legend.position = "bottom")
  
  # Subset data for labels
  sub2 <- if (!is.character(labels)) {
    subset(sub, adj_pval >= labels)
  } else {
    subset(sub, id %in% labels | term %in% labels)
  }
  
  # Add labels
  g <- g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = term), size = 4)
  
  # General theme settings
  g <- g + theme(axis.text = element_text(size = 14), 
                 axis.line = element_line(colour = "grey80"),
                 axis.ticks = element_line(colour = "grey80"), 
                 axis.title = element_text(size = 14, face = "bold"),
                 panel.background = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_line(colour = "grey80"), 
                 plot.background = element_blank())
  
  return(g)
}