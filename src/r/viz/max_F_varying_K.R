# 6. Create the Plot
p <- ggplot() +
  # Line 1: Max_F (Increased size for visibility)
  geom_line(data = spline_f, aes(x = x, y = y, color = "Max F-statistic"), 
            size = 1.5, alpha = 0.8) +
  
  # Line 2: n_sig_states (Increased size)
  geom_line(data = df, aes(x = K, y = n_sig_states * scale_factor, color = "Sig. Clusters (FDR)"), 
            size = 1.5, linetype = "dashed", alpha = 0.8) +
  
  # Highlight Maximum F (Larger point)
  geom_point(data = df[df$Max_F == max(df$Max_F),], aes(x = K, y = Max_F), 
             color = "#e74c3c", size = 6, shape = 18) +
  
  scale_color_manual(name = NULL, 
                     values = c("Max F-statistic" = "#2c3e50", 
                                "Sig. Clusters (FDR)" = "#27ae60")) +
  
  # Axis Configuration
  scale_x_continuous(breaks = seq(min(df$K), max(df$K), by = 1)) +
  scale_y_continuous(
    name = "Maximum F-statistic",
    n.breaks = 10,
    sec.axis = sec_axis(~./scale_factor, 
                        name = "Number of significant clusters (FDR)",
                        breaks = c(0, 1, 2)) 
  ) +
  
  # Styling - ADJUSTED TEXT SIZES
  labs(title = "", x = "K (Number of Clusters)") +
  theme_minimal() +
  theme(
    # Axis Titles (X and Y labels)
    axis.title.x = element_text(size = 18, margin = margin(t = 10), face = "bold"),
    axis.title.y = element_text(size = 18, margin = margin(r = 10), face = "bold", color = "#2c3e50"),
    axis.title.y.right = element_text(size = 18, margin = margin(l = 10), color = "#27ae60", face = "bold"),
    
    # Axis Ticks (The numbers)
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(size = 14, face = "bold", color = "#2c3e50"),
    axis.text.y.right = element_text(size = 14, face = "bold", color = "#27ae60"),
    
    # Tick marks themselves
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    
    # Legend text
    legend.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    
    panel.grid.minor = element_blank()
  )

# 7. Show and Save
print(p)
# Using a higher width/height ratio often helps with panel alignment later
ggsave("./figures/collapse_figure/Maximum_F_Plot.png", plot = p, width = 12, height = 7, dpi = 600)