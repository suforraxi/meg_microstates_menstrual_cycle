library(patchwork)
library(magick)
library(ggplot2)
library(grid)

# 1. Define paths
path7 <- "/Users/matte/Desktop/git_rep/women_cycle/src/r/figure/Microstate_7_Plot_across_phases.tiff"
path2 <- "/Users/matte/Desktop/git_rep/women_cycle/src/r/figure/Microstate_2_Plot_across_phases.tiff"

# 2. Read images
img7 <- magick::image_read(path7)
img2 <- magick::image_read(path2)

# 3. Convert to graphical objects
g7 <- grid::rasterGrob(img7)
g2 <- grid::rasterGrob(img2)

# 4. Create the patchwork with sub-titles
# We use wrap_elements and add a theme-based title to each
p1 <- wrap_elements(g7) + 
  ggtitle("Microstate map 7") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

p2 <- wrap_elements(g2) + 
  ggtitle("Microstate map 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# 5. Combine and annotate
final_plot <- p1 + p2 + 
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# 6. Display
print(final_plot)

# Optional: Save with high resolution for your paper
ggsave("/Users/matte/Desktop/git_rep/women_cycle/src/r/figure/Combined_Maps_violin_7_2.png", 
        final_plot, width = 12, height = 6, dpi = 600)