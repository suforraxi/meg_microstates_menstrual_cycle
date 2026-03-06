# 1. Load libraries
library(patchwork)
library(png)
library(grid)
library(ggplot2)

# 2. Define path
img_folder <- "/Users/matte/Desktop/git_rep/meg_microstates_menstrual_cycle/src/r/figures/collapse_figure/"

# 3. Helper Function
prep_img <- function(path, title_text) {
  img <- readPNG(path)
  g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + 
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    labs(title = title_text) +
    theme_void() + 
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(5, 5, 5, 5, "pt"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
}

# 4. Load the images
p1 <- prep_img(paste0(img_folder, "Maximum_F_Plot.png"), "Maximum F-statistic and significant clusters")
p3 <- prep_img(paste0(img_folder, "microstate_map_00.png"), "HDM 0 of the optimal K=13 solution")
p4 <- prep_img(paste0(img_folder, "0_heat_map_matrices.png"), "Correlation HDM 0 (K=13) with HDMs across K")
p5 <- prep_img(paste0(img_folder, "microstate_map_01.png"), "HDM 1 of the optimal K=13 solution")
p6 <- prep_img(paste0(img_folder, "1_heat_map_matrices.png"), "Correlation HDM 1 (K=13) with HDMs across K")

# 5. Define the Custom Design for Rows 2 and 3
# 'A' is the map, 'B' is the matrix, '#' is empty space
# This creates a 3x2 grid where the matrix spans the entire second column
row_design <- "
  # B
  A B
  # B
"

# 6. Assemble Sub-Grids
# Sub-grid for HDM 0
row2 <- p4 + p3  + 
  plot_layout(design = row_design, widths = c(1, 3))

# Sub-grid for HDM 1
row3 <- p6 + p5 + 
  plot_layout(design = row_design, widths = c(1, 3))

# 7. Final Assembly
# We stack them vertically and apply auto-tagging
final_figure <- (p1 / row2 / row3) + 
  plot_layout(heights = c(1.8, 2.1, 2.1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    plot.margin = margin(10, 10, 10, 10, "pt")
  )

# 8. Save
ggsave(paste0(img_folder, "Combined_Microstates_Figure.png"), 
       plot = final_figure, width = 16, height = 26, dpi = 600)