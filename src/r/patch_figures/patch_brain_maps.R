# 1. Load necessary libraries
# install.packages(c("patchwork", "magick", "ggplot2"))
library(patchwork)
library(magick)
library(ggplot2)

# 2. Set the path and list the files
img_path <- "/Users/matte/Desktop/microstate_plots_abs/"
# We list only .png files and sort them to ensure Microstate 1 comes before 10
files <- list.files(img_path, pattern = "\\.png$", full.names = TRUE)
files <- sort(files) # Standard alphanumeric sort

# 3. Read images and convert them to patchwork-compatible objects
# We use wrap_elements to treat the image as a plot layer
img_list <- lapply(files, function(f) {
  img <- image_read(f)
  wrap_elements(panel = grid::rasterGrob(img))
})


# 1. Create a 3x4 grid for the first 12
#top_grid <- wrap_plots(img_list[1:12], ncol = 4)

# 2. Create a row for the 13th plot with spacers to push it to the middle
# Using 1 spacer on each side of the plot in a 3-column sub-layout
#bottom_row <- plot_spacer() + img_list[[13]] + plot_spacer() + plot_layout(widths = c(1, 1, 1))

# 3. Stack them
#final_grid_centered <- (top_grid / bottom_row) + 
#  plot_layout(heights = c(3, 1)) # The top takes 3 parts, bottom takes 1

# Define a helper to extract images by their "Microstate Number" 
# (Assuming img_list[1] is MS1, img_list[2] is MS2, etc.)
ms <- function(n) img_list[[n]]

# 2. Construct each row individually
# Row 1 (4 images): Naturally fills the width
row1 <- ms(1) + ms(11) + ms(8) + ms(5) + plot_layout(ncol = 4)

# Row 2 (3 images): Centered using spacers
# We use a 5-column logic [spacer, img, img, img, spacer] to force the 3 maps into the center
row2 <- ms(12) + ms(2) + ms(3) + ms(9)  + plot_layout(ncol = 4)

# Row 3 (3 images): Centered
row3 <- plot_spacer() + ms(10) + ms(4)  + plot_spacer() + plot_layout(ncol = 4, widths = c(0.5, 1, 1, 0.5))

# Row 4 (3 images): Centered
row4 <- plot_spacer() + ms(6) + ms(7) + ms(13) + plot_spacer() + plot_layout(ncol = 5, widths = c(0.5, 1, 1, 1, 0.5))

# 3. Assemble the full figure
final_composite <- (row1 / row2 / row3 / row4) + 
  plot_annotation(
    title = "",
    #tag_levels = list(c("1", "11", "8", "5", "2", "9", "3", "10", "4", "12", "5", "7", "13"))
  )

# 4. Save
ggsave("/Users/matte/Desktop/git_rep/women_cycle/src/r/figure/FINAL_DESIGN_CENTERED_4_4_2.tiff", 
       plot = final_composite, width = 12, height = 14, dpi = 600, compression = "lzw")

print(final_composite)
