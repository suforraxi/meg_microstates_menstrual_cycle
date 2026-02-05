library(patchwork)
library(magick)
library(ggplot2)
library(grid)

# 1. Setup paths and filenames
base_path <- "./figures/violins/"
files <- c(
  "Autonomy_Plot_across_phases.tiff",
  "EnvironmentalMastery_Plot_across_phases.tiff",
  "PersonalGrowth_Plot_across_phases.tiff",
  "PositiveRelationswithOthers_Plot_across_phases.tiff",
  "PurposeinLife_Plot_across_phases.tiff",
  "Self_Acceptance_Plot_across_phases.tiff"
)

# 2. Function to process each image into a patchwork-ready object
process_img <- function(file_name) {
  # Create a clean title (e.g., "Personal Growth")
  clean_title <- gsub("_Plot_across_phases.tiff", "", file_name)
  clean_title <- gsub("_", " ", clean_title)
  
  # Read and convert
  img <- magick::image_read(paste0(base_path, file_name))
  grob <- grid::rasterGrob(img)
  
  # Wrap and add title
  wrap_elements(grob) + 
    ggtitle(clean_title) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}

# 3. Process all images into a list
plot_list <- lapply(files, process_img)

# 4. Combine with patchwork (3 columns x 2 rows)
final_wellbeing_plot <- wrap_plots(plot_list) + 
  plot_layout(ncol = 3, nrow = 2) +
  plot_annotation(
    title = "Wellbeing Subscales Across Menstrual Cycle Phases",
    tag_levels = 'A', # Adds A, B, C, D... to panels automatically
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

# 5. Display
print(final_wellbeing_plot)

# 6. Save as high-resolution TIFF for publication
ggsave("./figures/Combined_Wellbeing_Subscales_3x2.tiff", 
       final_wellbeing_plot, width = 15, height = 10, dpi = 600, compression = "lzw")