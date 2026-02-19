library(patchwork)
library(magick)
library(ggplot2)
library(grid)

# 1. Setup paths for the new files
base_path <- "./figures/"
loocv_files <- c(
  "LOOCV_Observed_vs_Predicted_PersonalGrowth.tiff",
  "LOOCV_Residuals_PersonalGrowth.tiff"
)

# 2. Re-using your logic to process images
# I've slightly modified the title cleanup for these specific names
process_loocv <- function(file_name) {
  # Create a clean title (e.g., "Observed vs Predicted")
  clean_title <- gsub("LOOCV_", "", file_name)
  clean_title <- gsub("_PersonalGrowth.tiff", "", clean_title)
  clean_title <- gsub("_", " ", clean_title)
  
  img <- magick::image_read(paste0(base_path, file_name))
  grob <- grid::rasterGrob(img)
  
  wrap_elements(grob) + 
    ggtitle(clean_title) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
}

# 3. Process the two images
loocv_plots <- lapply(loocv_files, process_loocv)
# Modifiche per ridurre lo spazio tra i grafici e ingrandire le etichette 'a' e 'b'

final_loocv_plot <- wrap_plots(loocv_plots) + 
  plot_layout(ncol = 1, nrow = 2) & # Usiamo '&' per applicare il tema a tutti i subplot
  theme(
    # Ingrandisce le etichette delle sottobonifere (a, b)
    plot.tag = element_text(size = 22, face = "bold"),
    # Riduce i margini di ogni singolo plot per avvicinarli
    # I valori sono: top, right, bottom, left
    plot.margin = margin(t = 2, r = 5, b = 2, l = 5) 
  )

# Aggiungiamo le annotazioni finali
final_loocv_plot <- final_loocv_plot + 
  plot_annotation(
    title = "", 
    tag_levels = 'a',
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

# Per ridurre ulteriormente lo spazio verticale tra i due plot nel layout
final_loocv_plot <- final_loocv_plot + 
  plot_layout(guides = 'collect') & 
  theme(plot.spacing = unit(0.1, "lines")) # Riduce lo spazio vuoto tra i pannelli


# 5. Display
print(final_loocv_plot)

# 6. Save (Adjusting dimensions for a vertical 2x1 layout)
ggsave(paste0(base_path, "Combined_LOOCV_Diagnostics_PersonalGrowth.tiff"), 
       final_loocv_plot, width = 8, height = 12, dpi = 600, compression = "lzw")