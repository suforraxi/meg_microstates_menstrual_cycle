library(patchwork)
library(magick)
library(ggplot2)
library(grid)

library(patchwork)
library(magick)
library(ggplot2)
library(grid)

# 1. Definizione percorsi
path7_v <- "./figures/violins/Microstate_z_7_Plot_across_phases.tiff"
path2_v <- "./figures/violins/Microstate_Z_2_Plot_across_phases.tiff"
path7_m <- "./figures/brains_7_2/microstate_map_07.png"
path2_m <- "./figures/brains_7_2/microstate_map_02.png"

# 2. Funzione per leggere e processare le immagini (con auto-trim)
prepare_img <- function(path, title) {
  img <- magick::image_read(path)
  img <- magick::image_trim(img) # Toglie i bordi bianchi extra
  grob <- grid::rasterGrob(img, interpolate = TRUE)
  
  wrap_elements(grob) + 
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
}

# ... (Parti 1 e 2 del tuo script rimangono invariate) ...

# 2. Funzione aggiornata per minimizzare i margini interni
prepare_img <- function(path, title) {
  img <- magick::image_read(path)
  img <- magick::image_trim(img) 
  grob <- grid::rasterGrob(img, interpolate = TRUE)
  
  wrap_elements(grob) + 
    ggtitle(title) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 5)),
      # Riduciamo a zero i margini del singolo pannello
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
}

# 3. Creazione dei pannelli (senza titoli per massimizzare lo spazio)
p_violin7 <- prepare_img(path7_v, "")
p_violin2 <- prepare_img(path2_v, "")
p_map7    <- prepare_img(path7_m, "")
p_map2    <- prepare_img(path2_m, "")

# 4. Assemblaggio 2x2 "stretto"
final_2x2 <- (p_violin7 + p_violin2) / (p_map7 + p_map2) + 
  # Regoliamo le altezze: 2.5 per i violini (dettagliati) e 1 per le mappe (piatte)
  plot_layout(heights = c(2., 1)) & 
  theme(
    # Riduce lo spazio tra le righe (valore negativo per "avvicinare")
    plot.spacing = unit(-1, "cm"), 
    # Ingrandisce le etichette A, B, C, D
    plot.tag = element_text(size = 24, face = "bold"),
    # Assicura che non ci siano margini extra nei subplot
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

# Aggiunta tag_levels separatamente per evitare conflitti di tema
final_2x2 <- final_2x2 + plot_annotation(tag_levels = 'A')
 
# 5. Salvataggio con altezza ridotta
# Ridurre height da 12 a 9 o 10 è il modo più efficace per eliminare il bianco sopra/sotto
ggsave("./figures/Combined_2x2_Tight.tiff", 
       final_2x2, width = 12, height = 9, dpi = 600, compression = "lzw")

print(final_2x2)
# 
# # 1. Define paths
# path7 <- "./figures/violins/Microstate_z_7_Plot_across_phases.tiff"
# path2 <- "./figures/violins/Microstate_Z_2_Plot_across_phases.tiff"
# map7  <- "./figures/brains/microstate_map_07.png"
# map2  <- "./figures/brains/microstate_map_02.png"
# 
# # 2. Read images
# img7 <- magick::image_read(path7)
# img2 <- magick::image_read(path2)
# 
# # 3. Convert to graphical objects
# g7 <- grid::rasterGrob(img7)
# g2 <- grid::rasterGrob(img2)
# 
# # 4. Create the patchwork with sub-titles
# # We use wrap_elements and add a theme-based title to each
# p1 <- wrap_elements(g7) + 
#   ggtitle("Microstate map 7") + 
#   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
# 
# p2 <- wrap_elements(g2) + 
#   ggtitle("Microstate map 2") + 
#   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
# 
# # 5. Combine and annotate
# final_plot <- p1 + p2 + 
#   plot_layout(ncol = 2) +
#   plot_annotation(
#     title = "",
#     theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
#   )
# 
# # 6. Display
# print(final_plot)
# 
# # Optional: Save with high resolution for your paper
# ggsave("./figures/Combined_Maps_violin_7_2.png", 
#         final_plot, width = 12, height = 6, dpi = 600)