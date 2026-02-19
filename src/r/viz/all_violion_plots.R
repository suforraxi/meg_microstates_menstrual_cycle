library(ggplot2)

# Caricamento dataset
data <- read.csv("./data/merged_data_microstates_for_R.csv")

data$Session <- as.factor(data$Session)
data$sub <- as.factor(data$sub)

# --- 1. SETTINGS ---
plot_vars <- list(
  "Microstate_z_7"                = "Microstate z-score occurrences",
  "Microstate_z_2"                = "Microstate z-score occurrences",
  "Hormones_PC1_within"          = "PC1 loads of Hormones",
  "WELLBEING"                    = "Total Wellbeing Score",
  "Autonomy"                     = "Autonomy Score",
  "EnvironmentalMastery"         = "Environmental Mastery Score",
  "PersonalGrowth"               = "Personal Growth Score",
  "PositiveRelationswithOthers"  = "Positive Relations Score",
  "PurposeinLife"                = "Purpose in Life Score",
  "Self.Acceptance"              = "Self-Acceptance Score"
)

out_dir <- "./figures/violins/"

# --- 2. LOOP ---
for (var_name in names(plot_vars)) {
  
  p <- ggplot(data, aes(x = Session, y = .data[[var_name]], fill = Session)) +
    geom_violin(alpha = 0.3, trim = FALSE) +
    geom_line(aes(group = sub), alpha = 0.2, color = "grey") +
    geom_point(alpha = 0.5) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    
    scale_x_discrete(labels = c("1" = "Early Follicular", 
                                "2" = "Peri-ovulatory", 
                                "3" = "Mid-luteal")) +
    
    theme_classic() +
    labs(
      x = "Phase", 
      y = plot_vars[[var_name]],
      fill = "Phase", # Titolo per la legenda
      title = ""
    ) +
    
    # --- MODIFICHE AL TEMA: GRASSETTO E DIMENSIONI ---
    theme(
      # Titoli degli assi (X e Y)
      axis.title = element_text(size = 22, face = "bold", color = "black"),
      # Testo dei tick (le etichette sulle scale)
      axis.text = element_text(size = 18, face = "bold", color = "black"),
      # Se vuoi mostrare la legenda, usa "right" o "bottom", altrimenti "none"
      legend.position = "none",
      # Spessore delle linee degli assi
      axis.line = element_line(linewidth = 1)
    )
  
  # --- 3. SAVE ---
  clean_filename <- paste0(gsub("\\.", "_", var_name), "_Plot_across_phases.tiff")
  
  ggsave(filename = paste0(out_dir, clean_filename), 
         plot = p, 
         width = 10, height = 7, units = "in", dpi = 600, compression = "lzw")
  
  cat("Saved:", clean_filename, "\n")
}

# library(ggplot2)
# 
# # Load the dataset
# data <- read.csv("./data/merged_data_microstates_for_R.csv")
# 
# data$Session <- as.factor(data$Session)
# # Ensure 'sub' is treated as a factor (random effect grouping variable)
# data$sub <- as.factor(data$sub)
# 
# #
# data$Session <- as.factor(data$Session)
# 
# # --- 1. SETTINGS ---
# # Define your variables and their corresponding Y-axis labels
# plot_vars <- list(
#   "Microstate_z_7"               = "Microstate z-score occurrences",
#   "Microstate_z_2"               = "Microstate z-score occurrences",
#   "Hormones_PC1_within"          = "PC1 loads of Hormones",
#   "WELLBEING"                    = "Total Wellbeing Score",
#   "Autonomy"                     = "Autonomy Score",
#   "EnvironmentalMastery"         = "Environmental Mastery Score",
#   "PersonalGrowth"               = "Personal Growth Score",
#   "PositiveRelationswithOthers"  = "Positive Relations Score",
#   "PurposeinLife"                = "Purpose in Life Score",
#   "Self.Acceptance"              = "Self-Acceptance Score"
# )
# 
# # Output directory
# out_dir <- "./figures/violins/"
# 
# # --- 2. LOOP ---
# for (var_name in names(plot_vars)) {
#   
#   p <- ggplot(data, aes(x = Session, y = .data[[var_name]], fill = Session)) +
#     # Violin and connecting lines
#     geom_violin(alpha = 0.3, trim = FALSE) +
#     geom_line(aes(group = sub), alpha = 0.2, color = "grey") +
#     geom_point(alpha = 0.5) +
#     geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
#     
#     # Standard Discrete Labels
#     scale_x_discrete(labels = c("1" = "Early Follicular", 
#                                 "2" = "Peri-ovulatory", 
#                                 "3" = "Mid-luteal")) +
#     
#     # Themes and Dynamic Labels
#     theme_classic() +
#     labs(
#       x = "Phase", 
#       y = plot_vars[[var_name]],
#       title = "" #paste("Distribution of", var_name, "across Phases")
#     ) +
#     theme(legend.position = "none")
# 
#   # --- 3. SAVE ---
#   # Generate clean filename (replacing dots with underscores for compatibility)
#   clean_filename <- paste0(gsub("\\.", "_", var_name), "_Plot_across_phases.tiff")
#   
#   ggsave(filename = paste0(out_dir, clean_filename), 
#          plot = p, 
#          width = 8, height = 6, units = "in", dpi = 600, compression = "lzw")
#   
#   cat("Saved:", clean_filename, "\n")
# }