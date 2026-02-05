library(ggplot2)
library(readr)
library(dplyr)

# 1. Configurazione percorsi e parametri
folder_res <- './data/'
agg_occurrences_f <- file.path(folder_res, 'aggregated_microstate_occurrences.csv')
mapOI <- 7

data <- read_csv(agg_occurrences_f)

# 1. Prepare data
data_mapOI <- data %>%
  filter(Microstate == mapOI) %>%
  mutate(
    sub = as.factor(sub),
    Session = as.factor(Session)
  )

# 2. Create the Bar Plot with Alpha and Standard Colors
p_alpha <- ggplot(data_mapOI, aes(x = sub, y = Occurrences, fill = Session)) +
  # Add alpha here (e.g., 0.6 for significant transparency)
  geom_bar(
    stat = "identity", 
    position = "stack",#position_dodge(width = 0.8), 
    width = 0.7,
    alpha = 0.6,      # <--- Change this value to adjust transparency
    #color = "black",  # Adding a black outline helps define bars when alpha is low
    #size = 0.05
  ) +
  
  scale_fill_discrete(
    name = "Phase", 
    labels = c("01" = "Early Follicular", 
               "02" = "Peri-ovulatory", 
               "03" = "Mid-luteal")
  ) +
  
  theme_minimal() +
  labs(
    title = "", #paste("Microstate", mapOI, "Occurrences by Subject"),
    x = "Subject",
    y = "Occurrences"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

print(p_alpha)

# 3. Export
#ggsave("micro.png", plot = p_alpha, width = 11, height = 5, dpi = 300)
ggsave("./figures/BAR_Microstate_map_7_across_phases.tiff", plot = p_alpha, 
       width = 8, height = 6, units = "in", dpi = 600, compression = "lzw")
