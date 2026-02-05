library(dplyr)

# --- 1. Load the demographic data ---
path_demo <- "./data/demographics.csv"
df_demo <- read.csv(path_demo)

# --- 2. Compute Mean and Standard Deviation ---
demo_summary <- df_demo %>%
  summarise(
    Mean_Age = mean(Age, na.rm = TRUE),
    SD_Age   = sd(Age, na.rm = TRUE),
    Mean_Edu = mean(Edu, na.rm = TRUE),
    SD_Edu   = sd(Edu, na.rm = TRUE)
  )

# --- 3. Print the results formatted to 2 decimal places ---
print(demo_summary %>% mutate(across(everything(), ~round(., 2))))