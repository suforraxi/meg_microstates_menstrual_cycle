library(lmerTest)
library(dplyr)

# --- 1. Load and Prep Data ---
data <- read.csv("./merged_data_microstates_for_R.csv")
data$Session <- as.factor(data$Session)
data$sub <- as.factor(data$sub)

# --- 2. Enhanced Function for ANOVA & Tukey Extraction ---
process_anova_results <- function(var_name, df) {
  
  # A. Fit Repeated Measures ANOVA
  formula_rm <- as.formula(paste(var_name, "~ Session + Error(sub/Session)"))
  anova_obj <- aov(formula_rm, data = df)
  anova_summ <- summary(anova_obj)
  
  # B. Extract ANOVA Main Effect (from the sub:Session stratum)
  # This stratum contains the effect of the menstrual cycle phase
  anova_main <- as.data.frame(anova_summ[[2]][[1]]) 
  
  anova_table <- data.frame(
    Variable      = var_name,
    Effect        = c("Session", "Residuals"),
    Df            = anova_main$Df,
    Sum_Sq        = round(anova_main$`Sum Sq`, 3),
    Mean_Sq       = round(anova_main$`Mean Sq`, 3),
    F_value       = round(anova_main$`F value`, 3),
    P_value       = round(anova_main$`Pr(>F)`, 3)
  )
  
  # C. Tukey Post-hoc
  formula_one_way <- as.formula(paste(var_name, "~ Session"))
  anova_one_way <- aov(formula_one_way, data = df)
  tukey_res <- TukeyHSD(anova_one_way)
  
  tukey_df <- as.data.frame(tukey_res$Session)
  tukey_df$Comparison <- rownames(tukey_df)
  tukey_df$Variable <- var_name
  
  # D. Replace Session Numbers with Clinical Names
  tukey_df$Comparison <- tukey_df$Comparison %>%
    gsub("1", "Early Follicular", .) %>%
    gsub("2", "Peri-ovulatory", .) %>%
    gsub("3", "Mid-luteal", .)
  
  # E. Formatting Tukey
  tukey_df <- tukey_df %>%
    mutate(across(c(diff, lwr, upr, `p adj`), ~ round(., 3))) %>%
    mutate(Significance = cut(`p adj`, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                              labels = c("***", "**", "*", "ns"))) %>%
    select(Variable, Comparison, diff, lwr, upr, `p adj`, Significance)
  
  return(list(anova = anova_table, tukey = tukey_df))
}

# --- 3. Execute ---
target_vars <- c("Hormones_PC1_within", "WELLBEING", "Microstate_z_7", "Microstate_z_2",
                 "PersonalGrowth", "Autonomy", "EnvironmentalMastery", 
                 "PositiveRelationswithOthers", "PurposeinLife", "Self.Acceptance")

all_processed <- lapply(target_vars, process_anova_results, df = data)

# --- 4. Combine and Export ---

# Table 1: Main Effects (F, Df, MS)
master_anova_table <- do.call(rbind, lapply(all_processed, function(x) x$anova))

# Table 2: Post-Hoc (Tukey)
master_tukey_table <- do.call(rbind, lapply(all_processed, function(x) x$tukey))

# Save files
write.csv(master_anova_table, 
          "./tables/Main_ANOVA_Effects.csv", 
          row.names = FALSE)

write.csv(master_tukey_table, 
          "./tables/PostHoc_Results_Summary_3digits.csv", 
          row.names = FALSE)

cat("Success! Two tables saved: 'Main_ANOVA_Effects.csv' and 'PostHoc_Results_Summary_3digits.csv'")