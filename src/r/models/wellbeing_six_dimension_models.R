library(lmerTest)
library(dplyr)
library(performance)
library(MuMIn) 
library(ggplot2)

# --- 1. DATA PREP & MERGING ---
df <- read.csv("./data/data_microstate_7_for_R.csv")
df_demo <- read.csv("./data/demographics.csv")

# Standardize and merge
df$sub <- as.character(df$sub)
df_demo$sub <- as.character(df_demo$sub)
df <- df %>% left_join(df_demo %>% select(sub, Age, Edu), by = "sub")
df$sub <- as.factor(df$sub)

wellbeing_subscales <- c("PersonalGrowth", "Autonomy", "EnvironmentalMastery", 
                         "PositiveRelationswithOthers", "PurposeinLife", "Self.Acceptance")
all_targets <- c("WELLBEING", wellbeing_subscales)

# --- 2. THE BIOLOGICAL LINK ---
# Investigating if the already computed Hormonal PC1 predicts Microstate activity
cat("\n--- STEP 1: BIOLOGICAL LINK (Hormones -> Microstates) ---\n")
m_bio_link <- lmer(Microstate_z ~ Hormones_PC1_within + (1 | sub), data = df)
print(summary(m_bio_link))

# --- 3. WELLBEING PERFORMANCE & FDR ---
cat("\n--- STEP 2: MODEL PERFORMANCE (AIC/BIC/LRT) ---\n")
performance_results <- data.frame()

for (target in all_targets) {
  f_base <- as.formula(paste(target, "~ Age + Edu + (1 | sub)"))
  f_full <- as.formula(paste(target, "~ Age + Edu + Hormones_PC1_within + Microstate_z + (1 | sub)"))
  
  m_base <- lmer(f_base, data = df, REML = FALSE)
  m_full <- lmer(f_full, data = df, REML = FALSE)
  
  lrt <- anova(m_base, m_full)
  r2_vals <- r2_nakagawa(m_full)
  
  performance_results <- rbind(performance_results, data.frame(
    Target = target,
    AIC_Base = AIC(m_base), AIC_Full = AIC(m_full),
    BIC_Base = BIC(m_base), BIC_Full = BIC(m_full), # Added BIC
    LRT_Chisq = lrt[2, "Chisq"],
    P_LRT_Uncorrected = lrt[2, "Pr(>Chisq)"],
    Marginal_R2 = r2_vals$R2_marginal,
    Conditional_R2 = r2_vals$R2_conditional
  ))
}

# FDR Correct by N=6 (Excluding Global WELLBEING)
sub_lrt <- performance_results %>% filter(Target != "WELLBEING") %>%
  mutate(P_LRT_FDR = p.adjust(P_LRT_Uncorrected, method = "fdr"))
global_lrt <- performance_results %>% filter(Target == "WELLBEING") %>% mutate(P_LRT_FDR = NA)

# Combine and print with 2 decimal places for scannability
final_table <- rbind(global_lrt, sub_lrt)
print(final_table %>% mutate(across(where(is.numeric), ~ round(., 2))))

# --- 4. PREDICTIVE VALIDATION (LOOCV) ---
loocv_lmer <- function(formula, data, out_dir = "./figures/") {
  # Extract the dependent variable name
  dependent_var <- all.vars(formula)[1]
  
  # Initialize results vectors
  predictions <- numeric(nrow(data)); residuals <- numeric(nrow(data))
  fitted_values <- numeric(nrow(data)); real_values <- numeric(nrow(data))
  r2_marginal <- numeric(nrow(data)); r2_conditional <- numeric(nrow(data))
  
  # LOOCV Loop
  for (i in 1:nrow(data)) {1
    train_data <- data[-i, ]; test_data <- data[i, ]
    model <- lmer(formula, data = train_data)
    
    predictions[i] <- predict(model, newdata = test_data, allow.new.levels = TRUE)
    residuals[i] <- test_data[[dependent_var]] - predictions[i]
    fitted_values[i] <- predictions[i]
    real_values[i] <- test_data[[dependent_var]]
    
    r2 <- r.squaredGLMM(model)
    r2_marginal[i] <- r2[1]; r2_conditional[i] <- r2[2]
  }
  
  mse <- mean((data[[dependent_var]] - predictions)^2)
  
  # Prep data for plotting
  results_df <- data.frame(
    Observed = real_values,
    Predicted = fitted_values,
    Residuals = residuals,
    Phase = factor(data$Session, 
                   levels = c("1", "2", "3"),
                   labels = c("Early Follicular", "Peri-ovulatory", "Mid-luteal"))
  )
  
  phase_colors <- c("Early Follicular" = "red", "Peri-ovulatory" = "green", "Mid-luteal" = "blue")
  
  # --- Plot 1: Residuals ---
  p1 <- ggplot(results_df, aes(x = Predicted, y = Residuals, color = Phase)) +
    geom_point(alpha = 0.3, size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = phase_colors, name = "Phase") +
    facet_wrap(~Phase) +
    labs(title = " ", #paste("LOOCV Residuals:", dependent_var),
         x = "Predicted Value", y = "Residual (Observed - Predicted)") +
    theme_minimal() + theme(legend.position = "bottom")
  
  # --- Plot 2: Observed vs Predicted ---
  p2 <- ggplot(results_df, aes(x = Predicted, y = Observed, color = Phase)) +
    geom_point(alpha = 0.3, size = 2.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = phase_colors, name = "Phase") +
    facet_wrap(~Phase) +
    labs(title = "",#paste("LOOCV Observed vs Predicted:", dependent_var), 
         x = "Predicted Value", y = "Observed Value") +
    theme_minimal() + theme(legend.position = "bottom")
  
  # --- SAVE FIGURES ---
  file_p1 <- paste0(out_dir, "LOOCV_Residuals_", dependent_var, ".tiff")
  file_p2 <- paste0(out_dir, "LOOCV_Observed_vs_Predicted_", dependent_var, ".tiff")
  
  ggsave(file_p1, plot = p1, width = 10, height = 6, dpi = 600, compression = "lzw")
  ggsave(file_p2, plot = p2, width = 10, height = 6, dpi = 600, compression = "lzw")
  
  print(p1); print(p2)
  cat("Figures saved to:", out_dir, "\n")
  
  return(list(predictions = predictions, mse = mse, r2_marginal = r2_marginal, r2_conditional = r2_conditional))
}


cat("\n--- STEP 3: CROSS-VALIDATION (Personal Growth) ---\n")
cv_results <- loocv_lmer(PersonalGrowth ~ Age + Edu + Microstate_z + Hormones_PC1_within + (1 | sub), df)

# Print the Mean Squared Error
cat("Mean Squared Error:", cv_results$mse, "\n")

# Print the average R² values
cat("Average Marginal R²:", mean(cv_results$r2_marginal), "\n")
cat("Average Conditional R²:", mean(cv_results$r2_conditional), "\n")


