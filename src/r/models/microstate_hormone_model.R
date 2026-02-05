# --- 0. LIBRARIES ---
library(lmerTest)
library(performance)

# --- 1. DATA SETUP ---
# Standardizing paths and factor types
data_path <- "./merged_data_microstates_for_R.csv"
df <- read.csv(data_path)

df$sub     <- as.factor(df$sub)
df$Session <- as.factor(df$Session)

# --- 2. STEP 1: MODEL SELECTION (ML) ---
# We use REML = FALSE to compare models with different fixed effects
m_null <- lmer(Microstate_z_7 ~ 1 + (1 | sub), data = df, REML = FALSE)
m_full <- lmer(Microstate_z_7 ~ Hormones_PC1_within + (1 | sub), data = df, REML = FALSE)

# Likelihood Ratio Test
model_comp <- anova(m_null, m_full)

# --- 3. STEP 2: FINAL ESTIMATION (REML) ---
# Refit the full model with REML = TRUE for unbiased parameter estimates
m_final <- lmer(Microstate_z_7 ~ Hormones_PC1_within + (1 | sub), data = df, REML = TRUE)

# Extract Coefficients and R-Squared
coeffs  <- as.data.frame(summary(m_final)$coefficients)
conf_int <- as.data.frame(confint(m_final, method = "Wald")) # 95% CI
r2_vals <- r2_nakagawa(m_final)

# --- 4. DATA EXTRACTION & FORMATTING ---
# Target the specific predictor row
res <- coeffs["Hormones_PC1_within", ]

# Store clean variables
beta_val <- round(res$Estimate, 3)
t_val    <- round(res$`t value`, 2)
df_val   <- round(res$df, 1)
p_val    <- res$`Pr(>|t|)`
ci_low   <- round(conf_int["Hormones_PC1_within", 1], 3)
ci_high  <- round(conf_int["Hormones_PC1_within", 2], 3)

# --- 5. CLEAN REPORTING ---
cat("\n================================================")
cat("\nSUMMARY OF BIOLOGICAL LINK (Hormones -> MS7)")
cat("\n================================================")
cat("\nModel Comparison (LRT): Chi-sq =", round(model_comp$Chisq[2], 3), 
    "| p =", round(model_comp$`Pr(>Chisq)`[2], 4))
cat("\nEffect Size (Marginal R2):", round(r2_vals$R2_marginal, 3))
cat("\n------------------------------------------------")
cat("\nBeta Estimate: ", beta_val, "[95% CI:", ci_low, ",", ci_high, "]")
cat("\nStatistics:    t(", df_val, ") =", t_val)
cat("\nP-value:       ", ifelse(p_val < 0.001, "< .001", round(p_val, 4)))
cat("\n------------------------------------------------\n")

# Formatted Sentence for Manuscript
manuscript_text <- sprintf(
  "Hormonal PC1 significantly predicted Microstate 7 occurrences (beta = %s, 95%% CI [%s, %s], t(%s) = %s, p %s).",
  beta_val, ci_low, ci_high, df_val, t_val, ifelse(p_val < 0.001, "< .001", paste("=", round(p_val, 4)))
)

cat("MANUSCRIPT TEXT:\n", manuscript_text, "\n")