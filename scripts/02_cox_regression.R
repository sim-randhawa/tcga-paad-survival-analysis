# ============================================================
# 02_cox_regression.R
# Cox Proportional Hazards Model — TCGA-PAAD
# Author: Simran Randhawa
# ============================================================

library(survival)
library(ggplot2)
library(tidyverse)

# ── 1. Load and prepare data ─────────────────────────────────
df <- read.csv("data/tcga_paad_clean.csv", stringsAsFactors = FALSE)
df <- df %>% filter(!is.na(survival_months) & !is.na(vital_status_binary))

# Create stage variable
df <- df %>%
  mutate(stage_simple = case_when(
    str_detect(diagnoses.ajcc_pathologic_stage, "Stage I") &
      !str_detect(diagnoses.ajcc_pathologic_stage, "Stage II|Stage III|Stage IV") ~ "Stage I",
    str_detect(diagnoses.ajcc_pathologic_stage, "Stage II") ~ "Stage II",
    str_detect(diagnoses.ajcc_pathologic_stage, "Stage III") ~ "Stage III",
    str_detect(diagnoses.ajcc_pathologic_stage, "Stage IV") ~ "Stage IV",
    TRUE ~ NA_character_
  ))

# Convert age to years
df$age_years <- as.numeric(df$demographic.age_at_index)

# ── 2. Fit Cox model ─────────────────────────────────────────
cox_df <- df %>%
  filter(!is.na(stage_simple) & !is.na(age_years) & !is.na(demographic.gender)) %>%
  mutate(stage_simple = factor(stage_simple, levels = c("Stage I","Stage II","Stage III","Stage IV")),
         gender = factor(demographic.gender))

cox_model <- coxph(Surv(survival_months, vital_status_binary) ~
                     stage_simple + age_years + gender,
                   data = cox_df)

cat("Cox Model Summary:\n")
print(summary(cox_model))

# ── 3. Forest plot of hazard ratios ──────────────────────────
cox_results <- as.data.frame(summary(cox_model)$conf.int)
cox_results$variable <- rownames(cox_results)
colnames(cox_results) <- c("HR", "exp_neg_coef", "lower_95", "upper_95", "variable")

cox_results <- cox_results %>%
  mutate(variable = recode(variable,
                           "stage_simpleStage II" = "Stage II vs I",
                           "stage_simpleStage III" = "Stage III vs I",
                           "stage_simpleStage IV" = "Stage IV vs I",
                           "age_years" = "Age (per year)",
                           "gendermale" = "Gender (Male vs Female)"
  ))

png("results/figures/06_cox_forest_plot.png", width=900, height=600, res=120)
ggplot(cox_results, aes(x = HR, y = reorder(variable, HR))) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), height = 0.2, color = "steelblue") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "Cox Proportional Hazards Model — TCGA-PAAD",
       subtitle = "Hazard Ratios with 95% Confidence Intervals",
       x = "Hazard Ratio (log scale)", y = "") +
  theme_classic(base_size = 12)
dev.off()

cat("Script 02 complete.\n")