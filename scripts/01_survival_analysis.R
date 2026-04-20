# ============================================================
# 01_survival_analysis.R
# TCGA-PAAD Survival Analysis
# Author: Simran Randhawa
# ============================================================

library(survival)
library(ggplot2)
library(tidyverse)
library(ggfortify)

# ── 1. Load clean data ───────────────────────────────────────
df <- read.csv("data/tcga_paad_clean.csv", stringsAsFactors = FALSE)
df <- df %>% filter(!is.na(survival_months) & !is.na(vital_status_binary))
cat("Patients loaded:", nrow(df), "\n")

# ── 2. Overall survival curve ────────────────────────────────
fit_overall <- survfit(Surv(survival_months, vital_status_binary) ~ 1, data = df)

png("results/figures/03_overall_survival.png", width=800, height=600, res=120)
plot(fit_overall,
     xlab = "Time (months)", ylab = "Survival Probability",
     main = "Overall Survival — TCGA-PAAD (n=185)",
     col = "steelblue", lwd = 2, conf.int = TRUE)
abline(v = summary(fit_overall)$table["median"], col="red", lty=2)
legend("topright", legend=paste("Median:", round(summary(fit_overall)$table["median"],1), "months"),
       col="red", lty=2)
dev.off()

cat("Median survival:", summary(fit_overall)$table["median"], "months\n")

# ── 3. Survival by stage ─────────────────────────────────────
df_stage <- df %>%
  mutate(stage_simple = case_when(
    str_detect(diagnoses.ajcc_pathologic_stage, "Stage I") & 
      !str_detect(diagnoses.ajcc_pathologic_stage, "Stage II|Stage III|Stage IV") ~ "Stage I",
    str_detect(diagnoses.ajcc_pathologic_stage, "Stage II") ~ "Stage II",
    str_detect(diagnoses.ajcc_pathologic_stage, "Stage III") ~ "Stage III",
    str_detect(diagnoses.ajcc_pathologic_stage, "Stage IV") ~ "Stage IV",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(stage_simple))

fit_stage <- survfit(Surv(survival_months, vital_status_binary) ~ stage_simple, data = df_stage)

png("results/figures/04_survival_by_stage.png", width=900, height=650, res=120)
plot(fit_stage,
     xlab = "Time (months)", ylab = "Survival Probability",
     main = "Survival by AJCC Stage — TCGA-PAAD",
     col = c("#E64B35","#4DBBD5","#00A087","#3C5488"), lwd = 2)
legend("topright", legend = c("Stage I","Stage II","Stage III","Stage IV"),
       col = c("#E64B35","#4DBBD5","#00A087","#3C5488"), lwd = 2)
dev.off()

# Log-rank test for stage
lr_stage <- survdiff(Surv(survival_months, vital_status_binary) ~ stage_simple, data = df_stage)
cat("Stage log-rank p-value:", 1 - pchisq(lr_stage$chisq, df=length(lr_stage$n)-1), "\n")

# ── 4. Survival by gender ────────────────────────────────────
fit_gender <- survfit(Surv(survival_months, vital_status_binary) ~ demographic.gender, data = df)

png("results/figures/05_survival_by_gender.png", width=900, height=650, res=120)
plot(fit_gender,
     xlab = "Time (months)", ylab = "Survival Probability",
     main = "Survival by Gender — TCGA-PAAD",
     col = c("#E64B35","#4DBBD5"), lwd = 2)
legend("topright", legend = c("Female","Male"),
       col = c("#E64B35","#4DBBD5"), lwd = 2)
dev.off()

cat("Script 01 complete. Figures saved.\n")