# ==============================================================================
# SCRIPT 05: BEHAVIOURAL STATISTICS AND VISUALISATION
# Description: This script processes the in vivo phenotypic data, calculating 
# LC50s, EC50s, survival models, and generating publication-ready figures.
# ==============================================================================

# --- 0. SETUP: LOAD REQUIRED PACKAGES AND DEFINE THEMES ---

# Uncomment and run the following line if you need to install missing packages:
# install.packages(c("tidyverse", "survival", "survminer", "emmeans", "drc", "brglm2", "ARTool", "effectsize", "ggtext", "ggrepel", "multcomp"))

library(tidyverse)
library(survival)
library(survminer)
library(emmeans)
library(drc)
library(brglm2)
library(ARTool)
library(effectsize)
library(ggtext)
library(ggrepel)
library(multcomp)

# Global custom aesthetics for plots
sex_colors <- c("Male" = "#2E86AB", "Female" = "#A23B72")

biostatsquid_theme <- theme(
  plot.title = element_text(size = rel(1.4), face = "bold", hjust = 0.5),
  plot.title.position = "plot",
  plot.background = element_rect(fill = NA, colour = 'white'),
  panel.background = element_rect(fill = 'white'),
  panel.grid.major.y = element_line(colour = 'lightgray'),
  panel.grid.minor.y = element_line(colour = 'lightgray'),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  axis.line = element_line(colour = 'black', linewidth = 0.8),
  axis.ticks = element_line(colour = 'black', linewidth = 0.8),
  axis.text = element_text(colour = "black", face = "bold", size = rel(1.1)),
  axis.title = element_text(size = rel(1.2), face = "bold"),
  legend.title = element_text(face = "bold", size = rel(1.2)),
  legend.text = element_text(size = rel(1.1)),
  legend.key = element_blank()
)

# ==============================================================================
# --- 1. SEX-SPECIFIC ACUTE TOXICITY ---
# ==============================================================================
print("================= 1. ACUTE TOXICITY (LC50) =================")

# Load Data
Phenanthrene <- read_csv("Data/Phenanthrene.csv")

# Fit the LL.2 model using drc
Phe.LL.2.1 <- drm(Dead/Exposed ~ Conc, Sex, weights = Exposed, data = Phenanthrene,
                  type = "binomial", fct = LL.2())

print(summary(Phe.LL.2.1))

# Extract LC50s
print("--- LC50 Estimates ---")
ED(Phe.LL.2.1, 50, interval = "delta")

# Compare LC50s between sexes
print("--- Formal Comparison of LC50s ---")
EDcomp(Phe.LL.2.1, c(50, 50))

# Plot Acute Toxicity Curves
tiff("Figure1_Acute_Toxicity.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
cb_palette <- c("#2E86AB", "#A23B72") # Matched to global sex colors
plot(Phe.LL.2.1, broken = TRUE, ylim = c(0, 1),  xlim = c(0,5000),
     legendPos = c(4, 0.95),   
     xlab = "Concentration (µg/L)", ylab = "Probability of Mortality",
     col = cb_palette, lty = c(1, 2), lwd = 3,                 
     main = "Phenanthrene Exposure Exhibits Sex-Specific Acute Toxicity",
     cex.axis = 1.2, cex.lab = 1.3, cex.main = 1.3, font.main = 2, font.axis = 2, adj = 0)

# Add LC50 lines
LC50_male <- 190.26  
LC50_female <- 111.27  
text(LC50_male, 0.55, labels = paste("Male:", round(LC50_male, 2)), pos = 4, col = cb_palette[1], cex = 1.2)
text(LC50_female, 0.45, labels = paste("Female:", round(LC50_female, 2)), pos = 2, col = cb_palette[2], cex = 1.2)
abline(v = LC50_male, lty = 1, col = cb_palette[1])  
abline(v = LC50_female, lty = 5, col = cb_palette[2])  
abline(h = 0.5, lty = 5, col = "grey50")  
dev.off()


# ==============================================================================
# --- 2. SEX-SPECIFIC FEEDING SUPPRESSION ---
# ==============================================================================
print("================= 2. FEEDING SUPPRESSION =================")

# Load Data
feeding_data <- read_csv("Data/feeding_rate_data_NEW.csv")

# 2A. ARTool Non-Parametric ANOVA
feeding_data_factors <- feeding_data %>%
  mutate(Sex = as.factor(Sex), Treatment = factor(Treatment, levels = c("Seawater", "DMSO_control", "PHE_10_ug_l", "PHE_50_ug_l", "PHE_100_ug_l")))

model_art <- art(Feeding_Rate ~ Sex * Treatment, data = feeding_data_factors)
print("--- ARTool ANOVA Results ---")
print(anova(model_art))

print("--- Partial Eta Squared ---")
print(eta_squared(model_art, partial = TRUE))

# 2B. Plot Enhanced Boxplot
boxplot_enhanced <- ggplot(feeding_data_factors, aes(x = Treatment, y = Feeding_Rate, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.7, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), aes(color = Sex), shape = 16) +
  scale_fill_manual(values = sex_colors) +
  scale_color_manual(values = sex_colors) +
  scale_x_discrete(labels = c("Seawater", "DMSO", "PHE 10", "PHE 50", "PHE 100")) +
  labs(
    title = "Effect of PHE on Feeding Rate in *Parhyale hawaiensis*",
    subtitle = expression(paste("Interaction: ", italic("F")[4, 110], " = 3.98, ", italic("p"), " = 0.005")),
    x = "Treatment Group (μg/L)", y = "Feeding Rate (mg / Parhyale / hour)"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "sans"), 
    plot.title = element_markdown(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 15),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.position = "top", legend.title = element_text(size = 13, face = "bold")
  )
ggsave("Figure2a_Feeding_Boxplot.tiff", plot = boxplot_enhanced, width = 8, height = 6, dpi = 300, compression = "lzw")

# 2C. Dose-Response Modelling (LL.3)
feeding_data_numeric <- feeding_data %>%
  mutate(Concentration = case_when(
    Treatment %in% c("Seawater", "DMSO_control") ~ 0,
    Treatment == "PHE_10_ug_l" ~ 10,
    Treatment == "PHE_50_ug_l" ~ 50,
    Treatment == "PHE_100_ug_l" ~ 100))

final_unified_model <- drm(Feeding_Rate ~ Concentration, curveid = Sex, data = feeding_data_numeric, fct = LL.3())

print("--- EC50 Estimates from Unified Model ---")
ED(final_unified_model, 50, interval = "delta")
print("--- Formal EC50 Comparison and Sensitivity Ratio CI ---")
EDcomp(final_unified_model, c(50, 50), operator = "/", interval = "delta")

# 2D. Enhanced Dose-Response Curve Plot
ec50_labels_final <- data.frame(
  Sex = c("Male", "Female"), Concentration = c(50.19, 19.15),
  Feeding_Rate = predict(final_unified_model, data.frame(Concentration = 0, Sex = "Male")) / 2,
  label = c("EC50 = 50.2", "EC50 = 19.1")
)

sensitivity_plot_enhanced <- ggplot(feeding_data_numeric, aes(x = Concentration, y = Feeding_Rate, color = Sex)) +
  geom_jitter(aes(shape = Sex), size = 3, alpha = 0.7) +
  geom_smooth(method = "drm", method.args = list(fct = LL.3()), se = FALSE, linewidth = 1.2) +
  geom_segment(data = ec50_labels_final, aes(x = Concentration, y = 0, xend = Concentration, yend = Feeding_Rate), linetype = "dashed", linewidth = 0.8) +
  geom_segment(data = ec50_labels_final, aes(x = 0, y = Feeding_Rate, xend = Concentration, yend = Feeding_Rate), linetype = "dashed", linewidth = 0.8) +
  geom_text_repel(data = ec50_labels_final, aes(label = label), size = 6, fontface = "bold", nudge_y = 0.05, segment.color = 'grey50') +
  scale_color_manual(values = sex_colors) +
  labs(
    title = "Dose-Response Curves for PHE by Sex",
    subtitle = "Females are more sensitive to PHE than males (Sensitivity Ratio = 2.62)",
    x = "PHE Concentration (μg/L)", y = "Feeding Rate (mg / Parhyale / hour)"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    text = element_text(family = "sans"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 15),
    axis.title = element_text(size = 15, face = "bold"),
    legend.position = "top"
  )
ggsave("Figure2b_Sensitivity_Curves.tiff", plot = sensitivity_plot_enhanced, width = 8, height = 6, dpi = 300, compression = "lzw")


# ==============================================================================
# --- 3. CHEMOSENSORY RESPONSE IMPAIRMENTS ---
# ==============================================================================
print("================= 3. CHEMOSENSORY IMPAIRMENTS =================")

# Load Data
Feed_chemo <- read_csv("Data/Feed_chemosensation_NEW.csv")
Feed_chemo <- Feed_chemo %>%
  mutate(
    Response_Time_sec = ifelse(Response_Time_sec > 1800, 1800, Response_Time_sec),
    Sex = as.factor(Sex),
    Treatment = factor(Treatment, levels = c("Seawater", "DMSO", "PHE_10", "PHE_50", "PHE_100")),
    Group = factor(paste0(Treatment, ".", substr(Sex, 1, 1))),
    SurvObj = Surv(Response_Time_sec, Response_status)
  )

# 3A. Kaplan-Meier Model & Log-Rank Test
km_fit_chemo <- survfit(SurvObj ~ Group, data = Feed_chemo)
print("--- Overall Log-Rank Test (Chemosensory) ---")
surv_pvalue(km_fit_chemo, data = Feed_chemo)

# 3B. Impairment Probability (Bias-Reduced GLM)
brglm_chemo <- glm(Impairment_status ~ Treatment * Sex, data = Feed_chemo, family = binomial("logit"), method = "brglmFit")
emm_chemo <- as.data.frame(emmeans(brglm_chemo, ~ Treatment * Sex, type = "response"))

# Clean columns for plot
names(emm_chemo)[names(emm_chemo) == "asymp.LCL"] <- "lower.CL"
names(emm_chemo)[names(emm_chemo) == "asymp.UCL"] <- "upper.CL"

emm_chemo <- emm_chemo %>%
  mutate(
    Treatment_Label = case_when(
      Treatment == "PHE_10" ~ "PHE 10", Treatment == "PHE_50" ~ "PHE 50", Treatment == "PHE_100" ~ "PHE 100",
      TRUE ~ as.character(Treatment)),
    GroupInteraction = fct_reorder(interaction(Treatment_Label, Sex, sep = ", "), as.numeric(interaction(Sex, Treatment)))
  )

dot_whisker_chemo <- ggplot(emm_chemo, aes(x = prob, y = GroupInteraction, colour = Sex)) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.3, linewidth = 0.8) +
  geom_point(size = 4) + scale_color_manual(values = sex_colors) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Predicted Probability of Impairment (%)", y = NULL,
    title = "PHE Exposure Increases Probability of Chemosensory Impairment",
    subtitle = "Impairment defined as response latency > 300 seconds"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12),
    legend.position = "top", axis.text = element_text(face = "bold"), axis.title = element_text(face = "bold")
  )
ggsave("Figure3_Chemo_Impairment.tiff", plot = dot_whisker_chemo, width = 8, height = 7, dpi = 300, compression = "lzw")


# ==============================================================================
# --- 4. REPRODUCTIVE LATENCY IMPAIRMENT ---
# ==============================================================================
print("================= 4. REPRODUCTIVE LATENCY =================")

# Load Data
precop_data <- read_csv("Data/Precop_data.csv")
precop_data <- precop_data %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Seawater", "DMSO", "PHE_10", "PHE_50", "PHE_100")),
    SurvObj = Surv(Response_Time_hours, Pairing_status)
  )

# 4A. Cox Proportional Hazards Model
cox_model <- coxph(SurvObj ~ Treatment, data = precop_data)
print("--- Cox Proportional Hazards Model Summary ---")
print(summary(cox_model))

print("--- Proportional Hazards Assumption Test ---")
print(cox.zph(cox_model))

# 4B. Kaplan-Meier Plot for Precopulatory Pairing
km_fit_precop <- survfit(SurvObj ~ Treatment, data = precop_data)

Kaplan_Plot_Precop <- ggsurvplot(
  km_fit_precop, data = precop_data, conf.int = FALSE, pval = TRUE, risk.table = TRUE,
  surv.median.line = "hv", legend.title = "Treatment", palette = "jco",
  xlab = "Time to Amplexus Formation (hours)", ylab = "Probability of Not Paired",
  break.time.by = 24, risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
  tables.theme = theme_cleantable(),
  ggtheme = biostatsquid_theme + theme(legend.position = "right")
)
# Save arranged KM plot
arranged_precop <- arrange_ggsurvplots(list(Kaplan_Plot_Precop), print = FALSE, ncol = 1, nrow = 2, heights = c(0.7, 0.3))
ggsave("Figure4a_Precop_KM.tiff", plot = arranged_precop, width = 11, height = 8.5, dpi = 300, compression = "lzw")

# 4C. Impairment EC50 (Dose-Response for reproductive failure)
bmd_data <- precop_data %>%
  group_by(Treatment) %>%
  summarise(impaired = sum(Impairment_status == 1), total = n(), .groups = 'drop') %>%
  mutate(dose = case_when(
    Treatment %in% c("Seawater", "DMSO") ~ 0, Treatment == "PHE_10" ~ 10, Treatment == "PHE_50" ~ 50, Treatment == "PHE_100" ~ 100))

bmd_model_3p <- drm(impaired / total ~ dose, weights = total, data = bmd_data, fct = L.3(names=c("Slope", "Upper", "ED50")), type = "binomial")

print("--- EC50 Estimation (Reproductive Impairment) ---")
ED(bmd_model_3p, 50, type = "absolute", interval = "delta")

# 4D. Bias-Reduced GLM for Impairment Visualisation
impairment_model_precop <- glm(Impairment_status ~ Treatment, data = precop_data, family = binomial("logit"), method = "brglmFit")
pred_probs_precop <- as.data.frame(emmeans(impairment_model_precop, ~ Treatment, type = "response"))

dot_whisker_precop <- ggplot(pred_probs_precop, aes(x = prob, y = fct_rev(Treatment))) +
  geom_errorbarh(aes(xmin = asymp.LCL, xmax = asymp.UCL), height = 0.3, linewidth = 0.8, colour = "#2E86AB") +
  geom_point(size = 4, color = "#A23B72") +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Predicted Probability of Impairment (%)", y = NULL,
    title = "PHE Exposure Increases Probability of Reproductive Impairment",
    subtitle = "Impairment defined as pairing latency > 48 hours"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12),
    axis.text = element_text(face = "bold"), axis.title = element_text(face = "bold")
  )
ggsave("Figure4b_Precop_Impairment.tiff", plot = dot_whisker_precop, width = 8, height = 6, dpi = 300, compression = "lzw")

print("================= SCRIPT COMPLETE =================")
