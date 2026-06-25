# Loading the database

germinacion <- read.csv("~/InternshipAlternaria/germinacion.csv")

View(germinacion)

#Cleaning the N/A's and adding the Total Spores count

raw_data <- read.csv("germinacion.csv", stringsAsFactors = FALSE)

cleaned_data <- raw_data %>%
  filter(!is.na(germinated) & !is.na(no.germinated) & treatment != "") %>%
  mutate(
    germinated = as.numeric(germinated),
    no_germinated = as.numeric(no.germinated),
    total_spores = germinated + no.germinated,
    pct_germination = (germinated / total_spores) * 100
  )

## Percentage of Inhibition of Germination (PIG)
```
Based on Kos K. 2002 

% of germination (%G) = (germinated / germinated + no.germinated) x 100

PIG = (Positive Control %G - Treatment %G / Positive Control %G ) x 100

´´´
# Calculate the control PIG for each set and integrate them into the data base

control_rates <- cleaned_data %>%
  filter(treatment == "C+") %>%
  group_by(inoculation, treat.set) %>%
  summarise(control_pct_mean = mean(pct_germination, na.rm = TRUE), .groups = 'drop')

analyzed_data <- cleaned_data %>%
  left_join(control_rates, by = c("inoculation", "treat.set")) %>%
 
# Calculating the PIG of each treatment
  # Excluding the positive control itself from downstream fungicide efficiency testing

    mutate(
    PIG = ((control_pct_mean - pct_germination) / control_pct_mean) * 100
  ) %>%
   filter(treatment != "C+")


anova_model <- lm(PIG ~ treatment + inoculation + treat.set, 
                  data = analyzed_data, 
                  weights = total_spores)

par(mfrow = c(2, 2))
plot(anova_model)


print("--- POST-HOC PAIRWISE COMPARISONS (TUKEY) ---")
post_hoc <- emmeans(anova_model, specs = pairwise ~ treatment, adjust = "tukey")
print(post_hoc$contrasts)


summary_table <- analyzed_data %>%
  group_by(treatment) %>%
  summarise(
    Mean_PIG = mean(PIG, na.rm = TRUE),
    SD_PIG = sd(PIG, na.rm = TRUE),
    Total_Spores_Evaluated = sum(total_spores)
  ) %>%
  arrange(desc(Mean_PIG))

print("--- COMPOUND EFFICIENCY SUMMARY RANKING ---")
print(summary_table)



##graphs

paired_subs <- c("aa", "ac", "al", "at", "am", "aasc", "as", "af")

dose_data <- analyzed_data %>%
  filter(gsub("0\\.5$|5$", "", treatment) %in% paired_subs) %>%
  mutate(
    Concentration = ifelse(grepl("0\\.5$", treatment), "0.5 % ", " 5 %"),
    Substance = gsub("0\\.5$|5$", "", treatment)
  )

othersubs <- analyzed_data %>%
  filter(!treatment %in% dose_data$treatment)

print("--- Statistical Analysis for First Group ---")

dose_model <- lm(PIG ~ Substance * Concentration + inoculation, data = dose_data, weights = total_spores)
print(Anova(dose_model, type = "III"))




print(graph_1)

graph_2 <- ggplot(othersubs, aes(x = reorder(treatment, -PIG), y = PIG)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.shape = 16) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "darkgreen", size = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey40") +
  labs(title = "Fungicidal Efficiency: Standalone Compounds & Controls",
       x = "Treatment / Compound Code", y = "Percent Inhibition of Germination (PIG %)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(graph_2)









```{r}
#| message: false
#| warning: false

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest) # Ensures p-values and denominator df are calculated
library(emmeans)

# =========================================================================
# 1. CLEANING & RECODING DATA
# =========================================================================

cleaned_data <- raw_data %>%
  filter(!is.na(germinated) & !is.na(no.germinated) & treatment != "") %>%
  mutate(
    germinated = as.numeric(germinated),
    no.germinated = as.numeric(no.germinated),
    total_spores = germinated + no.germinated, 
    pct_germination = (germinated / total_spores) * 100 
  )

# Calculate Positive Control Means per run set
control_rates <- cleaned_data %>%
  filter(treatment == "C+") %>% 
  group_by(inoculation, treat.set) %>%
  summarise(control_pct_mean = mean(pct_germination, na.rm = TRUE), .groups = 'drop')

# Calculate PIG and assign clean groups
processed_data <- cleaned_data %>%
  left_join(control_rates, by = c("inoculation", "treat.set")) %>% 
  mutate(
    PIG = ((control_pct_mean - pct_germination) / control_pct_mean) * 100
  ) %>%
  mutate(
    substance_clean = case_when(
      treatment %in% c("ac0.5", "ac5") ~ "Acetic acid", # <-- FIXED TYPO HERE
      treatment %in% c("al0.5", "al5") ~ "Citric acid",
      treatment %in% c("at0.5", "at5") ~ "Lactic acid",
      treatment %in% c("am0.5", "am5") ~ "Tartaric acid",
      treatment %in% c("aa0.5", "aa5") ~ "Malic acid",
      treatment %in% c("as0.5", "as5") ~ "Ascorbic acid",
      treatment %in% c("af0.5", "af5") ~ "Sorbic acid",
      treatment %in% c("aasc0.5", "aasc5") ~ "Formic acid",
      treatment == "sl5" ~ "Whey",
      treatment == "lp1" ~ "Milk powder",
      treatment == "ldv10" ~ "Skimmed milk",
      treatment == "lsv5" ~ "Semi-skimmed milk",
      treatment == "lev3" ~ "Whole milk",
      treatment == "lsc5" ~ "Goat milk",
      treatment == "yn1" ~ "Plain yogurt",
      treatment == "ke1" ~ "Kefir",
      treatment == "CaOH2" ~ "Calcium hydroxide",
      treatment == "NaHCO3" ~ "Sodium bicarbonate",
      treatment == "KOH" ~ "Potassium hydroxide",
      treatment == "NH4NO3" ~ "Ammonium nitrate",
      treatment == "KNO3" ~ "Potassium nitrate",
      treatment == "NH4SO4" ~ "Ammonium sulfate",
      treatment == "CuSO4" ~ "Copper sulfate",
      treatment == "CaCO3" ~ "Calcium carbonate",
      TRUE ~ treatment
    ),
    concentration_clean = case_when(
      grepl("0.5", treatment) ~ "Low (0.5%)",
      grepl("5", treatment) & treatment != "sl5" ~ "High (5%)",
      TRUE ~ "Single Dose"
    )
  )

# =========================================================================
# 2. MODEL-ESTIMATED THRESHOLDS FOR CONTROLS
# =========================================================================

model_controls <- lmer(PIG ~ treatment + (1 | inoculation), 
                       weights = total_spores, 
                       data = processed_data %>% filter(treatment %in% c("do", "C-")))

em_controls <- as.data.frame(emmeans(model_controls, ~ treatment, level = 0.90))

dodine_thresh <- em_controls %>% 
  filter(treatment == "do") %>%
  rename(mean_PIG = emmean, ymin = lower.CL, ymax = upper.CL)

c_minus_thresh <- em_controls %>% 
  filter(treatment == "C-") %>%
  rename(mean_PIG = emmean, ymin = lower.CL, ymax = upper.CL)

# =========================================================================
# 3. SPLIT DATASETS FOR EXPERIMENTAL SUBSTANCES
# =========================================================================

experimental_data <- processed_data %>% 
  filter(!treatment %in% c("do", "C-", "C+"))

# Dataset 1: Dual Concentration (8 Organic Acids)
data_analysis_1 <- experimental_data %>%
  filter(substance_clean %in% c("Acetic acid", "Citric acid", "Lactic acid", 
                                "Tartaric acid", "Malic acid", "Ascorbic acid", 
                                "Sorbic acid", "Formic acid"))

# Dataset 2: Single-Dose (Includes Whey and remaining substances)
data_analysis_2 <- experimental_data %>%
  filter(!substance_clean %in% data_analysis_1$substance_clean) # <-- FIXED TYPO HERE

# =========================================================================
# 4. EXECUTE MODELS & EXTRACT ESTIMATED MARGINAL MEANS
# =========================================================================

# Analysis 1 Model (Interaction Model)
model_1 <- lmer(PIG ~ substance_clean * concentration_clean + (1 | inoculation), 
                weights = total_spores, data = data_analysis_1)
em_1 <- emmeans(model_1, ~ substance_clean * concentration_clean, level = 0.90)
df_plot_1 <- as.data.frame(em_1)

# Analysis 2 Model (Single-Factor Model)
model_2 <- lmer(PIG ~ substance_clean + (1 | inoculation), 
                weights = total_spores, data = data_analysis_2)
em_2 <- emmeans(model_2, ~ substance_clean, level = 0.90)
df_plot_2 <- as.data.frame(em_2)




# ==========================================
# PLOT 1: DUAL CONCENTRATION SUBSTANCES
# ==========================================
ggplot(df_plot_1, aes(x = substance_clean, y = emmean, color = concentration_clean)) +
  # Translucent ribbon for Dodine Control
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = dodine_thresh$ymin, ymax = dodine_thresh$ymax), 
            fill = "aquamarine4", alpha = 0.15, inherit.aes = FALSE) +
  geom_hline(yintercept = dodine_thresh$mean_PIG, color = "aquamarine4", linetype = "dashed", size = 0.8) +
  
  # Translucent ribbon for C- (Negative Control)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = c_minus_thresh$ymin, ymax = c_minus_thresh$ymax), 
            fill = "brown2", alpha = 0.15, inherit.aes = FALSE) +
  geom_hline(yintercept = c_minus_thresh$mean_PIG, color = "brown2", linetype = "dashed", size = 0.8) +
  
  # Model data points and 90% CI Error Bars (Dedged to prevent overlapping)
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, position = position_dodge(0.4), size = 0.8) +
  geom_point(position = position_dodge(0.4), size = 3) +
  
  # Clean labels for the reference thresholds
  annotate("text", x = 0.6, y = dodine_thresh$mean_PIG + 3, label = "Dodine Threshold", color = "aquamarine4", fontface = "bold", hjust = 0, size = 3.5) +
  annotate("text", x = 0.6, y = c_minus_thresh$mean_PIG - 3, label = "C- Threshold", color = "brown2", fontface = "bold", hjust = 0, size = 3.5) +
  
  labs(
    title = "Inhibition Percentage (PIG) - Dual Concentration Substances",
    subtitle = "Points represent Model-Adjusted Estimated Marginal Means; Error bars indicate 90% Confidence Intervals",
    x = "Organic Acid Tested",
    y = "Percentage of Inhibition of Germination (PIG %)",
    color = "Concentration Tier"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # <-- FIXED TYPO HERE
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

# ==========================================
# PLOT 2: SINGLE CONCENTRATION SUBSTANCES
# ==========================================
ggplot(df_plot_2, aes(x = reorder(substance_clean, emmean), y = emmean)) +
  # Translucent ribbon for Dodine Control
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = dodine_thresh$ymin, ymax = dodine_thresh$ymax), 
            fill = "aquamarine4", alpha = 0.15, inherit.aes = FALSE) +
  geom_hline(yintercept = dodine_thresh$mean_PIG, color = "aquamarine4", linetype = "dashed", size = 0.8) +
  
  # Translucent ribbon for C- (Negative Control)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = c_minus_thresh$ymin, ymax = c_minus_thresh$ymax), 
            fill = "brown2", alpha = 0.15, inherit.aes = FALSE) +
  geom_hline(yintercept = c_minus_thresh$mean_PIG, color = "brown2", linetype = "dashed", size = 0.8) +
  
  # Model data points and 90% CI Error Bars
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.25, color = "navy", size = 0.8) +
  geom_point(color = "navy", size = 3) +
  
  # Clean labels for reference thresholds
  annotate("text", x = 0.6, y = dodine_thresh$mean_PIG + 3, label = "Dodine Threshold", color = "aquamarine4", fontface = "bold", hjust = 0, size = 3.5) +
  annotate("text", x = 0.6, y = c_minus_thresh$mean_PIG - 3, label = "C- Threshold", color = "brown2", fontface = "bold", hjust = 0, size = 3.5) +
  
  labs(
    title = "Inhibition Percentage (PIG) - Single Dose Treatments",
    subtitle = "Substances sorted dynamically from least effective to most effective (90% Confidence Intervals)",
    x = "Substance / Dairy / Inorganic Salt",
    y = "Percentage of Inhibition of Germination (PIG %)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # <-- FIXED TYPO HERE
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )





Global Bio-Assay Efficacy Matrix
Comprehensive Mixed-Effects Analysis of All Screening Compounds & Organic Acids
Unified Treatment Group	Performance vs. Commercial Fungicide (do)	Performance vs. Untreated Growth (C-)
Δ PIG (%)	p-value	Δ PIG (%)	p-value
Calcium carbonate	3.870	0.9997	−0.468	1.0000
Potassium nitrate	−2.173	1.0000	−6.512	0.9894
Ammonium nitrate	−2.451	1.0000	−6.790	0.9862
Malic acid (0.5%)	−2.720	1.0000	−7.058	0.9854
Sorbic acid (0.5%)	−3.609	0.9998	−7.947	0.9643
Lactic acid (0.5%)	−9.924	0.8872	−14.262	0.5183
Copper sulfate	−11.149	0.7983	−15.487	0.3873
Potassium hydroxide	−11.920	0.7336	−16.258	0.3191
Citric acid (0.5%)	−13.127	0.6398	−17.465	0.2460
Citric acid (5.0%)	−13.648	0.5831	−17.987	0.2056
Ascorbic acid (0.5%)	−15.857	0.3539	−20.195	0.0883
Semi-skimmed milk	−18.063	0.1917	−22.402	0.0357
Lactic acid (5.0%)	−18.431	0.2033	−22.770	0.0420
Ascorbic acid (5.0%)	−18.827	0.1457	−23.165	0.0249
Sorbic acid (5.0%)	−20.003	0.0951	−24.341	0.0142
Formic acid (0.5%)	−20.164	0.0894	−24.503	0.0131
Tartaric acid (5.0%)	−21.651	0.0490	−25.990	0.0060
Acetic acid (0.5%)	−22.920	0.0587	−27.258	0.0090
Calcium hydroxide	−25.685	7.12 × 10−3	−30.023	5.77 × 10−4
Formic acid (5.0%)	−26.951	3.58 × 10−3	−31.289	2.56 × 10−4
Tartaric acid (0.5%)	−29.420	8.40 × 10−4	−33.758	4.75 × 10−5
Sodium bicarbonate	−30.430	4.46 × 10−4	−34.769	2.30 × 10−5
Malic acid (5.0%)	−31.057	3.97 × 10−4	−35.395	2.11 × 10−5
Whole milk	−34.530	4.09 × 10−5	−38.868	1.52 × 10−6
Acetic acid (5.0%)	−36.502	1.53 × 10−5	−40.840	5.73 × 10−7
Plain yogurt	−36.649	2.58 × 10−5	−40.988	1.04 × 10−6
Goat milk	−36.954	4.91 × 10−6	−41.293	1.34 × 10−7
Kefir	−41.662	1.09 × 10−6	−46.000	3.27 × 10−8
Ammonium sulfate	−42.389	4.92 × 10−8	−46.727	8.78 × 10−10
Skimmed milk	−43.653	1.81 × 10−8	−47.991	2.69 × 10−10
Milk powder	−43.869	2.79 × 10−8	−48.208	4.64 × 10−10
Whey	−53.762	7.64 × 10−13	−58.100	9.68 × 10−14




