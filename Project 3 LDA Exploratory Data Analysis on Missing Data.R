library(haven)      # Dades SAS
library(dplyr)      # Manipulació Dades
library(tidyr)      # Transformacions (pivot_longer)
library(ggplot2)    # Gràfics
library(naniar)     # Missing Data 
library(VIM)        # Visualització Patrons
library(nlme)       # LMM 
library(lme4)       # GLMM 
library(geepack)    # GEE i WGEE 
library(mice)       # Sensitivity Analysis
library(corrplot)   # Matriu Correlació
library(scales)     # Etiquetes llegibles


#================================================================#

# DATASET

alz <- read_sas("C:/Users/win11/Documents/Guillem/Erasmus/Assignatures/Longitudinal Data Analysis/Project-3-LDA-Personal/alzheimer25.sas7bdat")

# Factor variables

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)
alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))

# Bins

alz$cdrsb_bin0 <- ifelse(alz$cdrsb0 > 10, 1, 0)
alz$cdrsb_bin1 <- ifelse(alz$cdrsb1 > 10, 1, 0)
alz$cdrsb_bin2 <- ifelse(alz$cdrsb2 > 10, 1, 0)
alz$cdrsb_bin3 <- ifelse(alz$cdrsb3 > 10, 1, 0)
alz$cdrsb_bin4 <- ifelse(alz$cdrsb4 > 10, 1, 0)
alz$cdrsb_bin5 <- ifelse(alz$cdrsb5 > 10, 1, 0)
alz$cdrsb_bin6 <- ifelse(alz$cdrsb6 > 10, 1, 0)

# Baseline

alz$ab_base <- alz$abpet0
alz$tau_base <- alz$taupet0
alz$cdrsb_base <- alz$cdrsb0
alz$bprs_base <- alz$bprs0

#================================================================#

# LONGITUDINAL DATASET

alz_long <- alz %>%
  pivot_longer(
    cols = matches(".*\\d+$"),     # Selecciona qualsevol columna que acabi en número
    names_to = c(".value", "time"),
    names_pattern = "(.*)(\\d+)$"  
  ) %>%
  mutate(
    time = as.numeric(time),       
    id = as.factor(patid)         
  ) %>%
  arrange(id, time)

head(alz_long)
str(alz_long)

#================================================================#

# EXPLORATORY ANALYSIS ON MISSING DATA

# General Behavior

alz %>%
  select(matches("\\d+$")) %>% 
  vis_miss() +
  labs(title = "Missingness Pattern") +
  theme(axis.text.x = element_text(angle = 90))

alz %>%
  select(matches("cdrsb\\d+$")) %>% 
  vis_miss() +
  labs(title = "Missingness Pattern on cdrsb") +
  theme(axis.text.x = element_text(angle = 90))

alz_long %>%
  select(cdrsb, bprs, abpet, taupet) %>% 
  gg_miss_var(show_pct = TRUE) + 
  labs(title = "% Missing Data",
       y = "Missings' Percentage")

cat("The missingness pattern is equal in every longitudinal variable.")

# Dropout

dropout_table <- alz_long %>%
  group_by(time) %>%
  summarise(
    N_Total = n(),
    Observed_CDRSB = sum(!is.na(cdrsb)),
    Missing_CDRSB = sum(is.na(cdrsb)),
    Pct_Missing = round(mean(is.na(cdrsb)) * 100, 2)
  )

print(dropout_table)

ggplot(dropout_table, aes(x = time, y = Observed_CDRSB)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(size = 3) +
  labs(title = "Number of observed patients evolution",
       x = "Year", y = "Observed Patients") +
  theme_bw()

# Monotoneness

check_intermittent <- alz_long %>%
  select(id, time, cdrsb) %>%
  arrange(id, time) %>%
  group_by(id) %>%
  summarise(
    is_missing = list(is.na(cdrsb)),
    pattern = paste(as.integer(!is.na(cdrsb)), collapse = "")
  )
table(check_intermittent$pattern)

cat("We conclude that missingness follows a monotone pattern.")

# Missingness per group

# Trial
alz %>%
  select(trial, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = trial) +
  labs(title = "% Missings depending on Trial",
       y = "% Observed Patients",
       x = "Trial")

# Sex
alz %>%
  select(sex, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = sex) +
  labs(title = "% Missings depending on Sex",
       y = "% Observed Patients",
       x = "Sex")

# 0 - men, 1 - women

# Age
alz %>%
  select(age, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = age) +
  labs(title = "% Missings depending on Age",
       y = "% Observed Patients",
       x = "Age")

# Edu
alz %>%
  select(edu, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = edu) +
  labs(title = "% Missings depending on Edu",
       y = "% Observed Patients",
       x = "Edu")

# BMI
alz %>%
  select(bmi, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = bmi) +
  labs(title = "% Missings depending on BMI",
       y = "% Observed Patients",
       x = "BMI")

# Inkomen
alz %>%
  select(inkomen, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = inkomen) +
  labs(title = "% Missings depending on Inkomen",
       y = "% Observed Patients",
       x = "Inkomen")

# Job
alz %>%
  select(job, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = job) +
  labs(title = "% Missings depending on Job",
       y = "% Observed Patients",
       x = "Job")

# ADL
alz %>%
  select(adl, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = adl) +
  labs(title = "% Missings depending on ADL",
       y = "% Observed Patients",
       x = "ADL")

# WZC
alz %>%
  select(wzc, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = wzc) +
  labs(title = "% Missings depending on WZC",
       y = "% Observed Patients",
       x = "WZC")

# 0 - home, 1 - residence

### COMMON PLOT ###
library(patchwork)

scala_uniforme <- scale_fill_viridis_c(limits = c(0, 100), name = "% Missing")

# 1. Sex
p1 <- alz %>%
  select(sex, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = sex) +
  scala_uniforme +
  labs(title = "Sex", y = "% Observed Patients", x = "")

# 2. Age
p2 <- alz %>%
  select(age, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = age) +
  scala_uniforme +
  labs(title = "Age", y = "% Observed Patients", x = "")

# 3. Edu
p3 <- alz %>%
  select(edu, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = edu) +
  scala_uniforme +
  labs(title = "Edu", y = "% Observed Patients", x = "")

# 4. BMI
p4 <- alz %>%
  select(bmi, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = bmi) +
  scala_uniforme +
  labs(title = "BMI", y = "% Observed Patients", x = "")

# 5. Inkomen
p5 <- alz %>%
  select(inkomen, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = inkomen) +
  scala_uniforme +
  labs(title = "Inkomen", y = "% Observed Patients", x = "")

# 6. Job
p6 <- alz %>%
  select(job, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = job) +
  scala_uniforme +
  labs(title = "Job", y = "% Observed Patients", x = "")

# 7. ADL
p7 <- alz %>%
  select(adl, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = adl) +
  scala_uniforme +
  labs(title = "ADL", y = "% Observed Patients", x = "")

# 8. WZC
p8 <- alz %>%
  select(wzc, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = wzc) +
  scala_uniforme +
  labs(title = "WZC", y = "% Observed Patients", x = "")

# Final Combination
grafico_finale <- (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

# Stampa
print(grafico_finale)

###

# MEAN PROFILES BY DROPOUT TIME

dropout_info <- alz_long %>%
  group_by(id) %>%
  summarise(
    first_missing = ifelse(any(is.na(cdrsb)), 
                           as.character(min(time[is.na(cdrsb)])),
                           "Complete")
  )

alz_patterns <- left_join(alz_long, dropout_info, by = "id")

# CDRSB
alz_patterns %>%
  group_by(time, first_missing) %>%
  summarise(mean_cdrsb = mean(cdrsb, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = time, y = mean_cdrsb, color = as.factor(first_missing))) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "CDRSB Mean's evolution depending on dropout time",
       x = "Year",
       y = "CDRSB Mean",
       color = "Dropout Time") +
  theme_bw()

# BPRS
alz_patterns %>%
  group_by(time, first_missing) %>%
  summarise(mean_bprs = mean(bprs, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = time, y = mean_bprs, color = as.factor(first_missing))) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "BPRS Mean's evolution depending on dropout time",
       x = "Year",
       y = "BPRS Mean",
       color = "Dropout Time") +
  theme_bw()

# TAUPET
alz_patterns %>%
  group_by(time, first_missing) %>%
  summarise(mean_taupet = mean(taupet, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = time, y = mean_taupet, color = as.factor(first_missing))) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "TAUPET Mean's evolution depending on dropout time",
       x = "Year",
       y = "TAUPET Mean",
       color = "Dropout Time") +
  theme_bw()

# ABPET
alz_patterns %>%
  group_by(time, first_missing) %>%
  summarise(mean_abpet = mean(abpet, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = time, y = mean_abpet, color = as.factor(first_missing))) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "ABPET Mean's evolution depending on dropout time",
       x = "Year",
       y = "ABPET Mean",
       color = "Dropout Time") +
  theme_bw()

cat("This may suggest us that we do not have a MCAR case,")
cat("Because patients who drop out have higher means than the completers.")

# DROPOUT DEPENDING ON PREVIOUS VALUE

alz_long_lag <- alz_long %>%
  group_by(id) %>%
  mutate(
    is_missing_next = lead(is.na(cdrsb)) 
  ) %>%
  filter(time < max(time)) 

ggplot(alz_long_lag, aes(x = as.factor(time), y = cdrsb, fill = is_missing_next)) +
  geom_boxplot() +
  labs(title = "Evidence for MAR in CDRSB: Observed vs Missing at next visit",
       x = "Year",
       y = "CDRSB",
       fill = "Missing at t+1?") +
  scale_fill_manual(values = c("skyblue", "salmon"), na.value="grey") +
  theme_minimal()

ggplot(alz_long_lag, aes(x = as.factor(time), y = bprs, fill = is_missing_next)) +
  geom_boxplot() +
  labs(title = "Evidence for MAR in BPRS: Observed vs Missing at next visit",
       x = "Year",
       y = "BPRS",
       fill = "Missing at t+1?") +
  scale_fill_manual(values = c("skyblue", "salmon"), na.value="grey") +
  theme_minimal()

ggplot(alz_long_lag, aes(x = as.factor(time), y = taupet, fill = is_missing_next)) +
  geom_boxplot() +
  labs(title = "Evidence for MAR in TAUPET: Observed vs Missing at next visit",
       x = "Year",
       y = "TAUPET",
       fill = "Missing at t+1?") +
  scale_fill_manual(values = c("skyblue", "salmon"), na.value="grey") +
  theme_minimal()

ggplot(alz_long_lag, aes(x = as.factor(time), y = abpet, fill = is_missing_next)) +
  geom_boxplot() +
  labs(title = "Evidence for MAR in ABPET: Observed vs Missing at next visit",
       x = "Year",
       y = "ABPET",
       fill = "Missing at t+1?") +
  scale_fill_manual(values = c("skyblue", "salmon"), na.value="grey") +
  theme_minimal()

cat("This suggest MAR, because the red boxes are higher than the blue ones.")
cat("The patients that are worse leave the study.")

# DROPOUT MODEL

alz_dropout_model <- alz_long %>%
  group_by(id) %>%
  mutate(
    is_missing = is.na(cdrsb),          
    prev_cdrsb = lag(cdrsb),            
    prev_bprs = lag(bprs) ,
    prev_taupet = lag(taupet),
    prev_abpet = lag(abpet)
  ) %>%
  filter(time > 0) %>% 
  ungroup()

model_dropout <- glm(is_missing ~ prev_cdrsb + prev_bprs + prev_taupet + prev_abpet +
                                  trial + sex + age + edu + bmi + inkomen + job + adl_num + wzc + time, 
                     data = alz_dropout_model, 
                     family = binomial)

summary(model_dropout)

# Interpretation
cat("Since we have that the p-value of prev_cdrsb <2e-16***, we have evidence of MAR.")
cat("Since the estimate for prev_cdrsb = 0.1209 > 0, patients with higher cdrsb tend to dropout more.")

# MCAR TEST

mcar_cdrsb <- mcar_test(alz %>% select(matches("cdrsb\\d+$")))
mcar_bprs <- mcar_test(alz %>% select(matches("bprs\\d+$")))
mcar_taupet <- mcar_test(alz %>% select(matches("taupet\\d+$")))
mcar_abpet <- mcar_test(alz %>% select(matches("abpet\\d+$")))

print(mcar_cdrsb)
print(mcar_bprs)
print(mcar_taupet)
print(mcar_abpet)

cat("Since the p-values of these tests are ~ 0, we conclude that we do not have a MCAR case.")

# COMPARE WITH BASELINE: COMPLETERS VS DROPOUTS 

completers_id <- alz_long %>%
  filter(time == 6 & !is.na(cdrsb)) %>%
  pull(patid)

alz_comparison <- alz %>%
  mutate(
    Status = ifelse(patid %in% completers_id, "Completer", "Dropout")
  )

alz_comparison %>%
  group_by(Status) %>%
  summarise(
    N = n(),
    CDRSB_Base_Mean = mean(cdrsb0, na.rm=TRUE),
    BPRS_Base_Mean = mean(bprs0, na.rm=TRUE),
    TAUPET_Base_Mean = mean(taupet0, na.rm=TRUE),
    ABPET_Base_Mean = mean(abpet0, na.rm=TRUE),
    Sex_Mean = mean(as.numeric(sex)-1, na.rm = TRUE),
    Age_Mean = mean(age, na.rm=TRUE),
    BMI_Mean = mean(bmi, na.rm = TRUE),
    Inkomen_Mean = mean(inkomen, na.rm = TRUE),
    Job_Mean = mean(as.numeric(job)-1, na.rm = TRUE),
    ADL_Mean = mean(adl_num, na.rm = TRUE),
    WZC_Mean = mean(as.numeric(wzc)-1, na.rm = TRUE)

  )

cat("The descriptive comparison between Dropouts (N=742) and Completers (N=511)")
cat("reveals that the dropout group had a higher mean Age (76.4 vs. 66.7 years)")
cat("and higher baseline BPRS scores (83.0 vs. 64.4),")
cat("while conversely presenting a lower mean baseline CDRSB (6.45 vs. 7.14) and lower ADL scores (6.04 vs. 10.5).")

# SPAGHETTI PLOTS

set.seed(123)
random_ids <- sample(unique(alz_long$patid), 30)

sample_data <- alz_long %>% 
  filter(patid %in% random_ids)

# CDRSB
ggplot(sample_data, aes(x = time, y = cdrsb_bin, group = patid, color = as.factor(patid))) +
  geom_line(alpha = 0.7) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") + 
  labs(title = "Trajectories of CDRSB",
       y = "CDRSB", x = "Year")

ggplot(sample_data, aes(x = time, y = cdrsb_bin, group = patid)) +
  geom_line(alpha = 0.2) +  
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "red", size = 1.2) + # La tendència
  labs(title = "Individual Trajectories + Mean",
       y = "CDRSB_bin", x = "Year") +
  theme_bw()

# BPRS
ggplot(sample_data, aes(x = time, y = bprs, group = patid, color = as.factor(patid))) +
  geom_line(alpha = 0.7) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") + 
  labs(title = "Trajectories of BPRS",
       y = "BPRS", x = "Year")

# CORRELATION MATRIX

# CDRSB
cor_data <- alz %>%
  select(matches("cdrsb_bin\\d+$")) %>% 
  na.omit() 

cor_matrix <- round(cor(cor_data), 2)
print(cor_matrix)

corrplot(cor_matrix, 
         title = "CDRSB Correlation Matrix",
         method = "number", type = "upper")

# BPRS
cor_data <- alz %>%
  select(matches("bprs\\d+$")) %>% 
  na.omit() 

cor_matrix <- round(cor(cor_data), 2)
print(cor_matrix)

corrplot(cor_matrix, 
         title = "BPRS Correlation Matrix",
         method = "number", type = "upper")

# DISTRIBUTION

# CDRSB
ggplot(alz_long, aes(x = as.factor(cdrsb_bin), fill = as.factor(cdrsb_bin))) +
  geom_bar() +                                
  facet_wrap(~time) + 
  labs(title = "CDRSB_bin", 
       y = "Number of Patients",
       x = "CDRSB_bin",
       fill = "State") +
  theme_bw()

# BPRS
ggplot(alz_long, aes(x = bprs)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  facet_wrap(~time) + 
  labs(title = "BPRS's Distribution", y = "Number of Patients") +
  theme_bw()


dropout_check <- alz_long %>%
  mutate(observed = !is.na(cdrsb)) %>%  
  group_by(id) %>%
  arrange(time) %>%
  summarise(
    first_missing = ifelse(any(!observed),
                           min(time[!observed]),
                           NA),
    non_monotone = ifelse(
      is.na(first_missing),
      FALSE,
      any(observed & time > first_missing)
    ),
    .groups = "drop"
  )

cat("Table of Non-Monotone patterns (FALSE = Monotone, TRUE = Intermittent):")
print(table(dropout_check$non_monotone))


# Evolution of CDRSB by Dropout Time
cdrsb_bin_analysis <- alz_long %>%
  group_by(id) %>%
  mutate(dropout_time = ifelse(all(!is.na(cdrsb)), 
                               7,  # <--- Als Completers els diem "7" (o el que vulguis)
                               min(time[is.na(cdrsb)]))) %>% # Als altres, el moment que marxen
  ungroup()

cdrsb_prop_group <- cdrsb_bin_analysis %>%
  group_by(dropout_time, time) %>%
  summarise(prop_positive = mean(cdrsb_bin == 1, na.rm = TRUE), .groups = "drop")

# Plot
ggplot(cdrsb_prop_group, aes(x = time, y = prop_positive, color = factor(dropout_time))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_discrete(
    name = "Dropout Time", 
    labels = c("1", "2", "3", "4", "5", "6", "Completers")
  ) +
  labs(
    x = "Visit (Year)",
    y = "Proportion with CDRSB > 10",
    color = "Dropout Time",
    title = "Evolution of Binary CDRSB Outcome by Dropout Time"
  ) +
  theme_minimal()


# Observed vs. Missing - Dichotomized CDRSB
cdrsb_bin_prop <- alz_long_lag %>%
  group_by(time, is_missing_next) %>%
  summarise(prop_positive = mean(cdrsb_bin == 1, na.rm = TRUE), .groups = "drop")

# Plot
ggplot(cdrsb_bin_prop, aes(x = as.factor(time), y = prop_positive, color = is_missing_next, group = is_missing_next)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("skyblue", "salmon"), labels = c("Observed", "Missing")) +
  labs(
    x = "Year",
    y = "Proportion with CDRSB > 10",
    color = "Status at t+1",
    title = "Evidence for MAR: Binary CDRSB (Observed vs Missing Next)"
  ) +
  theme_minimal()


