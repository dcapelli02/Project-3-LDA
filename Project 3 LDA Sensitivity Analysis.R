# ==============================================================================
# 0. SETUP LIBRARIES AND DATA
# ==============================================================================
library(haven)
library(sas7bdat)
library(dplyr)
library(tidyr)
library(mice)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(purrr)
library(writexl)

# --- DATA LOADING ---
# (Modify the path if necessary)
alz <- read_sas("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.sas7bdat")

# --- BASIC PRE-PROCESSING ---
alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)
alz$cdrsb_base <- alz$cdrsb0
alz$n_obs_data <- rowSums(!is.na(alz[, paste0("bprs", 0:6)]))

# ==============================================================================
# 1. DEFINITION OF ANALYSIS FUNCTIONS
# ==============================================================================

# --- FUNCTION A: STANDARD MICE + SHIFT (Classic Delta Adjustment) ---
run_standard_shift <- function(data_input, shift_value) {
  
  min_bprs <- 24
  max_bprs <- 168
  
  # MICE Setup
  ini <- mice(data_input, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  meth[paste0("bprs", 1:6)] <- "norm"
  post <- ini$post
  
  # Apply shift via post-processing
  for (v in paste0("bprs", 1:6)) {
    cmd <- paste0("imp[[j]][, i] <- imp[[j]][, i] + ", shift_value)
    # Add boundary check within the same command string
    post[v] <- paste0(cmd, "; imp[[j]][, i] <- pmin(pmax(imp[[j]][, i], ", min_bprs, "), ", max_bprs, ")")
  }
  
  # Imputation
  imp <- mice(data_input, m = 20, method = meth, post = post, 
              seed = 1234, ridge = 0.001, printFlag = FALSE)
  
  # Create Long Dataset
  alz_long <- complete(imp, action = "long") %>%
    pivot_longer(
      cols = matches("^(bprs|)[0-6]$"),
      names_to = c(".value", "TIME"),
      names_pattern = "([a-z]+)([0-6])"
    ) %>%
    mutate(year = as.numeric(TIME)) %>%
    group_by(.imp) %>%
    mutate(
      age_std = (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE),
      bmi_std = (bmi - mean(bmi, na.rm = TRUE)) / sd(bmi, na.rm = TRUE)
    ) %>%
    ungroup()
  
  print(summary(alz_long$bprs))
  
  # LMER Model
  fit <- alz_long %>%
    group_by(.imp) %>%
    do(model = lmer(bprs ~ trial + age_std + edu + bmi_std + job + adl_num + 
                      wzc + cdrsb_base + year + age_std:year + edu:year + job:year + 
                      adl_num:year + cdrsb_base:year + (1 + year | .id), 
                    data = ., REML = TRUE,
                    control = lmerControl(optimizer = "bobyqa")))
  
  return(pool(as.mira(as.list(fit$model))))
}

# --- FUNCTION B: PATTERN MIXTURE + SHIFT (Gold Standard) ---
run_pattern_shift <- function(data_input, shift_value) {
  
  min_bprs <- 24
  max_bprs <- 168
  
  # Create Pattern variable (last observed time)
  data_working <- data_input
  data_working$pattern <- apply(data_working[, paste0("bprs", 0:6)], 1, function(x) max(which(!is.na(x))) - 1)
  data_working$pattern <- as.factor(data_working$pattern)
  
  # MICE Setup
  ini <- mice(data_working, maxit = 0, printFlag = FALSE)
  pred <- ini$predictorMatrix
  meth <- ini$method
  post <- ini$post
  
  vars_bprs <- paste0("bprs", 0:6)
  
  # Configure Pattern + Post Processing
  for (i in 2:length(vars_bprs)) {
    current_var <- vars_bprs[i]
    
    # 1. Use Pattern as predictor
    pred[current_var, "pattern"] <- 1 
    
    # 2. Remove future variables (Temporal logic)
    future_vars <- vars_bprs[(i+1):length(vars_bprs)]
    future_vars <- future_vars[!is.na(future_vars)]
    if (length(future_vars) > 0) pred[current_var, future_vars] <- 0
    
    # 3. Add SHIFT + LIMITS (MODIFIED)
    cmd <- paste0("imp[[j]][, i] <- imp[[j]][, i] + ", shift_value)
    post[current_var] <- paste0(cmd, "; imp[[j]][, i] <- pmin(pmax(imp[[j]][, i], ", min_bprs, "), ", max_bprs, ")")
  }
  
  # Imputation
  imp <- mice(data_working, m = 20, method = meth, predictorMatrix = pred, post = post,
              seed = 1234, ridge = 0.001, printFlag = FALSE)
  
  # Create Long Dataset
  alz_long <- complete(imp, action = "long") %>%
    pivot_longer(
      cols = matches("^(bprs|)[0-6]$"),
      names_to = c(".value", "TIME"),
      names_pattern = "([a-z]+)([0-6])"
    ) %>%
    mutate(year = as.numeric(TIME)) %>%
    group_by(.imp) %>%
    mutate(
      age_std = (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE),
      bmi_std = (bmi - mean(bmi, na.rm = TRUE)) / sd(bmi, na.rm = TRUE)
    ) %>%
    ungroup()
  
  print(summary(alz_long$bprs))
  
  # LMER Model
  fit <- alz_long %>%
    group_by(.imp) %>%
    do(model = lmer(bprs ~ trial + age_std + edu + bmi_std + job + adl_num + 
                      wzc + cdrsb_base + year + age_std:year + edu:year + job:year + 
                      adl_num:year + cdrsb_base:year + (1 + year | .id), 
                    data = ., REML = TRUE,
                    control = lmerControl(optimizer = "bobyqa")))
  
  return(pool(as.mira(as.list(fit$model))))
}

# ==============================================================================
# 2. RUNNING THE ANALYSIS (LOOP)
# ==============================================================================

shifts_to_test <- c(0, 2, 4, 6)

# --- LOOP 1: Standard Analysis (Classic MICE + Shift) ---
scenario_names_std <- paste0("STD_Shift_", shifts_to_test)
results_std <- lapply(shifts_to_test, function(s) {
  message(paste(">>> Running STANDARD Analysis | Shift:", s))
  run_standard_shift(alz, s)
})
names(results_std) <- scenario_names_std

# --- LOOP 2: Pattern Analysis (Pattern Mixture + Shift) ---
scenario_names_pat <- paste0("PAT_Shift_", shifts_to_test)
results_pat <- lapply(shifts_to_test, function(s) {
  message(paste(">>> Running PATTERN Analysis | Shift:", s))
  run_pattern_shift(alz, s)
})
names(results_pat) <- scenario_names_pat

# ==============================================================================
# 3. FINAL TABLE CREATION AND EXPORT (MODIFIED FOR P-VALUE)
# ==============================================================================

# Helper function to extract clean data
extract_data <- function(pool_list, scenario_names) {
  map_df(scenario_names, function(nome) {
    # summary() on a pool object returns estimate, std.error and p.value
    summary(pool_list[[nome]], conf.int = TRUE) %>%
      mutate(scenario = nome)
  })
}

# Extraction from previous loops
df_std <- extract_data(results_std, scenario_names_std)
df_pat <- extract_data(results_pat, scenario_names_pat)

# Binding
mega_comparison <- bind_rows(df_std, df_pat)

# Excel Formatting: Estimate (sd=... - p=...)
# Note: Added logic to write "<.001" if p-value is very small
final_export_table <- mega_comparison %>%
  mutate(
    p_format = ifelse(p.value < 0.001, "<.001", sprintf("%.3f", p.value)),
    value = sprintf("%.3f (sd=%.3f - p=%s)", estimate, std.error, p_format)
  ) %>%
  dplyr::select(term, scenario, value) %>%
  # Order columns: All STD first, then all PAT
  mutate(scenario = factor(scenario, levels = c(scenario_names_std, scenario_names_pat))) %>%
  pivot_wider(names_from = scenario, values_from = value)

# Preview visualization
print(head(final_export_table))

# Saving
write_xlsx(final_export_table, "Sensitivity_Analysis_COMPLETE_With_Pvalues.xlsx")

message("--- ANALYSIS COMPLETED: Excel File (with P-Values) successfully created ---")
