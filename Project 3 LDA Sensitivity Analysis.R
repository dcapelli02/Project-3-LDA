# ==============================================================================
# 0. SETUP LIBRERIE E DATI
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

# --- CARICAMENTO DATI ---
# (Modifica il percorso se necessario)
alz <- read_sas("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.sas7bdat")

# --- PRE-PROCESSING BASE ---
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
# 1. DEFINIZIONE DELLE FUNZIONI DI ANALISI
# ==============================================================================

# --- FUNZIONE A: STANDARD MICE + SHIFT (Delta Adjustment classico) ---
run_standard_shift <- function(data_input, shift_value) {
  
  min_bprs <- 24
  max_bprs <- 168
  
  # Setup MICE
  ini <- mice(data_input, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  meth[paste0("bprs", 1:6)] <- "norm"
  post <- ini$post
  
  # Applica lo shift tramite post-processing
  for (v in paste0("bprs", 1:6)) {
    cmd <- paste0("imp[[j]][, i] <- imp[[j]][, i] + ", shift_value)
    # Aggiungiamo il controllo dei limiti nella stessa stringa di comando
    post[v] <- paste0(cmd, "; imp[[j]][, i] <- pmin(pmax(imp[[j]][, i], ", min_bprs, "), ", max_bprs, ")")
  }
  
  # Imputazione
  imp <- mice(data_input, m = 20, method = meth, post = post, 
              seed = 1234, ridge = 0.001, printFlag = FALSE)
  
  # Creazione Dataset Long
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
  
  # Modello LMER
  fit <- alz_long %>%
    group_by(.imp) %>%
    do(model = lmer(bprs ~ trial + age_std + edu + bmi_std + job + adl_num + 
                      wzc + cdrsb_base + year + age_std:year + edu:year + job:year + 
                      adl_num:year + cdrsb_base:year + (1 + year | .id), 
                    data = ., REML = TRUE,
                    control = lmerControl(optimizer = "bobyqa")))
  
  return(pool(as.mira(as.list(fit$model))))
}

# --- FUNZIONE B: PATTERN MIXTURE + SHIFT (Gold Standard) ---
run_pattern_shift <- function(data_input, shift_value) {
  
  min_bprs <- 24
  max_bprs <- 168
  
  # Crea variabile Pattern (ultimo tempo osservato)
  data_working <- data_input
  data_working$pattern <- apply(data_working[, paste0("bprs", 0:6)], 1, function(x) max(which(!is.na(x))) - 1)
  data_working$pattern <- as.factor(data_working$pattern)
  
  # Setup MICE
  ini <- mice(data_working, maxit = 0, printFlag = FALSE)
  pred <- ini$predictorMatrix
  meth <- ini$method
  post <- ini$post
  
  vars_bprs <- paste0("bprs", 0:6)
  
  # Configura Pattern + Post Processing
  for (i in 2:length(vars_bprs)) {
    current_var <- vars_bprs[i]
    
    # 1. Usa Pattern come predittore
    pred[current_var, "pattern"] <- 1 
    
    # 2. Rimuovi futuri (Logica temporale)
    future_vars <- vars_bprs[(i+1):length(vars_bprs)]
    future_vars <- future_vars[!is.na(future_vars)]
    if (length(future_vars) > 0) pred[current_var, future_vars] <- 0
    
    # 3. Aggiungi SHIFT + LIMITI (MODIFICATO)
    cmd <- paste0("imp[[j]][, i] <- imp[[j]][, i] + ", shift_value)
    post[current_var] <- paste0(cmd, "; imp[[j]][, i] <- pmin(pmax(imp[[j]][, i], ", min_bprs, "), ", max_bprs, ")")
  }
  
  # Imputazione
  imp <- mice(data_working, m = 20, method = meth, predictorMatrix = pred, post = post,
              seed = 1234, ridge = 0.001, printFlag = FALSE)
  
  # Creazione Dataset Long
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
  
  # Modello LMER
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
# 2. ESECUZIONE DELLE ANALISI (LOOP)
# ==============================================================================

shifts_da_testare <- c(0, 2, 4, 6)

# --- LOOP 1: Analisi Standard (MICE classico + Shift) ---
nomi_scenari_std <- paste0("STD_Shift_", shifts_da_testare)
risultati_std <- lapply(shifts_da_testare, function(s) {
  message(paste(">>> Running STANDARD Analysis | Shift:", s))
  run_standard_shift(alz, s)
})
names(risultati_std) <- nomi_scenari_std

# --- LOOP 2: Analisi Pattern (Pattern Mixture + Shift) ---
nomi_scenari_pat <- paste0("PAT_Shift_", shifts_da_testare)
risultati_pat <- lapply(shifts_da_testare, function(s) {
  message(paste(">>> Running PATTERN Analysis | Shift:", s))
  run_pattern_shift(alz, s)
})
names(risultati_pat) <- nomi_scenari_pat

# ==============================================================================
# 3. CREAZIONE TABELLA FINALE E EXPORT (MODIFICATA PER P-VALUE)
# ==============================================================================

# Funzione helper per estrarre dati puliti
extract_data <- function(pool_list, scenario_names) {
  map_df(scenario_names, function(nome) {
    # summary() su un oggetto pool restituisce estimate, std.error e p.value
    summary(pool_list[[nome]], conf.int = TRUE) %>%
      mutate(scenario = nome)
  })
}

# Estrazione dai loop precedenti
df_std <- extract_data(risultati_std, nomi_scenari_std)
df_pat <- extract_data(risultati_pat, nomi_scenari_pat)

# Unione
mega_confronto <- bind_rows(df_std, df_pat)

# Formattazione per Excel: Stima (sd=... - p=...)
# Nota: Ho aggiunto una logica per scrivere "<.001" se il p-value è molto piccolo
tabella_finale_export <- mega_confronto %>%
  mutate(
    p_format = ifelse(p.value < 0.001, "<.001", sprintf("%.3f", p.value)),
    valore = sprintf("%.3f (sd=%.3f - p=%s)", estimate, std.error, p_format)
  ) %>%
  dplyr::select(term, scenario, valore) %>%
  # Ordiniamo le colonne: Prima tutti gli STD, poi tutti i PAT
  mutate(scenario = factor(scenario, levels = c(nomi_scenari_std, nomi_scenari_pat))) %>%
  pivot_wider(names_from = scenario, values_from = valore)

# Visualizzazione anteprima
print(head(tabella_finale_export))

# Salvataggio
write_xlsx(tabella_finale_export, "Sensitivity_Analysis_COMPLETE_With_Pvalues.xlsx")

message("--- ANALISI COMPLETATA: File Excel (con P-Values) creato con successo ---")

