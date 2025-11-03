# ---- ABOUT ----
# Title of the manuscript : 
#    Comparing correlates of daily tobacco use and COVID-19 vaccination refusal 
#    in the working-age Hungarian population
# R-Code      : anonymized
# Closed      : 3rd of November, 2025  


# ---- PACKAGES AND ENVIRONMENT ----
#preparation
  rm(list = ls())
  gc(reset = TRUE)

#packages
  library(data.table)
  library(dplyr)
  
  library(car)
  library(VGAM)
  library(PRROC) 
  

# ---- LOADING INPUT DATA ----
  C19 <- as.data.table(readRDS("C19.RDS")) 

  yvarnames = paste0("y", 1:2)
  xvarnames = paste0("x", 1:21)
 
  
# ---- Descriptive statistics ----
  table2 <- summary(C19)
  print("Table 2.")
  print(table2) 


# ---- Data Imputation ----
#y1,y2 are not imputed, edu coding of x3,x4 considered
  C19_imp <- C19 %>%
    filter(!is.na(y1) & !is.na(y2)) %>%
    mutate(edu_group = case_when(
      x3 == 0 & x4 == 0 ~ "alapfoku",
      x3 == 1 & x4 == 0 ~ "kozepfoku",
      x3 == 1 & x4 == 1 ~ "felsofoku",
    )) %>%
    group_by(edu_group) %>%
    mutate(across(all_of(xvarnames), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    ungroup() %>%
    select(-edu_group) %>%
    select(-sorsz) #remove uid
#remove C19, not to analyze accidentaly
  rm(C19)


# ---- Scaling, for comparable AMEs ----
  C19_imp <- C19_imp %>%
    mutate(across(.cols = setdiff(names(.)[sapply(., is.numeric)], c("y1", "y2")),
                  .fns = ~ as.numeric(scale(.))))  


# ---- Multivariate modeling ----  
#y1~y2
  print("Phi(y1~y2)=")
  cor(C19_imp$y1, C19_imp$y2)
  
  print("JI(y1~y2)=")
  jaccard_index <- function(x, y) { sum(x & y) / sum(x | y) }
  jaccard_index(C19_imp$y1, C19_imp$y2)

#VIF  
  lm_vif_y1 <- glm(y1 ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
      x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + x21, 
                data = C19_imp, 
                family = binomial)
  
  lm_vif_y2 <- glm(y2 ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
      x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + x21, 
                   data = C19_imp, 
                   family = binomial)
  
  vif_y1 <- car::vif(lm_vif_y1)
  vif_y2 <- car::vif(lm_vif_y2)
  
  print("VIFs")
  summary(cbind(vif_y1, vif_y2))
  
# def y
  y_mat <- as.matrix(C19_imp[, c("y1", "y2")])

# def null_mod
  null_mod <-  vglm(
    y_mat ~ 1,
    family = binomialff(multiple.responses = TRUE),
    data = C19_imp
  )
  print("Null modell:")
  summary(null_mod)
  

#full mod
  fit_mv_logit_vgam <- VGAM::vglm(
    y_mat ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 
      x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 + x21,
    family = binomialff(multiple.responses = TRUE),
    data = C19_imp
  )
  print("Full modell:")
  summary(fit_mv_logit_vgam)

#aORs
  se <- sqrt(diag(vcov(fit_mv_logit_vgam)))
  
  log_odds_vec <- coef(fit_mv_logit_vgam)
  lower <- log_odds_vec - 1.96 * se
  upper <- log_odds_vec + 1.96 * se
  
  aor_ci <- data.frame(
    aOR = round(exp(log_odds_vec), 2),
    aOR_lower = round(exp(lower), 2),
    aOR_upper = round(exp(upper), 2)
  )
  aor_ci$str=paste0(aor_ci$aOR,"(",aor_ci$aOR_lower,"–", aor_ci$aOR_upper,")")
  
  print("Adjusted Odds Ratios with 95% CI:")
  print(aor_ci)
  
#AMEs
  pred <- predict(fit_mv_logit_vgam, type = "response")  
  coefs <- coef(fit_mv_logit_vgam, matrix = TRUE) 
  avg_grad <- colMeans(pred * (1 - pred))
  ame_matrix <- sweep(coefs, 2, avg_grad, FUN = "*")
  ame_xvars <- ame_matrix[rownames(coefs) %in% xvarnames, , drop = FALSE]
  print("AMEs:")
  print(round(ame_xvars,3))

#LL
  ll_null <- logLik(null_mod)
  ll_full <- logLik(fit_mv_logit_vgam)
  
  df_null <- length(coef(null_mod))
  df_full <- length(coef(fit_mv_logit_vgam))
  
  LR_stat <- 2 * (as.numeric(ll_full) - as.numeric(ll_null))  
  df_diff <- df_full - df_null
  
  pval <- pchisq(LR_stat, df = df_diff, lower.tail = FALSE)
  
  print("Likelihood-ratio test:")
  print(paste("LogLik null model:", round(as.numeric(ll_null), 2)))
  print(paste("LogLik full model:", round(as.numeric(ll_full), 2)))
  print(paste("LR Chi-square:", round(LR_stat, 2)))
  print(paste("df:", df_diff))
  print(paste("p-value:", signif(pval, 3)))
  
#AIC
  print("AIC-null-modell:")
  print(round(AIC(null_mod),2))
  print("AIC-full-modell:")
  print(round(AIC(fit_mv_logit_vgam),2))
  
#ROC-AUC and PR-AUC
  pred_y1 <- pred[, 1]
  pred_y2 <- pred[, 2]

prroc_metrics <- function(scores, labels) {
  s0 <- scores[labels == 0]  # negatives
  s1 <- scores[labels == 1]  # positives

  # IMPORTANT: for ROC, pass POSITIVES FIRST
  roc <- PRROC::roc.curve(scores.class0 = s1, scores.class1 = s0, curve = TRUE)
  pr  <- PRROC::pr.curve(scores.class0 = s0, scores.class1 = s1, curve = TRUE)

  list(roc_auc = roc$auc, pr_auc = pr$auc.integral,
       pr_curve = pr, roc_curve = roc)
}

m1 <- prroc_metrics(pred_y1, C19_imp$y1)
cat("ROC-AUC for y1:", round(m1$roc_auc, 2), "\n")
cat("PR-AUC  for y1:", round(m1$pr_auc,  2), "\n")
cat("PR-AUC baseline (prevalence) for y1:", round(mean(C19_imp$y1), 2), "\n")

m2 <- prroc_metrics(pred_y2, C19_imp$y2)
cat("ROC-AUC for y2:", round(m2$roc_auc, 2), "\n")
cat("PR-AUC  for y2:", round(m2$pr_auc,  2), "\n")
cat("PR-AUC baseline (prevalence) for y2:", round(mean(C19_imp$y2), 2), "\n")

stop("Code run completed.")  
