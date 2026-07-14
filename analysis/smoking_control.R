## Smoking domain controls:
## (1) per-cohort reproducibility of the C3 classification-agreement pattern (G4 linchpin),
## (2) robustness of the former-smoker low-agreement to the threshold choice,
## (3) a threshold-FREE continuous disagreement metric (inter-signature rank spread) by group.
suppressMessages({library(dplyr)})
setwd("/Users/polina/Documents/BioInfStuff/methylation-biomarker-agreement")
s50 <- readRDS("data/scores_GSE50660.rds")
s42 <- readRDS("data/scores_GSE42861.rds")
sigs <- c("ahrr","epismoke","joehanes","epitob")   # all: lower = more exposure

norm_status <- function(df){
  df$grp <- dplyr::recode(as.character(df$smoking),
    "0"="never","1"="former","2"="current",
    "never"="never","ex"="former","former"="former","current"="current",
    "occasional"="occasional", .default=NA_character_)
  df
}
s50 <- norm_status(s50); s42 <- norm_status(s42)

## ---- (1)+(2) classification agreement, per cohort, threshold sweep ----
agree_tab <- function(df, mult){
  never <- df[df$grp=="never",]
  thr <- sapply(sigs, function(s) mean(never[[s]],na.rm=TRUE) - mult*sd(never[[s]],na.rm=TRUE))
  # classify "smoker" if score below threshold (lower=more exposure)
  cls <- sapply(sigs, function(s) df[[s]] < thr[s])
  nsmk <- rowSums(cls)                     # how many signatures call "smoker"
  df$allagree <- (nsmk==0 | nsmk==4)
  df %>% group_by(grp) %>% summarise(agree=mean(allagree), n=n(), .groups="drop")
}
cat("=== (1) Per-cohort C3 agreement at threshold never_mean - 2*SD ===\n")
cat("GSE50660:\n"); print(as.data.frame(agree_tab(s50,2)))
cat("GSE42861:\n"); print(as.data.frame(agree_tab(s42,2)))

cat("\n=== (2) Threshold robustness (never/current/former agreement, both cohorts) ===\n")
cat("mult | GSE50660 never curr former | GSE42861 never curr former\n")
getg <- function(a,g) { v<-a$agree[a$grp==g]; ifelse(length(v),v,NA) }
for(m in c(1,1.5,2,2.5,3)){
  a50 <- agree_tab(s50,m); a42 <- agree_tab(s42,m)
  cat(sprintf("%.1f  |   %.3f %.3f %.3f  |   %.3f %.3f %.3f\n", m,
      getg(a50,"never"),getg(a50,"current"),getg(a50,"former"),
      getg(a42,"never"),getg(a42,"current"),getg(a42,"former")))
}

## ---- (3) threshold-FREE: inter-signature rank spread per sample, by group ----
## rank each signature within cohort (rank1 = most exposed = lowest score), std across 4 sigs
rankspread <- function(df){
  Rk <- sapply(sigs, function(s) rank(df[[s]], ties.method="average"))
  df$spread <- apply(Rk,1,sd)
  df %>% group_by(grp) %>% summarise(mean_spread=mean(spread), n=n(), .groups="drop")
}
cat("\n=== (3) Threshold-free inter-signature rank spread (std of 4 ranks), per cohort ===\n")
cat("GSE50660:\n"); print(as.data.frame(rankspread(s50)))
cat("GSE42861:\n"); print(as.data.frame(rankspread(s42)))
cat("\n(If former-smoker spread is NOT elevated vs never/current here, the C3 U-shape is a\n",
    "threshold artifact; if it IS elevated and reproduces across cohorts, concentration is real.)\n")
