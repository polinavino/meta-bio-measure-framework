# Reproduce FASD per-signature AUCs (cited: buccal 0.68-0.79, blood van der Laan 0.93-0.96)
setwd("/Users/polina/Documents/BioInfStuff/methylation-biomarker-agreement")
x <- readRDS("data/fasd_all_scores.rds")
sigs <- c("portales","lussier161","lussier183","vanderlaan")   # first 3 buccal, last blood-native
auc <- function(score, label){                                 # label: TRUE=case(FASD)
  r <- rank(score); nc <- sum(label); nk <- sum(!label)
  (sum(r[label]) - nc*(nc+1)/2) / (nc*nk)                       # P(case score > control score)
}
for (co in c("GSE112987","GSE113012")) {
  d <- x[x$cohort==co,]; lab <- d$diagnosis=="FASD"
  cat(sprintf("\n%s: n=%d (%d FASD / %d control)\n", co, nrow(d), sum(lab), sum(!lab)))
  for (s in sigs) cat(sprintf("  %-12s AUC = %.3f\n", s, auc(d[[s]], lab)))
}
cat("\nCited: buccal (portales/lussier) 0.68-0.79; van der Laan (blood) 0.93 discovery / 0.96 replication\n")
