rm(list=ls())
setwd("E:/Academic/4-1/Thesis/Serum and Plasma Metabolomic/plasma/")
load("dem_t_test.RData")
load("dem_k_wallis.RData")
load("dem_wilcoxon.RData")

resAll = data.frame("Metabolites" = resKW$Metabolites, "tt.p.value" = resTT$p.value, "tt.adjusted.p.value" = resTT$Adjusted.p.value, "kw.p.value" = resKW$p.value, "kw.adjusted.p.value" = resKW$Adjusted.p.value, "mww.p.value" = resMWW$p.value, "mww.adjusted.p.value" = resMWW$Adjusted.p.value)

kept = (resAll$tt.adjusted.p.value < 0.05) | (resAll$kw.adjusted.p.value < 0.05) | (resAll$mww.adjusted.p.value < 0.05)

resAll = resAll[kept,]

write.csv(resAll, file="plasma_unioned_table.csv")
