rm(list=ls())

setwd("E:/Academic/4-1/Thesis/Serum and Plasma Metabolomic/plasma/")

data <- read.csv("plusmaLC.csv", row.names=1)

alpha <- 0.05

cnt <- 0

names <- vector('character')
pvals <- vector('numeric')

for (i in 1:nrow(data)){
  vec1 = as.numeric(data[i:i,1:41])
  vec2 = as.numeric(data[i:i,42:82])
  
  vec1 = log2(vec1)
  vec2 = log2(vec2)
  
  #  vec1 = vec1 - mean(vec1)
  #  vec2 = vec2 - mean(vec2)
  
  #  vec1 = scale(vec1)
  # vec2 = scale(vec2)
  
  tt = kruskal.test(list(vec1, vec2))
  names <- c(names, row.names(data)[i])
  pvals <- c(pvals, tt$p.value)
}

padjusted = p.adjust(pvals, method = "BH")

resKW = data.frame(row.names=names, "p-value"=pvals, "Adjusted p-value"=padjusted)
save(resKW,file="kwstats.rdata")

kept = padjusted < 0.05

res2 = data.frame("Metabolites"=names[kept], "p-value"=pvals[kept], "Adjusted p-value"=padjusted[kept])

cat(nrow(res2))

write.csv(resKW, file="plasma_kruskal_test_result_all.csv")
write.csv(res2, file="plasma_kruskal_test_result.csv")

kwMeta <- names[kept]
save(kwMeta, resKW, file = "dem_k_wallis.RData")
