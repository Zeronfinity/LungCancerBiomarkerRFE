setwd("E:/Academic/4-1/Thesis/Serum and Plasma Metabolomic/serum/")
rm(list=ls())

library(gplots)
library(viridis)

data <- read.csv("serumLC.csv", row.names=1)

load("dem_t_test.RData")
load("dem_k_wallis.RData")
load("dem_wilcoxon.RData")

DEMetabols = union(ttMeta, kwMeta)
DEMetabols = union(DEMetabols, wilMeta)

data = data[DEMetabols,]

names_up <- vector('character')
fc_up <- vector('numeric')
lfc_up <- vector('numeric')
names_down <- vector('character')
fc_down <- vector('numeric')
lfc_down <- vector('numeric')

data_ = as.matrix(data)

for (i in 1:nrow(data)){
  vec1 = as.numeric(data[i:i,1:41])
  vec2 = as.numeric(data[i:i,42:82])
  
  m1 = mean(vec1)
  m2 = mean(vec2)

  data_[i,] = (data_[i,]-mean(data_[i,]))/(max(data_[i,]) - min(data_[i,]))

  FC = m1/m2
  
  logFC = log2(FC)
  
  if (logFC > 0)
  {
    names_up <- c(names_up, row.names(data)[i])
    fc_up <- c(fc_up, FC)
    lfc_up <- c(lfc_up, logFC)    
  }
  else
  {
    names_down <- c(names_down, row.names(data)[i])
    fc_down <- c(fc_down, FC)
    lfc_down <- c(lfc_down, logFC)    
  }
}


resUp = data.frame("Metabolites"=names_up, "FC Value"=fc_up, "log2 of FC"=lfc_up)
resDown = data.frame("Metabolites"=names_down, "FC Value"=fc_down, "log2 of FC"=lfc_down)

concol = vector('character')

for (j in 1:ncol(data_))
{
  if (j <= 41)
    concol = c(concol, '#0F6C01')
  else 
    concol = c(concol, '#21B80A')  
}  

#heatmap.2(data_, main="Significant Metabolites of Serum Sample", trace="none", margins=c(4,14), density="none", 
 #         cexRow=1, cexCol=0.2, col=magma(100), ColSideColors=concol)
#legend(0.82,1.22,legend=c("Cancer","Control"),fill=c('#0F6C01','#21B80A'),cex=0.7)

write.csv(resUp, file="upregulated.csv")
write.csv(resDown, file="downregulated.csv")
