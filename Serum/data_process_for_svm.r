rm(list=ls())

setwd("E:/Academic/4-1/Thesis/Serum and Plasma Metabolomic/serum/")

load("dem_t_test.RData")
load("dem_k_wallis.RData")
load("dem_wilcoxon.RData")

DEMetabols = union(ttMeta, kwMeta)
DEMetabols = union(DEMetabols, wilMeta)

data <- read.csv("serumLC.csv", row.names=1)

names <- vector('character')
class_num <- vector('numeric')

features <- data.frame(matrix(nrow=82, ncol=nrow(data)))
class <- data.frame(matrix(nrow=82, ncol=1))

df <- data.frame(matrix(nrow=82, ncol=nrow(data)+2))

for (i in 1:nrow(data))
{
  names[i] = row.names(data)[i]
}

colnames(features) <- names
colnames(class) <- "Class"

colnames(df) <- c("Samples", "Label", names)

for (i in 1:ncol(data))
{
  features[i, ] <- data[1:nrow(data), i]
  
  if (i <= 41)
  {
    class[i, ] <- "disease"
    class_num = c(class_num, 0)
    df[i, ] <- c(colnames(data)[i], 0, data[1:nrow(data), i])
  }
  else 
  {
    class[i, ] <- "control"
    class_num = c(class_num, 1)
    df[i, ] <- c(colnames(data)[i], 1, data[1:nrow(data), i])
  }
}

df = df[,c("Samples","Label",DEMetabols)]

write.csv(df, file="data_for_ROC_MetaboAnalyst.csv")

save(features, class, file="features_class.RData")