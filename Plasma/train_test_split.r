library("caret")

setwd("E:/Academic/Thesis/Serum and Plasma Metabolomic/plasma/")

load("features_class.RData")

train_index <- createDataPartition(class[[1]], p=0.7, list=FALSE)


featuresTrain <- features[train_index, 1:ncol(features)]
classTrain <- class[train_index, 1:ncol(class)]
featuresTest <- features[-train_index, 1:ncol(features)]
classTest <- class[-train_index, 1:ncol(class)]

save(featuresTrain, featuresTest, classTrain, classTest, file="train_test_features_class.RData")