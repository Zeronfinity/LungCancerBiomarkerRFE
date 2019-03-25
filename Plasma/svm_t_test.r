library("e1071")
rm(list=ls())

setwd("E:/Academic/4-1/Thesis/Serum and Plasma Metabolomic/plasma/")

load("train_test_features_class.RData")
load("dem_t_test.RData")

featuresTrain = featuresTrain[,ttMeta]
featuresTest = featuresTest[,ttMeta]
classTrain <- as.factor(classTrain)

svm_model <- svm(classTrain ~ ., data = featuresTrain, kernel="radial")
summary(svm_model)

system.time(predicted <- predict(svm_model, featuresTest))

mat <- table(predicted, classTest)

accuracy1 <- sum(diag(mat)) / sum(mat)

system.time(svm_tune <- tune(svm, train.x=featuresTrain, train.y=as.factor(classTrain), 
                             kernel="radial", ranges=list(cost=2^(-1:3), gamma=c(2^(-7:0)))))

#print(svm_tune)


svm_model_after_tune <- svm(classTrain ~ ., data = featuresTrain, kernel="radial", cost=1, gamma=0.01)
#summary(svm_model_after_tune)

system.time(predicted <- predict(svm_model_after_tune,featuresTest))

mat = table(predicted, classTest)

accuracy2 <- sum(diag(mat)) / sum(mat)

print(accuracy1)
print(accuracy2)