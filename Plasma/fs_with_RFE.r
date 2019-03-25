rm(list=ls())

library("e1071")

setwd("E:/Academic/Thesis/Serum and Plasma Metabolomic/plasma/")

load("dem_t_test.RData")
load("dem_k_wallis.RData")
load("dem_wilcoxon.RData")

DEMetabols = union(ttMeta, kwMeta)
DEMetabols = union(DEMetabols, wilMeta)

load("features_class.RData")

features = features[,DEMetabols]

# ensure the results are repeatable
set.seed(7)
# load the library
library(mlbench)
library(caret)
# calculate correlation matrix
correlationMatrix <- cor(features)
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75, names=TRUE)
# print indexes of highly correlated attributes
print(highlyCorrelated)

# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
system.time(results <- rfe(features[,1:28], as.factor(class[,1]), sizes=c(1:10), rfeControl=control))
# summarize the results
print(results)
# list the chosen features
rfeRanked = predictors(results)
# plot the results
plot(results, type=c("g", "o"))

bestMetabols = rfeRanked[1:8]

#write.csv(bestMetabols, file="rfePlasmaMetabols.csv")

#save(bestMetabols, file="rfePlasmaMetabols.rdata")

load("rfePlasmaMetabols.rdata")

#library(Boruta)

#for (i in 1:41)
#  class[i,1] = 1
#for (i in 41:82)
#  class[i,1] = 0

#boruta.result <- Boruta(as.matrix(features[1:82,1:28]), as.factor(class[,1]), doTrace = 2, maxRuns=1000)

#print(boruta.result)

#final.boruta = TentativeRoughFix(boruta.result)

#print(final.boruta)

#Metabols = getSelectedAttributes(boruta.result, withTentative = F)

load("train_test_features_class.RData")
load("features_class.RData")

features = features[,bestMetabols]
featuresTrain = featuresTrain[,bestMetabols]
featuresTest = featuresTest[,bestMetabols]
classTrain <- as.factor(classTrain)

svm_model <- svm(classTrain ~ ., data = featuresTrain, kernel="radial")
summary(svm_model)

system.time(predicted <- predict(svm_model, featuresTest))

mat <- table(predicted, classTest)

accuracy1 <- sum(diag(mat)) / sum(mat)

#system.time(svm_tune <- tune(svm, train.x=featuresTrain, train.y=as.factor(classTrain), 
#                kernel="radial", ranges=list(cost=2^(-16:16), gamma=c(2^(-16:16)))))

#print(svm_tune)

#svm_model_after_tune <- svm(classTrain ~ ., data = featuresTrain, kernel="radial", cost=svm_tune$best.parameters$cost, gamma=svm_tune$best.parameters$gamma)
svm_model_after_tune <- svm(classTrain ~ ., data = featuresTrain, kernel="radial", cost=4, gamma=0.00390625)
#summary(svm_model_after_tune)

system.time(predicted <- predict(svm_model_after_tune,featuresTest, decision.values = TRUE))

decision_values = attr(predicted,"decision.values")

mat = table(predicted, classTest)

accuracy2 <- sum(diag(mat)) / sum(mat)

print(accuracy1)
print(accuracy2)

sen = mat[2,2]/(mat[2,2]+mat[1,2])
spe = mat[1,1]/(mat[1,1]+mat[2,1])
print(sen)
print(spe)

# k fold

# predictedClass = matrix("", nrow = nrow(features), ncol = 1);
# predictedFvalue = matrix(0L, nrow = nrow(features), ncol = 1);
# 
# 
# numFold = 82
# foldSize = as.integer(nrow(features)/numFold)
# 
# for (k in  1:(numFold+1)){
#   train_index = vector("numeric")
#   test_index = vector("numeric")
#   for (i in 1:nrow(features))
#   {
#     foldno = as.integer((i-1)/foldSize) + 1
# 
#     if (foldno == k)
#     {
#       test_index = c(test_index, i)
#     }
#     else
#     {
#       train_index = c(train_index, i)
#     }
#   }
#   tempTrainClass = class[train_index,1:1];
#   tempTrainClass = as.factor(tempTrainClass)
#   tempTrainData = features[train_index, 1:ncol(features)];
#   tempTestData=features[-train_index, 1:ncol(features)];
# 
#   svm_model <- svm(tempTrainClass ~ ., data = tempTrainData, kernel="radial", cost=4, gamma=0.00390625)
#   predicted <- predict(svm_model, tempTestData, decision.values = TRUE)
# 
#   decision_values2 = attr(predicted,"decision.values")
# 
#   for (i in 1:length(test_index))
#   {
#     ti = test_index[i]
#     predictedClass[ti, 1] = as.character(predicted[i]);
#     predictedFvalue[ti, 1] = decision_values2[i];
#   }
# }
# 
# predClass= vector(mode = "numeric", length = nrow(features));
# predFval = vector(mode = "numeric", length = nrow(features));
# 
# for (i in 1:length(predictedClass))
# {
#   if (i > 41)
#       {predClass[i] = 1}
#   predFval[i] = predictedFvalue[i]
# }

library(pROC)

roc1 <- roc(classTest,
            decision_values[,1],
            # arguments for ci
            ci=TRUE, boot.n=10000,
            # arguments for plot
            plot=TRUE, smooth=TRUE, auc.polygon=TRUE, max.auc.polygon=FALSE, grid=TRUE,
            print.auc=TRUE, show.thres=TRUE)

# roc1 <- roc(predClass,
#             predFval,
#             # arguments for ci
#             ci=TRUE, boot.n=10000,
#             # arguments for plot
#             plot=TRUE, smooth=TRUE, auc.polygon=TRUE, max.auc.polygon=FALSE, grid=TRUE,
#             print.auc=TRUE, show.thres=TRUE)



# library(ROCR)
# 
# #pred_obj = prediction(predictedFvalue, class)
# pred_obj = prediction(decision_values, classTest)
# 
# 
# # Calculating Area under Curve
# perf_val <- performance(pred_obj,"auc")
# perf_val
# 
# # Calculating True Positive and False Positive Rate
# perf_val <- performance(pred_obj, "tpr", "fpr")
# 
# # Plot the ROC curve
# plot(perf_val, col = "green", lwd = 1.5)
