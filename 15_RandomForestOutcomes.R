########################################
#' Try random forest to predict group
#' Karine Thai - 2023/08/19
########################################

######################
#' 0 - Load required packages
######################

library(randomForest)
library(dplyr)
library(caret)

######################
#' 1 - Load table of MFIs and proportions with clinical data
######################

load('Output/RData/PropMFIMetaData.RData')
load('Output/RData/PropTableAll.RData')

######################
#' 3 - Subset Outcome 5 yrs group
######################

PropTableOutcome5yrs <- PropTableAll %>%
  filter(Group == "MS" & Cohort == "Outcome5yrs")

PropMFIOutcome5yrs <- PropMFIMetaData %>%
  filter(Group == "MS" & Cohort == "Outcome5yrs")

######################
#' 4 - Randomly subset to have equal number or HC and MS
######################

# Subset the table for HC and MS separately
subset_EDA <- PropTableOutcome5yrs %>% filter(Progression_5yrs == "Progressor")
subset_NEDA <- PropTableOutcome5yrs %>% filter(Progression_5yrs == "Stable")

subset_EDA <- PropMFIOutcome5yrs %>% filter(NEDA5yrs == "EDA")
subset_NEDA <- PropMFIOutcome5yrs %>% filter(NEDA5yrs == "NEDA")

# Randomly sample an equal number of rows from both subsets

num_samples <- min(nrow(subset_EDA), nrow(subset_NEDA))
set.seed(1)
subset_EDA <- subset_EDA %>% sample_n(num_samples)
set.seed(1)
subset_NEDA <- subset_NEDA %>% sample_n(num_samples)

# Combine the subsets back together
subPropTableOutcome5yrs <- bind_rows(subset_EDA, subset_NEDA)

subPropMFIOutcome5yrs <- bind_rows(subset_EDA, subset_NEDA)

######################
#' 5 - Change columns to factors
######################

Factor_col=c('Cohort','Group','Sex','MSTypeSampling','Smoking','TreatmentSampling','Progression_5yrs','ClinicalActivity5yrs', 'RadiologicalActivity5yrs', 'NEDA5yrs', 'Severity5yrs')
for (fac in Factor_col){
  PropTableOutcome5yrs[,fac]=as.factor(PropTableOutcome5yrs[,fac])
}

for (fac in Factor_col){
  subPropMFIOutcome5yrs[,fac]=as.factor(subPropMFIOutcome5yrs[,fac])
}

######################
#' 6 - Compute random forest
######################

set.seed(1)
rf = randomForest(NEDA5yrs~., data=subPropMFIOutcome5yrs[,c(1:170,173:176,210)], ntree = 2500, na.action = na.omit,  proximity=T)
# rf = randomForest(Group~., data=subPropTableUntreatedAllHC[,c(1:35, 37:39)], ntree = 500, na.action = na.omit, proximity=T)

rf
#Determine ntrees
oob.error.data = data.frame(
  Trees=rep(1:nrow(rf$err.rate), times=3),
  Type=rep(c('OOB', 'NEDA', 'EDA'), each=nrow(rf$err.rate)),
  Error=c(rf$err.rate[,'OOB'],
          rf$err.rate[,'NEDA'],
          rf$err.rate[,'EDA'])
)
png(file = "Output/Plots/nTree_RandomForestGroupMFIandProp.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
ggplot(data=oob.error.data, aes(x=Trees,y=Error))+geom_line(aes(color=Type))
dev.off()
#Determine mtry
oob.values = vector(length=50)
set.seed(1)
for (i in 1:50) {
  temp.model = randomForest(NEDA5yrs~., data=subPropMFIOutcome5yrs[,c(1:170,173:176,210)], ntree = 2250, na.action = na.omit, mtry=i, proximity=T)
  oob.values[i] = temp.model$err.rate[nrow(temp.model$err.rate),1]
}
set.seed(1)
for (i in 1:20) {
  temp.model = randomForest(Group~., data=subPropMFIOutcome5yrs[,c(1:170,173:176,181,182,186,187,192,210)], ntree=1500,mtry=i)
  oob.values[i] = temp.model$err.rate[nrow(temp.model$err.rate),1]
}
png(file = "Output/Plots/mtry_RandomForestGroupMFIandProp.png",
    width = 8,
    height = 6,
    units = 'in',
    res = 600)
plot(oob.values, xlab = 'mtry', xaxp = c(1,20,19))
dev.off()

#Compute RF
set.seed(1)
rf = randomForest(NEDA5yrs~., data=subPropMFIOutcome5yrs[,c(1:170,173:176,210)], ntree = 2250, mtry = 30, proximity=T)
# rf = randomForest(Group~., data=subPropTableUntreatedAllHC[,c(1:35, 37:39)], ntree = 500, mtry = 6, na.action = na.omit, proximity=T)

rf

#Look at most important variable
varImpPlot(rf)

#Make mds plot
distance.matrix = dist(1-rf$proximity)
mds.stuff = cmdscale(distance.matrix, eig = T, x.ret = T)
mds.var.per = round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
mds.values = mds.stuff$points
mds.data = data.frame(Sample=rownames(mds.values),
                      X=mds.values[,1],
                      Y=mds.values[,2],
                      Status=subPropMFIOutcome5yrs$NEDA5yrs)
png(file = "Output/Plots/MDSplot_RandomForestGroupMFIandProp.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) + geom_point(aes(color=Status)) + theme_bw() + theme(text = element_text(size = 20))
dev.off()

#Look at most important variable
png(file = "Output/Plots/VarImp_RandomForestGroupMFIandProp.png",
    width = 10,
    height = 10,
    units = 'in',
    res = 600)
varImpPlot(rf)
dev.off()

######################
#' 7 - Training and validation dataset
######################
set.seed(1)
ind <- sample(2, nrow(subPropMFIUntreatedAllHC), replace = TRUE, prob = c(0.7, 0.3))
train <- subPropMFIUntreatedAllHC[ind==1,]
test <- subPropMFIUntreatedAllHC[ind==2,]
train$Group = factor(train$Group, levels = c('MS', 'HC'))
test$Group = factor(test$Group, levels = c('MS', 'HC'))
set.seed(1)
rf = randomForest(Group~., data=train[,c(1:170, 172)], ntree = 3000, na.action = na.omit, proximity=T)

#Determine ntrees
oob.error.data = data.frame(
  Trees=rep(1:nrow(rf$err.rate), times=3),
  Type=rep(c('OOB', 'HC', 'MS'), each=nrow(rf$err.rate)),
  Error=c(rf$err.rate[,'OOB'],
          rf$err.rate[,'HC'],
          rf$err.rate[,'MS'])
)
png(file = "Output/Plots/nTree_RandomForestGroupMFIandPropTrain.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
ggplot(data=oob.error.data, aes(x=Trees,y=Error))+geom_line(aes(color=Type))
dev.off()
#Determine mtry
oob.values = vector(length=20)
set.seed(1)
for (i in 1:20) {
  temp.model = randomForest(Group~., data=train[,c(1:170, 172)], ntree=2000,mtry=i)
  oob.values[i] = temp.model$err.rate[nrow(temp.model$err.rate),1]
}

png(file = "Output/Plots/mtry_RandomForestGroupMFIandPropTrain.png",
    width = 8,
    height = 6,
    units = 'in',
    res = 600)
plot(oob.values, xlab = 'mtry', xaxp = c(1,20,19))
dev.off()
set.seed(1)

rf = randomForest(Group~., data=train[,c(1:170, 172)], ntree = 2000, mtry = 19, na.action = na.omit, proximity=T)

p1 <- predict(rf, train)
confusionMatrix(p1, train$Group)

p2 = predict(rf, test)
confusionMatrix(p2, test$Group)

#Make mds plot
distance.matrix = dist(1-rf$proximity)
mds.stuff = cmdscale(distance.matrix, eig = T, x.ret = T)
mds.var.per = round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
mds.values = mds.stuff$points
mds.data = data.frame(Sample=rownames(mds.values),
                      X=mds.values[,1],
                      Y=mds.values[,2],
                      Status=train$Group)
png(file = "Output/Plots/MDSplot_RandomForestGroupMFIandPropTrain.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) + geom_point(aes(color=Status)) + theme_bw() + theme(text = element_text(size = 20))
dev.off()

#Logistic regression
model <- glm(NEDA5yrs~., data=subPropMFIOutcome5yrs[,c(1:170, 210)], family = binomial)
summary(model)
# Make predictions
probabilities <- model %>% predict(test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")
# Model accuracy
mean(predicted.classes == test$Group)

classifier = svm(formula = NEDA5yrs ~ .,
                 data = subPropMFIOutcome5yrs[,c(1:170, 210)],
                 type = 'C-classification',
                 kernel = 'radial')
y_pred = predict(classifier, newdata = test[,c(1:170, 172)], type = 'decision')
y_train_pred = predict(classifier, newdata = train[,c(1:170, 172)])
cm = table(test[,172], y_pred)
cm2 = table(train[,172], y_train_pred )
y_pred

correct_pred <- y_pred == test$Group
table(correct_pred)

prop.table(table(correct_pred))

varImp(bayes2)
bayes <- naiveBayes(Group ~ ., data = train[,c(1:170, 172)])
bayes

# Predicting on test data'
bayes_pred <- predict(bayes, newdata = test[,c(1:170,172)])

# Confusion Matrix
cm <- table(test$Group, bayes_pred)
cm

# Model Evaluation
confusionMatrix(cm)
varImp(bayes)


# Fitting KNN Model 
# to training dataset
model = train(train[,c(1:170)],train$Group,'nb',trControl=trainControl(method='cv',number=10))
train_scale = scale(train[,c(1:170)])
test_scale = scale(test[,c(1:170)])
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train$Group,
                      k = 1)
classifier_knn

# Confusiin Matrix
cm <- table(test$Group, classifier_knn)
cm

# Model Evaluation - Choosing K
# Calculate out of Sample error
misClassError <- mean(classifier_knn != test$Group)
print(paste('Accuracy =', 1-misClassError))

# K = 3
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train$Group,
                      k = 3)
misClassError <- mean(classifier_knn != test$Group)
print(paste('Accuracy =', 1-misClassError))

# K = 5
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train$Group,
                      k = 5)
misClassError <- mean(classifier_knn != test$Group)
print(paste('Accuracy =', 1-misClassError))

# K = 7
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train$Group,
                      k = 7)
misClassError <- mean(classifier_knn != test$Group)
print(paste('Accuracy =', 1-misClassError))

# K = 15
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train$Group,
                      k = 15)
misClassError <- mean(classifier_knn != test$Group)
print(paste('Accuracy =', 1-misClassError))

# K = 19
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train$Group,
                      k = 19)
misClassError <- mean(classifier_knn != test$Group)
print(paste('Accuracy =', 1-misClassError))
