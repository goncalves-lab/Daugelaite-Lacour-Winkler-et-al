#this script trains classification of GCs cells in SN and S cells using scenic AUC scores

library(caret)
library(ggplot2)
#input data are produced in the script data_preparation_tfs.R 
load("auc_trainTransformed.Rdata")
load("auc_testTransformed.Rdata")


label_test <- auc_testTransformed$label

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=T,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#resampling functions
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

svmControl <- trainControl(method = "repeatedcv",
                           number = 10, repeats = 10,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           search = "random")

###############
#model fitting
###############
#different algorithms were tested, based on perfomance stats SVM with Linear Kernel was chosen
#SVM with Linear Kernel
set.seed(825)
auc_svm <- train(label ~ ., data = auc_trainTransformed, 
                 method = "svmLinear", 
                 trControl = svmControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 metric = "ROC",
                 tuneLength = 15)

svm_test_prob <- predict(auc_svm, newdata = auc_testTransformed, type = "prob")
svm_test <- predict(auc_svm, newdata = auc_testTransformed)

confusionMatrix(data = svm_test, reference = label_test)
postResample(pred = svm_test, obs = label_test)

#L2 Regularized Support Vector Machine (dual) with Linear Kernel
set.seed(825)
auc_svm3 <- train(label ~ ., data = auc_trainTransformed, 
                 method = "svmLinear3", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 metric = "Accuracy",
                 tuneLength = 15)

svml2_test <- predict(auc_svm3, newdata = head(auc_testTransformed))
confusionMatrix(data = svml2_test, reference = label_test)
postResample(pred = svm_test, obs = label_test)

#Extreme boosting
set.seed(825)
auc_xgb <- train(label ~ ., data = auc_trainTransformed, 
                  method = "xgbTree", 
                  trControl = svmControl,
                  ## This last option is actually one
                  ## for gbm() that passes through
                  verbose = FALSE,
                  metric = "ROC",
                  tuneLength = 15)

xgb_test <- predict(auc_xgb, newdata = head(auc_testTransformed))
confusionMatrix(data = xgb_test, reference = label_test)
postResample(pred = xgb_test, obs = label_test)

#Random forest
set.seed(825)
auc_rf <- train(label ~ ., data = auc_trainTransformed, 
                 method = "rf", 
                 trControl = svmControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 metric = "ROC",
                 tuneLength = 15)

rf_test <- predict(auc_rf, newdata = head(auc_testTransformed))
confusionMatrix(data = rf_test, reference = label_test)
postResample(pred = rf_test, obs = label_test)

#MLP (multi-layer percepton)
set.seed(825)
auc_mlp <- train(label ~ ., data = auc_trainTransformed, 
                method = "mlp", 
                trControl = svmControl,
                ## This last option is actually one
                ## for gbm() that passes through
                verbose = FALSE,
                metric = "ROC",
                tuneLength = 15)

mlp_test <- predict(auc_mlp, newdata = auc_testTransformed)
confusionMatrix(data = mlp_test, reference = label_test)
postResample(pred = mlp_test, obs = label_test)

#ELM (extreme learning machine)
set.seed(825)
auc_elm <- train(label ~ ., data = auc_trainTransformed, 
                 method = "elm", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 metric = "Accuracy",
                 tuneLength = 15)

elm_test <- predict(auc_elm, newdata = head(auc_testTransformed))
confusionMatrix(data = elm_test, reference = label_test)
postResample(pred = elm_test, obs = label_test)

#Naive Bayes

set.seed(825)
auc_nb <- train(label ~ ., data = auc_trainTransformed, 
                 method = "nb", 
                 trControl = svmControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 metric = "ROC",
                 tuneLength = 15)

nb_test <- predict(auc_nb, newdata = head(auc_testTransformed))
confusionMatrix(data = nb_test, reference = label_test)
postResample(pred = mlp_test, obs = label_test)

#LogitBoost (logistic regression)

set.seed(825)
auc_lb <- train(label ~ ., data = auc_trainTransformed, 
                method = "LogitBoost", 
                trControl = svmControl,
                ## This last option is actually one
                ## for gbm() that passes through
                verbose = FALSE,
                metric = "ROC",
                tuneLength = 15)

lb_test <- predict(auc_lb, newdata = head(auc_testTransformed))
confusionMatrix(data = lb_test, reference = label_test)
postResample(pred = lb_test, obs = label_test)

#LDA

set.seed(825)
auc_lda <- train(label ~ ., data = auc_trainTransformed, 
                method = "lda", 
                trControl = svmControl,
                ## This last option is actually one
                ## for gbm() that passes through
                verbose = FALSE,
                metric = "ROC",
                tuneLength = 15)

lda_test <- predict(auc_lda, newdata = head(auc_testTransformed))
confusionMatrix(data = lda_test, reference = label_test)
postResample(pred = lda_test, obs = label_test)

###################
#genetic algorithm

#SVM with linear kernel
svm_ga_ctrl <- gafsControl(functions = caretGA, method = "cv", number = 10)
set.seed(825)

svm_ga_search <- gafs(
  x = auc_trainTransformed[,-ncol(auc_trainTransformed)], 
  y = auc_trainTransformed$label,
  iters = 15, 
  gafsControl = svm_ga_ctrl,
  # now options to `train` for caretGA
  method = "svmLinear",
  trControl = trainControl(method = "cv", allowParallel = FALSE)
) 

label <- auc_trainTransformed$label
auc_trainTransformed_2 <- auc_trainTransformed[,colnames(auc_trainTransformed)%in%svm_ga_search$optVariables]
auc_trainTransformed_2$label <- label

set.seed(825)
auc_svm_2 <- train(label ~ ., data = auc_trainTransformed_2, 
                 method = "svmLinear", 
                 trControl = svmControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 metric = "ROC",
                 tuneLength = 15)

svm_test_2 <- predict(auc_svm_2, newdata = auc_testTransformed)

confusionMatrix(data = svm_test_2, reference = label_test)
postResample(pred = svm_test_2, obs = label_test)

#mlp

set.seed(825)

mlp_ga_search <- gafs(
  x = auc_trainTransformed[,-ncol(auc_trainTransformed)], 
  y = auc_trainTransformed$label,
  iters = 10, 
  gafsControl = svm_ga_ctrl,
  # now options to `train` for caretGA
  method = "mlp",
  trControl = trainControl(method = "cv", allowParallel = FALSE)
) 

auc_trainTransformed_2 <- auc_trainTransformed[,colnames(auc_trainTransformed)%in%mlp_ga_search$optVariables]
auc_trainTransformed_2$label <- label

set.seed(825)
auc_mlp_2 <- train(label ~ ., data = auc_trainTransformed_2, 
                   method = "mlp", 
                   trControl = svmControl,
                   ## This last option is actually one
                   ## for gbm() that passes through
                   verbose = FALSE,
                   metric = "ROC",
                   tuneLength = 15)

mlp_test_2 <- predict(auc_mlp_2, newdata = auc_testTransformed)

confusionMatrix(data = mlp_test_2, reference = label_test)
postResample(pred = mlp_test_2, obs = label_test)


perf_xgb <- data.frame(method=c("xgb"),
                  roc=c(auc_xgb$resample$ROC),
                  sens=c(auc_xgb$resample$Sens),
                  spec=c(auc_xgb$resample$Spec))

perf_rf <- data.frame(method=c("rf"),
                       roc=c(auc_rf$resample$ROC),
                       sens=c(auc_rf$resample$Sens),
                       spec=c(auc_rf$resample$Spec))
perf_mlp <- data.frame(method=c("mlp"),
                       roc=c(auc_mlp$resample$ROC),
                       sens=c(auc_mlp$resample$Sens),
                       spec=c(auc_mlp$resample$Spec))
perf_nb <- data.frame(method=c("nb"),
                       roc=c(auc_nb$resample$ROC),
                       sens=c(auc_nb$resample$Sens),
                       spec=c(auc_nb$resample$Spec))
perf_lr <- data.frame(method=c("lr"),
                       roc=c(auc_lb$resample$ROC),
                       sens=c(auc_lb$resample$Sens),
                       spec=c(auc_lb$resample$Spec))
perf_lda <- data.frame(method=c("lda"),
                       roc=c(auc_lda$resample$ROC),
                       sens=c(auc_lda$resample$Sens),
                       spec=c(auc_lda$resample$Spec))
perf_svm <- data.frame(method=c("svm"),
                       roc=c(auc_svm$resample$ROC),
                       sens=c(auc_svm$resample$Sens),
                       spec=c(auc_svm$resample$Spec))
perf <- rbind(perf_lda, perf_lr, perf_mlp, perf_nb, perf_rf, perf_svm, perf_xgb)

tgc_spec <- summarySE(perf, measurevar="spec", groupvars=c("method"))

ggplot(tgc_spec, aes(x=method, y=spec)) + 
  geom_errorbar(aes(ymin=spec-se, ymax=spec+se), width=.1) +
  geom_point()+ylim(0,1)+theme_classic()+ylab("Specificity")+xlab("")+
  theme(axis.text.x=element_text(size=20, angle = 90), axis.title=element_text(size=20),
        plot.title = element_text(size=10), axis.text.y=element_text(size=20))

tgc_sens <- summarySE(perf, measurevar="sens", groupvars=c("method"))

ggplot(tgc_sens, aes(x=method, y=sens)) + 
  geom_errorbar(aes(ymin=sens-se, ymax=sens+se), width=.1) +
  geom_point()+ylim(0,1)+theme_classic()+ylab("Sensitivity")+xlab("")+
  theme(axis.text.x=element_text(size=20, angle = 90), axis.title=element_text(size=20),
        plot.title = element_text(size=10), axis.text.y=element_text(size=20))


tgc_roc <- summarySE(perf, measurevar="roc", groupvars=c("method"))

ggplot(tgc_roc, aes(x=method, y=roc)) + 
  geom_errorbar(aes(ymin=roc-se, ymax=roc+se), width=.1) +
  geom_point()+ylim(0,1)+theme_classic()+ylab("ROC")+xlab("")+
  theme(axis.text.x=element_text(size=20, angle = 90), axis.title=element_text(size=20),
        plot.title = element_text(size=10), axis.text.y=element_text(size=20))

test_val =data.frame(method=c("xgb", "rf", "mlp", "nb", "lr","lda", "svm"),
                     perf=c(unname(postResample(pred = xgb_test, obs = label_test)[1]),
                            unname(postResample(pred = rf_test, obs = label_test)[1]),
                            unname(postResample(pred = mlp_test, obs = label_test)[1]),
                            unname(postResample(pred = nb_test, obs = label_test)[1]),
                            unname(postResample(pred = lb_test, obs = label_test)[1]),
                            unname(postResample(pred = lda_test, obs = label_test)[1]),
                            unname(postResample(pred = svm_test, obs = label_test)[1])))

ggplot(data=test_val, aes(x=method, y=perf)) +
  geom_bar(stat="identity")+theme_classic()+ylab("Accuracy")+xlab("")+
  theme(axis.text.x=element_text(size=20, angle = 90), axis.title=element_text(size=20),
        plot.title = element_text(size=10), axis.text.y=element_text(size=20))

save(auc_svm, file="auc_svm.Rdata")
save(auc_svm_2, file="auc_svm_2.Rdata")
save(auc_mlp, file="auc_mlp.Rdata")
save(auc_mlp_2, file="auc_mlp_2.Rdata")
