#this script trains and performs classification of GCs cells in SN and S cells using DEG expression
library(caret)

#produced in data_preparation_genes.R script
load("exprs_trainTransformed.Rdata")
load("exprs_testTransformed.Rdata")


label_test <- exprs_testTransformed$label

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

#SVM with Linear Kernel
#different algorithms were tested, based on perfomance stats SVM with Liner Kernel was chosen
set.seed(825)
exprs_svm <- train(label ~ ., data = exprs_trainTransformed, 
                 method = "svmLinear", 
                 trControl = svmControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 metric = "ROC",
                 tuneLength = 15)

svm_test_prob <- predict(exprs_svm, newdata = head(exprs_testTransformed), type = "prob")
svm_test <- predict(exprs_svm, newdata = head(exprs_testTransformed))

confusionMatrix(data = svm_test, reference = label_test)
postResample(pred = svm_test, obs = label_test)

#MLP (multi-layer percepton)
set.seed(825)
exprs_mlp <- train(label ~ ., data = exprs_trainTransformed, 
                 method = "mlp", 
                 trControl = svmControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 metric = "ROC",
                 tuneLength = 15)

mlp_test <- predict(exprs_mlp, newdata = head(exprs_testTransformed))
confusionMatrix(data = mlp_test, reference = label_test)
postResample(pred = mlp_test, obs = label_test)


perf_mlp <- data.frame(method=c("mlp"),
                       roc=c(exprs_mlp$resample$ROC),
                       sens=c(exprs_mlp$resample$Sens),
                       spec=c(exprs_mlp$resample$Spec))

perf_svm <- data.frame(method=c("svm"),
                       roc=c(exprs_svm$resample$ROC),
                       sens=c(exprs_svm$resample$Sens),
                       spec=c(exprs_svm$resample$Spec))

perf <- rbind(perf_mlp, perf_svm)

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

test_val =data.frame(method=c("mlp","svm"),
                     perf=c(unname(postResample(pred = mlp_test, obs = label_test)[1]),
                            unname(postResample(pred = svm_test, obs = label_test)[1])))
ggplot(data=test_val, aes(x=method, y=perf)) +
  geom_bar(stat="identity")+theme_classic()+ylab("Accuracy")+xlab("")+
  theme(axis.text.x=element_text(size=20, angle = 90), axis.title=element_text(size=20),
        plot.title = element_text(size=10), axis.text.y=element_text(size=20))

save(exprs_svm, file="exprs_svm.Rdata")
save(exprs_mlp, file="exprs_mlp.Rdata")

###################
#genetic algorithm

#SVM with linear kernel
svm_ga_ctrl <- gafsControl(functions = caretGA, method = "cv", number = 10)
set.seed(825)

svm_ga_search <- gafs(
  x = exprs_trainTransformed[,1:47], 
  y = exprs_trainTransformed$label,
  iters = 20, 
  gafsControl = svm_ga_ctrl,
  # now options to `train` for caretGA
  method = "svmLinear",
  trControl = trainControl(method = "cv", allowParallel = FALSE)
) 

label <- exprs_trainTransformed$label
exprs_trainTransformed_2 <- exprs_trainTransformed[,colnames(exprs_trainTransformed)%in%svm_ga_search$optVariables]
exprs_trainTransformed_2$label <- label

set.seed(825)
exprs_svm_2 <- train(label ~ ., data = exprs_trainTransformed_2, 
                   method = "svmLinear", 
                   trControl = svmControl,
                   ## This last option is actually one
                   ## for gbm() that passes through
                   verbose = FALSE,
                   metric = "ROC",
                   tuneLength = 15)

svm_test_2 <- predict(exprs_svm_2, newdata = head(exprs_testTransformed))

confusionMatrix(data = svm_test_2, reference = label_test)
postResample(pred = svm_test_2, obs = label_test)

#mlp

set.seed(825)

mlp_ga_search <- gafs(
  x = exprs_trainTransformed[,1:15], 
  y = exprs_trainTransformed$label,
  iters = 10, 
  gafsControl = svm_ga_ctrl,
  # now options to `train` for caretGA
  method = "mlp",
  trControl = trainControl(method = "cv", allowParallel = FALSE)
) 

exprs_trainTransformed_2 <- exprs_trainTransformed[,colnames(exprs_trainTransformed)%in%mlp_ga_search$optVariables]
exprs_trainTransformed_2$label <- label

set.seed(825)
exprs_mlp_2 <- train(label ~ ., data = exprs_trainTransformed_2, 
                   method = "mlp", 
                   trControl = svmControl,
                   ## This last option is actually one
                   ## for gbm() that passes through
                   verbose = FALSE,
                   metric = "ROC",
                   tuneLength = 15)

mlp_test_2 <- predict(exprs_mlp_2, newdata = head(exprs_testTransformed))

confusionMatrix(data = mlp_test_2, reference = label_test)
postResample(pred = mlp_test_2, obs = label_test)
save(exprs_svm_2, file="exprs_svm_2.Rdata")
save(exprs_mlp_2, file="exprs_mlp_2.Rdata")

