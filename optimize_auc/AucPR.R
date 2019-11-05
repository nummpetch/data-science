rm(list=ls())
cat("\f")
library(compiler)
library(glmnet)
library(pROC)
library(plsgenomics)
require(caret)
#library(preprocessCore)
require(MLmetrics)

# perform nested cross-validation
cvcv_lda <- function(X, y, nfolds, alpha){
  #Output: AUC and Cross Entropy
  nvar = rep(0, nfolds)                                  
  auc = rep(0, nfolds)                                   
  acc = rep(0, nfolds)                                   
  ce = rep(0, nfolds)                                    
  set.seed(ntime)                                        
  cvSplits <- createFolds(y, nfolds, returnTrain = TRUE) 
  grid=10^seq(10,-2,length=100)                          
  for(i in 1:nfolds){
    Xtrain = X[cvSplits[[i]],]                           # create training input
    ytrain = y[cvSplits[[i]]]                            # create training output
    Xpos = Xtrain[which(ytrain == 1),]                   # create training positive class input 
    Xneg = Xtrain[which(ytrain == 0),]                   # create training negative class input
    mu = colMeans(Xpos) - colMeans(Xneg)                 # calculate pairwise mean
    sigma = cov(Xpos) + cov(Xneg)                        # claculate covariance
    # perform cross-validate regularized linear regression without intercept
    cvob = cv.glmnet(sigma,mu,standardize=FALSE,family="gaussian",intercept=FALSE,alpha=alpha,nfolds=nfolds)
    pred = predict(cvob,X[-cvSplits[[i]],],s="lambda.1se",standardize=FALSE,type="response")
    
    auc[i] = auc(y[-cvSplits[[i]]], as.vector(pred))     
    acc[i] = mean((round(pred) == y[-cvSplits[[i]]]))    
    ce[i] = LogLoss(pred, as.numeric(y[-cvSplits[[i]]])) 
    #count variable selection
    if (length(dim(predict(cvob,type="nonzero")))==0){   
      nvar[i] = 0;                                       
    } else {                                             
      nvar[i] = nrow(predict(cvob,type="nonzero"))       
    }                                                    
  }
  out = c(mean(auc), mean(acc), mean(ce), mean(nvar))    
  return(out)                                            
}
cvcv_lda = cmpfun(cvcv_lda)


f_lda <- function(X, y, nfolds){
  #Input: Matrix X[1..nsample, 1..d], column vector y[1..nsample] and nfolds
  #Output: AUC and Cross Entropy
  auc = rep(0, nfolds)                                   
  acc = rep(0, nfolds)                                   
  ce = rep(0, nfolds)                                   
  set.seed(ntime)                                        
  cvSplits <- createFolds(y, nfolds, returnTrain = TRUE) 
  grid=10^seq(10,-3,length=100)                          
  for(i in 1:nfolds){
    result = cbind(X[cvSplits[[i]],], y[cvSplits[[i]]])
    dat <- as.data.frame(result)
    pre <- pls.lda(X[cvSplits[[i]],],y[cvSplits[[i]]],X[-cvSplits[[i]],],ncomp=3,nruncv=nfolds)
    pred <- pre$predclass
    pred <- as.character(pred)
    pred <- as.numeric(pred)
    
    auc[i] = auc(y[-cvSplits[[i]]], pred)     
    acc[i] = mean((round(pred) == y[-cvSplits[[i]]]))    
    ce[i] = LogLoss(pred, as.numeric(y[-cvSplits[[i]]])) 
                                             
  }
  out = c(mean(auc), mean(acc), mean(ce))    
  return(out)                                        
}
f_lda = cmpfun(f_lda)


#print("colon")
#load("colon.rda")
#load("C:/Users/teer/Downloads/colon.rda")
#load("~/Dropbox/3 work with kampol/Rcode/colon.rda")
#X = colon.x
#y = colon.y
#y = as.character(y)
#print("leuke")
#load("leukemia.rda")
#load("C:/Users/teer/Downloads/leukemia.rda")
#X = lx.original
#y = ly.original
#y = as.character(y)
print("prostate")
load("C:/Users/teer/Downloads/prostate.rda")
#load("~/Dropbox/3 work with kampol/Rcode/prostate.rda")
#load("prostate.rda")
X = prostate.x
y = prostate.y
y = as.character(y)
a <- table(y)
print(a)
#print("lym")
#load("lymphoma.rda")
#load("~/Dropbox/3 work with kampol/Rcode/lymphoma.rda")
#load("C:/Users/teer/Downloads/lymphoma.rda")
#X = lymphoma.x
#y = lymphoma.y
#y[c(which(y != 2))]=0
#y[c(which(y == 2))]=1
#y = as.character(y)
#load("~/Dropbox/3 work with kampol/Rcode/srbct.rda")
#X = srbct.x
#y = srbct.y
#y[c(which(y != 3))]=0
#y[c(which(y == 3))]=1
#y = as.character(y)
#load("brain.rda")
#load("C:/Users/teer/Downloads/brain.rda")
#X = brain.x
#y = brain.y
#y[c(which(y != 2))]=0
#y[c(which(y == 2))]=1
#y = as.character(y)
#print("breast")
#load("west.RData")
#load("C:/Users/teer/Downloads/west.RData")
#load("~/Dropbox/3 work with kampol/Rcode/west.RData")
#X = west$x
#X = as.matrix(X)
#y = west$y
#y = as.character(y)
#y[which(y=="negative")] <- "0"
#y[which(y=="positive")] <- "1"

#load("shipp.RData")
#load("C:/Users/teer/Downloads/shipp.RData")
#load("~/Dropbox/3 work with kampol/Rcode/shipp.RData")
#X = shipp$x
#X = as.matrix(X)
#y = shipp$y
#y = as.character(y)
#y[which(y=="DLBCL")] <- "0"
#y[which(y=="FL")] <- "1"


nfolds = 10   # number of fold in cross-validation
ntimes = 30 # number of repeat experiment
auc = Matrix(data=NA, nrow=ntimes, ncol=3)
result <- c()
#alpha = 1
for (i in 1:ntimes){
  ntime = i
  print(i)
  auc = f_lda(X, y, nfolds)
  result <- rbind(result,auc)
}
write.csv(result, file = "mydata2.csv",row.names=FALSE)
