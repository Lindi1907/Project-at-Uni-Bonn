library(grf)
library(ranger)
set.seed(1999)
dgp<-function(N,p,case){
  
  X <- matrix(rnorm(N * p),nrow = N,ncol = p)
  tau.true <- pmax(X[, 1], 0)
  
  if(case=="case 1"){
    W <- rbinom(N,1, p = 0.5)
    E <- X[, 2] + pmin(X[, 3], 0)
  }
  
  if(case=="case 2"){
    W<-rbinom(N,1,0.4+0.2*(X[,1]>0))
    E <- X[, 2] + pmin(X[, 3], 0)
  }
  
  if(case=="case 3"){
    W<-rbinom(N,1,0.5)
    E <- pmax(X[,1]+X[,2],X[,3],0)+pmax(X[,4]+X[,5],0)
  }
  
  if(case=="case 4"){
    W<-rbinom(N,1,0.4+0.2*(X[,1]>0))
    E <- pmax(X[,1]+X[,2],X[,3],0)+pmax(X[,4]+X[,5],0)
  }
  
  Y <- tau.true * W + E + rnorm(N)
  df<-data.frame(Y,W,tau.true,X)
  
  return(df)
}

CATE<-function(df){
  t<-df$t
  X<-df[,4:13]
  Y<-df$Y
  W<-df$W
  #Empirical CATE
  Emp.CATE<-mean(t)
  #causal forest CATE
  c.forest <- causal_forest(X, Y, W)
  cf.CATE<-average_treatment_effect(c.forest, target.sample = 'all')
  cf.CATE<-data.frame(cf.CATE)
  cf.CATE<-cf.CATE[1,]
  #random forest CATE
  ranger.model <- ranger(Y ~ ., data = data.frame(X, W, Y))
  ranger.cate <- function(ranger.model, X) {
    data_untreated <- data.frame(X, W = 0)
    data_treated <- data.frame(X, W = 1)
    mean(predict(ranger.model, data_treated)$predictions - predict(ranger.model, data_untreated)$predictions)
  }
  rf.CATE<-ranger.cate(ranger.model, X)
  
  CATE.tab<-rbind(Emp.CATE,cf.CATE,rf.CATE)
  
  return(CATE.tab)
}


simu.CATE<-function(p,num.simu,case){
  CATE.table<-matrix(NA,nrow=3,ncol=5,byrow = FALSE)
  rownames(CATE.table)<-c("Empirical Tau","Causal Forest","Traditional Regressor")
  colnames(CATE.table)<-c("numob=500","numob=1000","numob=5000","numob=10000","numob=20000")
  
  j<-0
  for (num.ob in c(500,1000,5000,10000,20000)) {
    j<-j+1
    CATE.rep<-matrix(NA,nrow=3,ncol=num.simu)
    for (i in 1:num.simu) {
      data.case<-dgp(N=num.ob,p,case)
      CATE.rep[,i]<-CATE(data.case)
    }
    CATE.table[1,j]<-mean(CATE.rep[1,])
    CATE.table[2,j]<-mean(CATE.rep[2,])
    CATE.table[3,j]<-mean(CATE.rep[3,])
  }
  return(CATE.table)
}

#Case 1
simu.CATE(10,100,case="case 1")
#Case 2
simu.CATE(10,100,case="case 2")
#Case 3
simu.CATE(10,100,case="case 3")
#Case 4
simu.CATE(10,100,case="case 4")

cf.selection<-function(df){
  Y<-df$Y
  W<-df$W
  X<-df[,4:13]
  cf.untrained<-causal_forest(X,Y,W)
  important.var<-variable_importance(cf.untrained)
  select.index<-which(important.var>=median(important.var))
  X.select<-X[,select.index]
  cf.trained<-causal_forest(X.select,Y,W)
  tau.hat<-predict(cf.trained)$predictions
  df$tau.hat<-tau.hat
  
  return(df)
}


simu.MSE<-function(p,num.simu){
  
  MSE.tau.hat<-matrix(NA,nrow=4,ncol=5)
  rownames(MSE.tau.hat)<-c("case 1","case 2","case 3","case 4")
  colnames(MSE.tau.hat)<-c("numob=500","numob=1000","numob=5000","numob=10000","numob=20000")
  
  #case 1
  j<-0
  for (num.ob in c(500,1000,5000,10000,20000)) {
    j<-j+1
    MSE<-c()
    for (i in 1:num.simu) {
      data.case<-dgp(N=num.ob,p,case="case 1")
      data.case<-cf.selection(data.case)
      MSE[i]<-mean((data.case$tau.true-data.case$tau.hat)^2)
    }
    MSE.tau.hat[1,j]<-mean(MSE)
  }
  
  #case 2
  j<-0
  for (num.ob in c(500,1000,5000,10000,20000)) {
    j<-j+1
    MSE<-c()
    for (i in 1:num.simu) {
      data.case<-dgp(N=num.ob,p,case="case 2")
      data.case<-cf.selection(data.case)
      MSE[i]<-mean((data.case$tau.true-data.case$tau.hat)^2)
    }
    MSE.tau.hat[2,j]<-mean(MSE)
  }
  
  #case 3
  j<-0
  for (num.ob in c(500,1000,5000,10000,20000)) {
    j<-j+1
    MSE<-c()
    for (i in 1:num.simu) {
      data.case<-dgp(N=num.ob,p,case="case 3")
      data.case<-cf.selection(data.case)
      MSE[i]<-mean((data.case$tau.true-data.case$tau.hat)^2)
    }
    MSE.tau.hat[3,j]<-mean(MSE)
  }
  
  #case 4
  j<-0
  for (num.ob in c(500,1000,5000,10000,20000)) {
    j<-j+1
    MSE<-c()
    for (i in 1:num.simu) {
      data.case<-dgp(N=num.ob,p,case="case 4")
      data.case<-cf.selection(data.case)
      MSE[i]<-mean((data.case$tau.true-data.case$tau.hat)^2)
    }
    MSE.tau.hat[4,j]<-mean(MSE)
  }
  
  return(MSE.tau.hat)
}

simu.MSE(10,100)
