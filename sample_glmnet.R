library(MASS)
library(glmnet)
set.seed(101)
n=1000
p = 1400
X = matrix(rnorm(n*p,seq(-2,2,length=n*p),1),ncol=p)
beta = rpois(p,3)/20
y = rbinom(dim(X)[1], 1,exp(1 + X%*%beta)/(1+exp(1 + X%*%beta)))
set.seed (101)
alpha   = 0.5 #mixing propotion
# mixing propotion measure can be changed and described as above
cv.out  = cv.glmnet(X,y,family="binomial",alpha=alpha,nfolds=10,keep=T) 
x0=matrix(rnorm(p),nrow=1)
predict(cv.out, newx=x0, s="lambda.1se") # link
predict(cv.out, newx=x0, s="lambda.1se",type="response")       # Fitted probabilities, estimates Pr[y=1]
predict(cv.out, newx=x0, s="lambda.1se",type="class")          # Predict class. If pr >0.5 then y = 1
predict(cv.out, newx=X[y==1,], s="lambda.1se",type="response") 
all.equal(as.character(round(predict(cv.out, newx=X, s="lambda.1se",type="response"))), 
          as.vector(predict(cv.out, newx=X, s="lambda.1se",type="class"))) # Must be iden