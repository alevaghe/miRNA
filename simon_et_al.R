rm(list = ls())
cat("\014")

library(plyr)
library(Coxnet)
library(survival)
library(survminer)
library(glmnet)
library(MASS)

#function for the sets dimensions
k.fold = function(dati, K) {
  #number of individuals
  n = nrow(dati)
  
  #remaining
  rem = n %% K
  
  #sets dimensions
  dimm = rep(floor(n / K), K)
  
  #allocating the remaining individuals
  if (rem != 0) {
    #random allocation of the remaining individuals
    a = sample(K, rem)
    dimm[a] = dimm[a] + 1
  }
  return(dimm)
}

#reading miRNA data
miRNA = as.data.frame(read.table("I:/IRST Gruppi/Biostatistica/Personali/Alessandro/miRNA/allmirna.txt"))
covariates = as.character(miRNA[, 1])
miRNA = as.data.frame(t(miRNA[, 2:ncol(miRNA)]))
names(miRNA) = covariates

#reading outcome data
df = as.data.frame(read.table("I:/IRST Gruppi/Biostatistica/Personali/Alessandro/miRNA/outcome_stadio_dummy.txt", header = TRUE))

#merging the data
df = cbind(df[, c(1, 4, 3, 6, 7)], miRNA[, - 1])
rownames(df) = NULL
colnames(df)[3] = 'status'

#sets dimensions
dimm = k.fold(df, nrow(df))

#splitting the data.frame by rows
v = NULL
for (ii in 1:(length(dimm) - 1)) {
  v = c(v, sum(dimm[1:ii]))
}
sets = split(df, cumsum(1:nrow(df) %in% (v + 1)))
names(sets) = 1:length(dimm)

#initializing objects for the cross-validation
cox.cv = cox.coeff = pred.indices = medians = vector('list')
hl = rep(0, nrow(df))

#saving start time
ptm = proc.time()

#cross-validation (simon et al. 2011)
for (ii in 1:length(dimm)) {
  #training and test split
  test = sets[[ii]]
  training = ldply(sets[c(1:length(dimm))[- ii]], data.frame)[, - 1]
  
  #matrices for the Cox model
  x = as.matrix(training[, 4:ncol(training)])
  y = as.matrix(training[, c(2, 3)])
  
  #Cox regularised model
  #set.seed(1234)
  cv.res = cv.glmnet(x, y, family = "cox", nfolds = 10, maxit = 10^7)
  cv.coeff = as.vector(coef(cv.res, s = "lambda.min"))
  
  #predictive index
  predictive = cv.coeff %*% t(x)
  
  #median predictive
  med = median(predictive)
  
  #identifying high and low risks
  hl[test[which(cv.coeff %*% t(as.matrix(test[, 4:ncol(training)])) > med), 1]] = 1
  
  #saving results
  cox.cv[[ii]] = cv.res
  cox.coeff[[ii]] = cv.coeff
  pred.indices[[ii]] = predictive
  medians[[ii]] = med
  
  #printing loop counter
  print(paste('individual = ', ii, sep = ''))
  
  #printing time 
  print(proc.time() - ptm)
  
  print('')
}

#initializing counter
count = rep(0, length(dimm))

#counting
for (jj in 1:length(dimm)) {
  count = count + (1:83 %in% which(as.vector(cox.coeff[[jj]]) != 0))
}

#not ever zero covariates
perc = count[count > 0] / 83
names(perc) = colnames(df)[which(count / 83 > 0) + 3]

