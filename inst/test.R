
#install.packages('/Users/ruilinli/myglmnet', repo=NULL,type='source')
library(pgenlibr)
library(myglmnet)
n = 100
p = 100

m = PlinkMatrix("/Users/ruilinli/plink-ng/toy_data.pgen", 1:n, 1:p)
m = actualize(m)
#plinktest(m$pgen, rnorm(100))
pgen <- pgenlibr::NewPgen("/Users/ruilinli/plink-ng/toy_data.pgen", pvar = NULL, sample_subset =1:n)
m2 = ReadList(pgen, 1:p, meanimpute =TRUE)
v = rep(1.0, p)
eta = double(n)
mytest(m@ptr, m@xim, n, p, v,eta)
beta = rnorm(p) * rbinom(p,1, 0.3)
y = m2 %*% beta
y = y/sd(y)
a1 = myglmnet(m, y, family='gaussian', standardize=F, intercept=F)
a1 = myglmnet(m, y, family='gaussian', standardize=F, intercept=F)
a1 = myglmnet(m, y, family='gaussian', standardize=F, intercept=F)

a2 = myglmnet(m2, y, family='gaussian', standardize=F, intercept=F)

y2 = rep(1.0, p)
y2[y<median(y)] = 0.0

a1 = myglmnet(m, y2, family='logistic', standardize=F, intercept=F)
a1 = myglmnet(m, y2, family='logistic', standardize=F, intercept=F)
a1 = myglmnet(m, y2, family='logistic', standardize=F, intercept=F)

a2 = myglmnet(m2, y2, family='logistic', standardize=F, intercept=F)



















stop('sufficient')
library(myglmnet)
library(glmnet)
library(microbenchmark)
n = 20
p = 5
 set.seed(1)

X = matrix(rnorm(n*p),n,p)
beta = rnorm(p) * rbinom(p, 1, 0.6)

y = X %*% beta
sdy = sd(y)
y = y /sdy

y2 = rep(1.0, n)
y2[y<median(y)] = 0.0
lam = c(0.3, 0.2, 0.15, 0.1, 0.05)

ref = glmnet(X, y2, family='binomial', lambda = lam, standardize = F, intercept = T)

test = myglmnet(X, y2, family='logistic', lambda = lam,standardize = F, intercept =T)

max(ref$beta[,10] - test$beta[,10])



ref = glmnet(X, y, family='gaussian', exclude=c(2L), standardize = F, intercept = T)

test = myglmnet(X, y, family='gaussian', exclude=c(2L), standardize = F, intercept = T)

## Performance profiling
library(microbenchmark)
library(myglmnet)
library(glmnet)
n = 1500
p = 3000
set.seed(1)
X = matrix(rnorm(n*p),n,p)
beta = rnorm(p) * rbinom(p, 1, 0.3)
y = X%*% beta
y = y/sd(y)
w = rep(1/n, n)
y2 = rep(1.0, n)
y2[y < median(y)] = 0.0
l1 = glmnet(X, y, family='gaussian', exclude=c(2L), weights = w)
l2 = myglmnet(X, y, family='gaussian', exclude=c(2L))

a1 = myglmnet(X, y2, family='logistic', exclude=c(2L))
a2 = glmnet(X, y2, family='binomial', exclude=c(2L))

microbenchmark(
  myglmnet(X, y, family='gaussian', exclude=c(2L)),
  glmnet(X, y, family='gaussian', exclude=c(2L)),
  times=3L
)

microbenchmark(
  myglmnet(X, y2, family='logistic', exclude=c(2L)),
  glmnet(X, y2, family='binomial', exclude=c(2L)),
  times=3L
)






#
# ref = glmnet(X, y, family='gaussian', intercept = F, standardize = F, weights = w)
# lambdar = ref$lambda
# # ref = glmnet(X, y, family='gaussian', lambda=lambdar, intercept = F, standardize = F)
# # tset = wrapper(X, y,lambda=lambdar)
#
# microbenchmark(glmnet(X, y, family='gaussian', lambda=lambdar, intercept = F, standardize = F,weights = w),
#                wrapper(X, y,lambda=lambdar), times=1L)



# X = matrix(rnorm(n*p),n,p)
# beta = rnorm(p) * rbinom(p, 1, 0.6)

# y = X %*% beta
# sdy = sd(y)
# y = y /sdy


# ref = glmnet(X, y, family='gaussian', intercept = T, standardize = F)
# lambdar = ref$lambda
# wrapper(X, y,lambda=lambdar)



# microbenchmark(
#   glmnet(X, y, family='gaussian',lambda = ref$lambda, intercept = T, standardize = F),
#   wrapper(X, y,lambda=lambdar),
#   times = 1L
# )

# ref2 = glmnet(X, y, family='gaussian',lambda = ref$lambda, intercept = T, standardize = F)
#
# #
# # profvis({
# #   ref2 = glmnet(X, y, family='gaussian',lambda = ref$lambda, intercept = F, standardize = F)
# # })
#
# wrapper(X, y,lambda=lambdar)
#
#
#
# ref =glmnet(X,y, family='gaussian', lambda=lambdar,intercept = T, standardize = F)
# ref$beta
#
# # ref = glmnet()
#
# # wrapper()
