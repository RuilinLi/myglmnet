library(myglmnet)
library(glmnet)
library(microbenchmark)
n = 3000
p = 5000
set.seed(1)

X = matrix(rnorm(n*p),n,p)
beta = rnorm(p) * rbinom(p, 1, 0.6)

y = X %*% beta
sdy = sd(y)
y = y /sdy


ref = glmnet(X, y, family='gaussian', intercept = T, standardize = F)
lambdar = ref$lambda

microbenchmark(
  glmnet(X, y, family='gaussian',lambda = ref$lambda, intercept = T, standardize = F),
  wrapper(X, y,lambda=lambdar),
  times = 1L
)

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
