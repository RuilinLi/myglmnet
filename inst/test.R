library(myglmnet)
library(glmnet)
library(microbenchmark)
n = 30
p = 5
set.seed(1)

X = matrix(rnorm(n*p),n,p)
beta = rnorm(p) * rbinom(p, 1, 0.3)

y = X %*% beta
sdy = sd(y)
y = y /sdy
w = rep(1/n, n)
y2 = rep(1.0, n)
y2[y < median(y)] = 0.0

ref = glmnet(X, y2, family='binomial', intercept = F, standardize = F)
lambdar = ref$lambda

wrapper(X, y2, lambdar)

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