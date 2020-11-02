formatR::tidy_file('/Users/ruilinli/myglmnet/R/hello.R', arrow=TRUE)
tools::package_native_routine_registration_skeleton('/Users/ruilinli/myglmnet', con='/Users/ruilinli/myglmnet/src/init.c')
devtools::document('/Users/ruilinli/myglmnet')
install.packages('/Users/ruilinli/myglmnet', repo=NULL, type='source')


library(myglmnet)
n = 6
p = 5
set.seed(1)

X = matrix(rnorm(n*p),n,p)
beta = rnorm(p) * rbinom(p, 1, 0.99)

y = X %*% beta
sdy = sd(y)
y = y /sdy

result = myglmnet(X, y, family='gaussian', exclude=c(2L), intercept=TRUE, standardize=FALSE)
