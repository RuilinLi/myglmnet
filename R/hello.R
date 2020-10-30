#' @useDynLib myglmnet, .registration = TRUE
#' @export
wrapper = function(x, y, lambda){
  no = nrow(x)
  ni = ncol(x)
  nlam = length(lambda)
  weight = rep(1/no, no)
  vp = rep(1.0, ni)
  lower.limits = rep(-100.0, ni)
  upper.limits = rep(100.0, ni)
  cl = rbind(lower.limits, upper.limits)
  ju = rep(1L, ni)
  .Call("test", x, y, lambda, weight, 1L, ju, vp, cl, PACKAGE = "myglmnet")
}
