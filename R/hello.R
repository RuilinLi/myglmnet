
#' @useDynLib myglmnet, .registration = TRUE
#' @export
wrapper <- function(x, y, lambda) {
    no <- nrow(x)
    ni <- ncol(x)
    nlam <- length(lambda)
    weight <- rep(1/no, no)
    vp <- rep(1, ni)
    lower.limits <- rep(-100, ni)
    upper.limits <- rep(100, ni)
    cl <- rbind(lower.limits, upper.limits)
    ju <- rep(1L, ni)
    # .Call('test', x, y, lambda, weight, 0L, ju, vp, cl, PACKAGE = 'myglmnet')
    NULL
}
#' @export
mytest = function(xptr, xim, no, ni, v, eta) {
    .Call('testplink',xptr, no, ni, xim, v, eta)
}
#tools::package_native_routine_registration_skeleton('/Users/ruilinli/myglmnet', con='/Users/ruilinli/myglmnet/src/init.c')
#devtools::document('/Users/ruilinli/myglmnet')
#install.packages('/Users/ruilinli/myglmnet', repo=NULL,type='source')

#' @export
myglmnet <- function(x, y, family = c("gaussian", "logistic"), weights = NULL, offset = NULL, 
    alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04), 
    lambda = NULL, standardize = TRUE, intercept = TRUE, thresh = 1e-07, dfmax = nvars + 
        1, pmax = min(dfmax * 2 + 20, nvars), exclude = NULL, penalty.factor = rep(1, 
        nvars), lower.limits = -Inf, upper.limits = Inf, maxit = 1e+05) {
    
    this.call <- match.call()
    ### Need to do this first so defaults in call can be satisfied
    np <- dim(x)
    ## check dims
    if (is.null(np) | (np[2] <= 1)) 
        stop("x should be a matrix with 2 or more columns")
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    family <- match.arg(family)
    if (alpha > 1) {
        warning("alpha >1; set to 1")
        alpha <- 1
    }
    if (alpha < 0) {
        warning("alpha<0; set to 0")
        alpha <- 0
    }
    alpha <- as.double(alpha)
    nlam <- as.integer(nlambda)
    y <- drop(y)  # we dont like matrix responses unless we need them
    if (is.null(weights)) 
        weights <- rep(1/nobs, nobs) else if (length(weights) != nobs) 
        stop(paste("number of elements in weights (", length(weights), ") not equal to the number of rows of x (", 
            nobs, ")", sep = "")) else weights <- weights/(sum(weights))
    dimy <- dim(y)
    nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
    if (nrowy != nobs) 
        stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (", 
            nobs, ")", sep = ""))
    vnames <- colnames(x)
    if (is.null(vnames)) 
        vnames <- paste("V", seq(nvars), sep = "")
    ne <- as.integer(dfmax)
    nx <- as.integer(pmax)
    if (is.null(exclude)) 
        exclude <- integer(0)
    if (any(penalty.factor == Inf)) {
        exclude <- c(exclude, seq(nvars)[penalty.factor == Inf])
        exclude <- sort(unique(exclude))
    }
    ju <- rep(1L, nvars)
    if (length(exclude) > 0) {
        jd <- match(exclude, seq(nvars), 0)
        if (!all(jd > 0)) 
            stop("Some excluded variables out of range")
        penalty.factor[jd] <- 1  #ow can change lambda sequence
        ju[jd] <- 0L
    }
    vp <- as.double(penalty.factor)
    internal.parms <- glmnet.control()
    
    # if (internal.parms$itrace) trace.it <- 1 else { if (trace.it) {
    # glmnet.control(itrace = 1) on.exit(glmnet.control(itrace = 0)) } } check on
    # limits
    if (any(lower.limits > 0)) {
        stop("Lower limits should be non-positive")
    }
    if (any(upper.limits < 0)) {
        stop("Upper limits should be non-negative")
    }
    lower.limits[lower.limits == -Inf] <- -internal.parms$big
    upper.limits[upper.limits == Inf] <- internal.parms$big
    if (length(lower.limits) < nvars) {
        if (length(lower.limits) == 1) 
            lower.limits <- rep(lower.limits, nvars) else stop("Require length 1 or nvars lower.limits")
    } else lower.limits <- lower.limits[seq(nvars)]
    if (length(upper.limits) < nvars) {
        if (length(upper.limits) == 1) 
            upper.limits <- rep(upper.limits, nvars) else stop("Require length 1 or nvars upper.limits")
    } else upper.limits <- upper.limits[seq(nvars)]
    cl <- rbind(lower.limits, upper.limits)
    if (any(cl == 0)) {
        ### Bounds of zero can mess with the lambda sequence and fdev; ie nothing happens
        ### and if fdev is not zero, the path can stop fdev <- glmnet.control()$fdev if
        ### (fdev != 0) { glmnet.control(fdev = 0) on.exit(glmnet.control(fdev = fdev)) }
        stop("This case has not been implemented")
    }
    storage.mode(cl) <- "double"
    ### end check on limits
    
    isd <- as.integer(standardize)
    intr <- as.integer(intercept)
    if (!missing(intercept) && family == "cox") 
        warning("Cox model has no intercept")
    
    # Don't have this yet jsd <- as.integer(standardize.response)
    thresh <- as.double(thresh)
    if (is.null(lambda)) {
        if (lambda.min.ratio >= 1) 
            stop("lambda.min.ratio should be less than 1")
        flmin <- as.double(lambda.min.ratio)
        ulam <- double(1)
    } else {
        flmin <- as.double(1)
        if (any(lambda < 0)) 
            stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    is.sparse <- FALSE
    ix <- jx <- NULL
    if (inherits(x, "sparseMatrix")) {
        stop("sparse matrices not implemented yet")
        ## Sparse case
        is.sparse <- TRUE
        x <- as(x, "CsparseMatrix")
        x <- as(x, "dgCMatrix")
        ix <- as.integer(x@p + 1)
        jx <- as.integer(x@i + 1)
        x <- as.double(x@x)
    }
    
    maxit <- as.integer(maxit)
    lmu <- integer(1)
    a0 <- double(nlam)
    ca <- double(nx * nlam)
    ia <- integer(nx)
    nin <- integer(nlam)
    devratio <- double(nlam)
    alm <- double(nlam)
    nlp <- integer(1)
    jerr <- integer(1)
    nulldev <- double(1)
    
    if (is.null(offset)) {
        offset <- double(0)
        has_offset <- 0L
    } else if (length(offset) != nobs) {
        stop("The length of offset must be the same as the number of observations")
    } else {
        has_offset <- 1L
    }
    mxitnr <- internal.parms$mxitnr
    if (family == "gaussian") {
        mxitnr <- 1
    }
    mxitnr <- as.integer(mxitnr)

    if(inherits(x, "PlinkMatrix")){
        x_list = list("Plink", np[1], np[2], x@ptr, x@xim)
    } else {
        x_list = list("Dense", np[1], np[2], x)
    }
    
    .Call("solve", alpha, x_list, y, weights, ju, vp, cl, nx, nlam, flmin, ulam, thresh, 
        isd, intr, maxit, lmu, a0, ca, ia, nin, devratio, alm, nlp, family, offset, 
        has_offset, mxitnr, nulldev, jerr)
    
    # if (trace.it) { if (relax) cat('Training Fit\n') pb <- createPB(min = 0, max =
    # nlam, initial = 0, style = 3) }
    fit <- list(jerr = jerr, a0 = a0, nin = nin, lmu = lmu, alm = alm, ca = ca, ia = ia)
    outlist <- getcoef(fit, nvars, nx, vnames)
    dev <- devratio[seq(lmu)]
    outlist <- c(outlist, list(dev.ratio = dev, nulldev = nulldev, npasses = nlp, 
        jerr = jerr, offset = has_offset))
    
    outlist$call <- this.call
    outlist$nobs <- nobs
    class(outlist) <- c(class(outlist), "glmnet")
    outlist
    
    # No need for this either kopt <- switch(match.arg(type.logistic), Newton = 0,
    # modified.Newton = 1) if (family == 'multinomial') { type.multinomial <-
    # match.arg(type.multinomial) if (type.multinomial == 'grouped') kopt <- 2 } kopt
    # <- as.integer(kopt)
    
}
