#include "R.h"
#include "R_ext/Print.h"
#include "Rinternals.h"
#include "glmnetMatrix.h"

void glmnetPath(double alpha, MatrixGlmnet *X, const double *y, const double *v,
                int intr, const int *ju, const double *vp, const double *cl,
                int nx, double thr, int maxit, const char *family,
                bool has_offset, double *offset, const double *lambdas,
                int nlambda);


#ifdef __cplusplus
extern "C" {
#endif

SEXP test(SEXP x2, SEXP y2, SEXP lambda2, SEXP v2, SEXP intr2, SEXP ju2,
          SEXP vp2, SEXP cl2) {
    int no = nrows(x2);
    int ni = ncols(x2);
    DenseM X(no, ni, REAL(x2));
    double *y = REAL(y2);
    double *lambda = REAL(lambda2);
    double *v = REAL(v2);
    int intr = asInteger(intr2);
    int *ju = INTEGER(ju2);
    double *vp = REAL(vp2);
    double *cl = REAL(cl2);
    int nlambda = length(lambda2);
    double alpha = 1.0;
    // Rprintf("cl lower is %f \n", cl[]);
    // Rprintf("cl upper is %f\n")
    glmnetPath(alpha, &X, y, v, intr, ju, vp, cl, ni, 1e-7, 10000, "logistic", false,
               nullptr, lambda, nlambda);
    return R_NilValue;
}
#ifdef __cplusplus
}
#endif