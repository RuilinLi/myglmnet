#include "R.h"
#include "R_ext/Print.h"
#include "Rinternals.h"
#include "glmfamily.h"
#include "glmnetMatrix.h"

void glmnetPath(double alpha, MatrixGlmnet *X, const double *y, const double *v,
                int intr, const int *ju, const double *vp, const double *cl,
                int nx, double thr, int maxit, GlmFamily *fam,
                bool has_offset, double *offset, const double *lambdas,
                int nlambda, int mxitnr, const double flmin);

GlmFamily *get_family(const char *family) {
    if (strcmp(family, "gaussian") == 0) {
        Rprintf("this one!\n");
        return new Gaussian();
    }
    if (strcmp(family, "logistic") == 0) {
        return new Logistic();
    }
    error("invalid family name\n");
}

#ifdef __cplusplus
extern "C" {
#endif

// SEXP test(SEXP x2, SEXP y2, SEXP lambda2, SEXP v2, SEXP intr2, SEXP ju2,
//           SEXP vp2, SEXP cl2) {
//     int no = nrows(x2);
//     int ni = ncols(x2);
//     DenseM X(no, ni, REAL(x2));
//     double *y = REAL(y2);
//     double *lambda = REAL(lambda2);
//     double *v = REAL(v2);
//     int intr = asInteger(intr2);
//     int *ju = INTEGER(ju2);
//     double *vp = REAL(vp2);
//     double *cl = REAL(cl2);
//     int nlambda = length(lambda2);
//     double alpha = 1.0;
//     // Rprintf("cl lower is %f \n", cl[]);
//     // Rprintf("cl upper is %f\n")
//     glmnetPath(alpha, &X, y, v, intr, ju, vp, cl, ni, 1e-7, 10000,
//     "logistic",
//                false, nullptr, lambda, nlambda);
//     return R_NilValue;
// }

SEXP solve(SEXP alpha2, SEXP x2, SEXP y2, SEXP weights2, SEXP ju2, SEXP vp2,
           SEXP cl2, SEXP nx2, SEXP nlam2, SEXP flmin2, SEXP ulam2,
           SEXP thresh2, SEXP isd2, SEXP intr2, SEXP maxit2, SEXP lmu2,
           SEXP a02, SEXP ca2, SEXP nin2, SEXP devratio2, SEXP alm2, SEXP nlp2,
           SEXP family2, SEXP offset2, SEXP has_offset2, SEXP mxitnr2) {
    // Create matrix object
    int no = nrows(x2);
    int ni = ncols(x2);
    DenseM X(no, ni, REAL(x2));

    // Create family object
    GlmFamily *fam = get_family(CHAR(STRING_ELT(family2, 0)));
    double flmin = asReal(flmin2);
    int nlambda = asInteger(nlam2);
    double *lambdas;


    double alpha = asReal(alpha2);
    double *y = REAL(y2);
    double *v = REAL(weights2);
    int intr = asInteger(intr2);
    int *ju = INTEGER(ju2);
    double *vp = REAL(vp2);
    double *cl = REAL(cl2);
    int nx = asInteger(nx2);
    double thr = asReal(thresh2);
    int maxit = asInteger(maxit2);
    bool has_offset = (bool)asInteger(has_offset2);
    double *offset = nullptr;
    if (has_offset) {
        offset = REAL(offset2);
    }
    int mxitnr = asInteger(mxitnr2);


    glmnetPath(alpha, &X, y, v, intr, ju, vp, cl, nx, thr, maxit, fam,
               has_offset, offset, lambdas, nlambda, mxitnr, flmin);

    delete fam;

    return R_NilValue;
}

#ifdef __cplusplus
}
#endif