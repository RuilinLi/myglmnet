#include "R.h"
#include "R_ext/Print.h"
#include "Rinternals.h"
#include "glmfamily.h"
#include "glmnetMatrix.h"

void glmnetPath(double alpha, MatrixGlmnet *X, const double *y, const double *v,
                int intr, const int *ju, const double *vp, const double *cl,
                int nx, double thr, int maxit, GlmFamily *fam, bool has_offset,
                double *offset, const double *lambdas, int nlambda, int mxitnr,
                const double flmin, int *lmu, double *a0, double *ca, int *ia,
                int *nin, double *devratio_vec, double *alm, int *nlp,
                double *nulldev, int *jerr);

GlmFamily *get_family(const char *family) {
    if (strcmp(family, "gaussian") == 0) {
        return new Gaussian();
    }
    if (strcmp(family, "logistic") == 0) {
        return new Logistic();
    }
    error("invalid family name\n");
}

void get_xmxs_dense(const double *x, const double *v, int *ju, double *xm,
                    double *xs, int no, int ni) {
    for (int j = 0; j < ni; ++j) {
        if (!ju[j]) {
            continue;
        }

        double ex = 0.0;
        double ex2 = 0.0;
        for (int i = 0; i < no; ++i) {
            ex += v[i] * x[j * no + i];
            ex2 += v[i] * x[j * no + i] * x[j * no + i];
        }
        // Assumes sum_{i} v_i = 1
        double variance = ex2 - ex * ex;
        if (variance <= 0) {
            ju[j] = 0;
            continue;
        }
        xm[j] = ex;
        xs[j] = sqrt(variance);
    }
}

void standardize(double *x, const int *ju, double *xm, double *xs, int intr,
                 int isd, int no, int ni) {
    if (intr) {
        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) {
                continue;
            }
            for (int i = 0; i < no; ++i) {
                x[j * no + i] -= xm[j];
            }
        }

        if (isd) {
            for (int j = 0; j < ni; ++j) {
                if (!ju[j]) {
                    continue;
                }
                for (int i = 0; i < no; ++i) {
                    x[j * no + i] /= xs[j];
                }
            }
        }
        return;
    }

    if (isd) {
        for (int j = 0; j < ni; ++j) {
            if (!ju[j]) {
                continue;
            }
            for (int i = 0; i < no; ++i) {
                x[j * no + i] /= xs[j];
            }
        }
    }
}

#ifdef __cplusplus
extern "C" {
#endif

SEXP solve(SEXP alpha2, SEXP x2, SEXP y2, SEXP weights2, SEXP ju2, SEXP vp2,
           SEXP cl2, SEXP nx2, SEXP nlam2, SEXP flmin2, SEXP ulam2,
           SEXP thresh2, SEXP isd2, SEXP intr2, SEXP maxit2, SEXP lmu2,
           SEXP a02, SEXP ca2, SEXP ia2, SEXP nin2, SEXP devratio2, SEXP alm2,
           SEXP nlp2, SEXP family2, SEXP offset2, SEXP has_offset2,
           SEXP mxitnr2, SEXP nulldev2, SEXP jerr2) {
    // Create matrix object
    int no = nrows(x2);
    int ni = ncols(x2);
    double *xptr = REAL(x2);
    int intr = asInteger(intr2);
    int isd = asInteger(isd2);
    int *ju = INTEGER(ju2);
    double *v = REAL(weights2);
    bool dup_x = (bool)(intr + isd);

    // always compute xm and xs, maybe won't be used
    double *xm = (double *)malloc(sizeof(double) * ni);
    double *xs = (double *)malloc(sizeof(double) * ni);

    get_xmxs_dense(REAL(x2), v, ju, xm, xs, no, ni);

    if (dup_x) {
        xptr = REAL(PROTECT(duplicate(x2)));
        standardize(xptr, ju, xm, xs, intr, isd, no, ni);
    }

    DenseM X(no, ni, xptr);

    // Create family object
    GlmFamily *fam = get_family(CHAR(STRING_ELT(family2, 0)));
    double flmin = asReal(flmin2);
    int nlambda = asInteger(nlam2);
    double *lambdas = REAL(ulam2);

    double alpha = asReal(alpha2);
    double *y = REAL(y2);

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

    // Things that needs to be modified
    int *lmu = INTEGER(lmu2);
    double *a0 = REAL(a02);
    double *ca = REAL(ca2);
    int *ia = INTEGER(ia2);
    int *nin = INTEGER(nin2);
    double *devratio = REAL(devratio2);
    double *alm = REAL(alm2);
    int *nlp = INTEGER(nlp2);
    int *jerr = INTEGER(jerr2);
    double *nulldev = REAL(nulldev2);

    glmnetPath(alpha, &X, y, v, intr, ju, vp, cl, nx, thr, maxit, fam,
               has_offset, offset, lambdas, nlambda, mxitnr, flmin, lmu, a0, ca,
               ia, nin, devratio, alm, nlp, nulldev, jerr);

    // scale parameters back if standardization happend
    if (isd) {
        for (int m = 0; m < (*lmu); ++m) {
            for (int k = 0; k < nin[m]; ++k) {
                ca[m * nx + k] /= xs[ia[k]];
            }
        }
    }

    if (intr) {
        for (int m = 0; m < (*lmu); ++m) {
            for (int k = 0; k < nin[m]; ++k) {
                a0[m] -= ca[m * nx + k] * xm[ia[k]];
            }
        }
    }

    delete fam;
    free(xm);
    free(xs);
    if (dup_x) {
        UNPROTECT(1);
    }

    return R_NilValue;
}

#ifdef __cplusplus
}
#endif