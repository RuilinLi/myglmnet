#include "R.h"
#include "R_ext/Print.h"
#include "Rinternals.h"
#include "glmfamily.h"
#include "glmnetMatrix.h"
// #include "gperftools/profiler.h"

void glmnetPath(double alpha, MatrixGlmnet *X, const double *y, const double *v,
                int intr, const int *ju, const double *vp, const double *cl,
                int nx, double thr, int maxit, GlmFamily *fam, bool has_offset,
                double *offset, const double *lambdas, int nlambda, int mxitnr,
                const double flmin, int *lmu, double *a0, double *ca, int *ia,
                int *nin, double *devratio_vec, double *alm, int *nlp,
                double *nulldev, int *jerr, double *a, int *iy, int *mm, int nino, int warm, double *residuals);

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

MatrixGlmnet *get_matrix(SEXP xptr, const char *mattype, int no, int ni,
                         int isd, int intr, int *ju, double *xm, double *xs,
                         const double *v, const double *xim, int ncov, const double *cov) {
    if (strcmp(mattype, "Dense") == 0) {
        double *x = REAL(xptr);
        get_xmxs_dense(x, v, ju, xm, xs, no, ni);
        standardize(x, ju, xm, xs, intr, isd, no, ni);
        MatrixGlmnet *result = new DenseM(no, ni, x);
        return result;
    }

    if (strcmp(mattype, "Plink") == 0) {
        uintptr_t *x = (uintptr_t *)R_ExternalPtrAddr(xptr);
        MatrixGlmnet *result = new PlinkMatrix(no, ni, x, xim, intr, ncov, cov);
        return result;
    }

    error("invalid matrix type\n");
}

#ifdef __cplusplus
extern "C" {
#endif

SEXP testplink(SEXP x2, SEXP no2, SEXP ni2, SEXP xim2, SEXP v2, SEXP eta2) {
    uintptr_t *xptr = (uintptr_t *)R_ExternalPtrAddr(x2);
    MatrixGlmnet *x =
        new PlinkMatrix(asInteger(no2), asInteger(ni2), xptr, REAL(xim2), 0, 0, nullptr);
    double *eta = REAL(eta2);
    double *v = REAL(v2);
    x->compute_eta(eta, v, 0.0,
                              false, nullptr);

    return R_NilValue;
}

SEXP solve(SEXP alpha2, SEXP x2, SEXP y2, SEXP weights2, SEXP ju2, SEXP vp2,
           SEXP cl2, SEXP nx2, SEXP nlam2, SEXP flmin2, SEXP ulam2,
           SEXP thresh2, SEXP isd2, SEXP intr2, SEXP maxit2, SEXP lmu2,
           SEXP a02, SEXP ca2, SEXP ia2, SEXP nin2, SEXP devratio2, SEXP alm2,
           SEXP nlp2, SEXP family2, SEXP offset2, SEXP has_offset2,
           SEXP mxitnr2, SEXP nulldev2, SEXP jerr2, SEXP beta02, SEXP iy2,
           SEXP mm2, SEXP nino2, SEXP warm2, SEXP residuals2) {
    // Create matrix object
    // ProfilerStart("/home/ruilinli/myglmnet/inst/prof.out");
    const char *mattype = CHAR(STRING_ELT(VECTOR_ELT(x2, 0), 0));
    int no = asInteger(VECTOR_ELT(x2, 1));
    int ni = asInteger(VECTOR_ELT(x2, 2));
    SEXP xptr = VECTOR_ELT(x2, 3);

    double *xim = nullptr;  // only for plink matrix
    double *cov = nullptr; // only for plink matrix
    int ncov = 0;
    if (strcmp(mattype, "Plink") == 0) {
        xim = REAL(VECTOR_ELT(x2, 4));
        ncov = asInteger(VECTOR_ELT(x2, 5));
        if(ncov > 0) {
            cov = REAL(VECTOR_ELT(x2, 6));
        }
    }


    int intr = asInteger(intr2);
    int isd = asInteger(isd2);
    int *ju = INTEGER(ju2);
    double *v = REAL(weights2);
    bool dup_x = (bool)(intr + isd);
    dup_x = dup_x && (strcmp(mattype, "Dense") == 0);
    if (dup_x) {
        xptr = PROTECT(duplicate(xptr));
    }

    // always compute xm and xs, maybe won't be used
    double *xm = (double *)malloc(sizeof(double) * ni);
    double *xs = (double *)malloc(sizeof(double) * ni);

    MatrixGlmnet *X =
        get_matrix(xptr, mattype, no, ni, isd, intr, ju, xm, xs, v, xim, ncov, cov);

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
    double *residuals = REAL(residuals2);

    double *beta0 = REAL(PROTECT(duplicate(beta02)));
    int *iy = INTEGER(iy2);
    int *mm = INTEGER(mm2);
    int nino = asInteger(nino2);
    int warm = asInteger(warm2);

    if(warm && isd) {
        for(int j = 0; j < ni; ++j) {
            beta0[j] *= xs[j];
        }
    }

    glmnetPath(alpha, X, y, v, intr, ju, vp, cl, nx, thr, maxit, fam,
               has_offset, offset, lambdas, nlambda, mxitnr, flmin, lmu, a0, ca,
               ia, nin, devratio, alm, nlp, nulldev, jerr, beta0, iy, mm, nino, warm, residuals);

    // scale parameters back if standardization happend
    if (isd && (strcmp(mattype, "Plink") != 0)) {
        for (int m = 0; m < (*lmu); ++m) {
            for (int k = 0; k < nin[m]; ++k) {
                ca[m * nx + k] /= xs[ia[k]];
            }
        }
    }

    if (intr)
    {
        if (strcmp(mattype, "Plink") == 0)
        {
            for (int m = 0; m < (*lmu); ++m)
            {
                for (int k = 0; k < nin[m]; ++k)
                {
                    a0[m] -= ca[m * nx + k] * xim[ia[k]];
                }
            }
        }
        else
        {
            for (int m = 0; m < (*lmu); ++m)
            {
                for (int k = 0; k < nin[m]; ++k)
                {
                    a0[m] -= ca[m * nx + k] * xm[ia[k]];
                }
            }
        }
    }

    delete fam;
    free(xm);
    free(xs);
    UNPROTECT(1);
    if (dup_x) {
        UNPROTECT(1);
    }
    delete X;
    //ProfilerStop();
    return R_NilValue;
}

#ifdef __cplusplus
}
#endif