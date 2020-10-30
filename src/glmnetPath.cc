#include <math.h>
#include <string.h>

#include "R.h"
#include "glmfamily.h"
#include "glmnetMatrix.h"

void wls_base(double alm0, double almc, double alpha, int m, int no, int ni,
              MatrixGlmnet *X, double *r, const double *v, int intr,
              const int *ju, const double *vp, const double *cl, int nx,
              double thr, int maxit, double *__restrict a, double *aint,
              double *__restrict g, int *__restrict ia, int *__restrict iy,
              int *iz, int *__restrict mm, int *nino, double *rsqc, int *nlp,
              double *__restrict xv, int *jerr);

GlmFamily *get_family(const char *family) {
    if (strcmp(family, "gaussian") == 0) {
        Rprintf("this one!\n");
        return new Gaussian();
    }
    if (strcmp(family, "logistic") == 0) {
        return new Logistic();
    }

    return nullptr;
}

void glmnetPath(double alpha, MatrixGlmnet *X, const double *y, const double *v,
                int intr, const int *ju, const double *vp, const double *cl,
                int nx, double thr, int maxit, const char *family,
                bool has_offset, double *offset, const double *lambdas,
                int nlambda) {
    GlmFamily *fam = get_family(family);
    int no = X->get_no();
    int ni = X->get_ni();
    double aint = 0;  // intercept

    double *r = (double *)malloc(sizeof(double) * no);
    double *w = (double *)malloc(sizeof(double) * no);
    double *z = (double *)malloc(sizeof(double) * no);
    double *eta = (double *)malloc(sizeof(double) * no);

    double *a = (double *)malloc(sizeof(double) * ni);
    double *g = (double *)malloc(sizeof(double) * ni);
    double *xv = (double *)malloc(sizeof(double) * ni);
    int *ia = (int *)malloc(sizeof(int) * ni);
    int *iy = (int *)malloc(sizeof(int) * ni);
    int *mm = (int *)malloc(sizeof(int) * ni);

    for (int i = 0; i < ni; ++i) {
        iy[i] = 0;
        mm[i] = 0;
        a[i] = 0.0;
    }

    int iz = 0;
    int nino = 0;
    int nlp = 0;
    int jerr = 0;
    double rsqc = 0;

    double nulldev =
        fam->null_deviance(y, v, intr, eta, has_offset, offset, &aint, no);

    int mxitnr = (strcmp(family, "gaussian") == 0) ? 1 : 25;
    Rprintf("maxiter is %d \n", mxitnr);
    double alm0 = 0;
    double prev_dev = nulldev;
    double current_dev;
    for (int m = 0; m < nlambda; ++m) {
        double almc = lambdas[m];
        // Rprintf("Currently lambda is %f\n", almc);
        // Rprintf("Currently nlp is %d\n", nlp);
        // Rprintf("Currently dev is %f\n", prev_dev);
        // Rprintf("Currently a3 is %f\n", a[2]);
        // Rprintf("Currently nino is %d\n", nino);
        // for (int l = 0; l < no; ++l) {
        //     Rprintf("eta[%d] is %f\n", l + 1, eta[l]);
        // }

        // for (int l = 0; l < ni; ++l) {
        //     Rprintf("beta[%d] is %f\n", l + 1, a[l]);
        // }

        // IRLS, z here is actually the weighted working response
        for (int j = 0; j < mxitnr; ++j) {
            fam->get_workingset(eta, y, v, w, z, no);
            wls_base(alm0, almc, alpha, m, no, ni, X, z, w, intr, ju, vp, cl,
                     nx, thr, maxit, a, &aint, g, ia, iy, &iz, mm, &nino, &rsqc,
                     &nlp, xv, &jerr);

            X->compute_eta(eta, a, aint, has_offset, offset);
            current_dev = fam->get_deviance(y, eta, v, no);

            double rel_diff =
                (prev_dev - current_dev) / fmax(fabs(prev_dev), 1.0);
            prev_dev = current_dev;
            if (fabs(rel_diff) < thr) {
                break;
            }
            // TODO update eta = X * a + aint + offset
            // Use relative change in deviance to determine IRLS convergence
            // Use devratio to determine early-stop in lambda iteration
        }
        double devratio = 1 - current_dev / nulldev;
        if (devratio > 0.999) {
            break;
        }
        alm0 = almc;
        // for (int l = 0; l < ni; ++l) {
        //     Rprintf("beta[%d] is %f\n", l + 1, a[l]);
        // }
        // Rprintf("Currently nlp is %d\n", nlp);

    }
    delete fam;  // Mixed C and C++ style heap allocation...
    free(r);
    free(a);
    free(g);
    free(w);
    free(z);
    free(eta);
    free(xv);
    free(ia);
    free(iy);
    free(mm);

    return;
}