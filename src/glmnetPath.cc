#include <math.h>
#include <string.h>
#include <assert.h>
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


void glmnetPath(double alpha, MatrixGlmnet *X, const double *y, const double *v,
                int intr, const int *ju, const double *vp, const double *cl,
                int nx, double thr, int maxit, GlmFamily *fam,
                bool has_offset, double *offset, const double *ulambdas,
                int nlambda, int mxitnr, const double flmin) {
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
        fam->null_deviance(y, v, r, intr, eta, has_offset, offset, &aint, no);

    Rprintf("maxiter is %d \n", mxitnr);

    // Compute max_lambda here instead of using user defined lambdas
    double *lambdas = (double*)malloc(sizeof(double) * nlambda);
    if(flmin < 1.0) {
        double max_lambda = X->max_grad(r, ju, vp);
        lambdas = (double*)malloc(sizeof(double) * nlambda);
        lambdas[0] = max_lambda;
        double ratio = pow(flmin, 1.0/(nlambda - 1));
        for(int i = 1; i < nlambda; ++i) {
            lambdas[i] = lambdas[i-1] * ratio;
            Rprintf("lambda[%d] is %f\n", i+1, lambdas[i]);
        }
    } else {
        for(int i = 0; i < nlambda; ++i) {
            lambdas[i] = ulambdas[i];
        }
    }

    double alm0 = 0;
    double prev_dev = nulldev;
    double current_dev;
    for (int m = 0; m < nlambda; ++m) {

        double almc = lambdas[m];

        // IRLS, z here is actually the weighted working response
        for (int j = 0; j < mxitnr; ++j) {
            fam->get_workingset(eta, y, v, w, z, no);

            wls_base(alm0, almc, alpha, m, no, ni, X, z, w, intr, ju, vp, cl,
                     nx, thr, maxit, a, &aint, g, ia, iy, &iz, mm, &nino, &rsqc,
                     &nlp, xv, &jerr);
            
            assert(jerr == 0);

            X->compute_eta(eta, a, aint, has_offset, offset);
            current_dev = fam->get_deviance(y, eta, v, no);

            double rel_diff =
                (prev_dev - current_dev) / fmax(fabs(prev_dev), 1.0);
            prev_dev = current_dev;
            if (fabs(rel_diff) < 1e-5) {
                break;
            }

        }
        double devratio = 1 - current_dev / nulldev;
        Rprintf("devratio is %f\n", devratio);
        if ((devratio > 0.999) || (nino > nx)) {
            break;
        }
        Rprintf("number of pass is %d\n", nlp);
        Rprintf("nino is %d\n", nino);

        alm0 = almc;
        // for (int l = 0; l < ni; ++l) {
        //     Rprintf("beta[%d] is %f\n", l + 1, a[l]);
        // }
        // Rprintf("Currently nlp is %d\n", nlp);
        // Rprintf("no is  %d\n", no);
        // Rprintf("ni is  %d\n", ni);
        // Rprintf("ia is  %d\n", ia[0]);
        // Rprintf("iy is  %d\n", iy[ni - 1]);
        // Rprintf("mm is  %d\n", mm[ni - 1]);
        // Rprintf("nino is  %d\n", nino);
        // Rprintf("iz is  %d\n", iz);
        // Rprintf("nlp is  %d\n", nlp);
        // Rprintf("jerr is  %d\n", jerr);
    }

    // Rprintf("no is  %d\n", no);
    // Rprintf("ni is  %d\n", ni);
    // Rprintf("ia is  %d\n", ia[0]);
    // Rprintf("iy is  %d\n", iy[ni-1]);
    // Rprintf("mm is  %d\n", mm[ni-1]);
    // Rprintf("nino is  %d\n", nino);
    // Rprintf("iz is  %d\n", iz);
    // Rprintf("nlp is  %d\n", nlp);
    // Rprintf("jerr is  %d\n", jerr);

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
    free(lambdas);

    return;
}