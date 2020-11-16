#include <assert.h>
#include <math.h>
#include <string.h>

#include "glmfamily.h"
#include "glmnetMatrix.h"

void wls_base(double alm0, double almc, double alpha, int m, int no, int ni,
              MatrixGlmnet *X, double *r, const double *v, int intr,
              const int *ju, const double *vp, const double *cl, int nx,
              double thr, int maxit, double *__restrict a, double *aint,
              double *__restrict g, int *__restrict ia, int *__restrict iy,
              int *iz, int *__restrict mm, int *nino, double *rsqc, int *nlp,
              double *__restrict xv, double * xm,int *jerr, int irls_iter, double rsum, double vsum);

void glmnetPath(double alpha, MatrixGlmnet *X, const double *y, const double *v,
                int intr, const int *ju, const double *vp, const double *cl,
                int nx, double thr, int maxit, GlmFamily *fam, bool has_offset,
                double *offset, const double *ulambdas, int nlambda, int mxitnr,
                const double flmin, int *lmu, double *a0, double *ca, int *ia,
                int *nin, double *devratio_vec, double *alm, int *nlp,
                double *nulldev_ptr, int *jerr, double *a, int *iy, int *mm, int nino, int warm, double *residuals) {
    int no = X->get_no();
    int ni = X->get_ni();
    double aint = 0;  // intercept

    double *w = (double *)malloc(sizeof(double) * no);
    double *z = (double *)malloc(sizeof(double) * no);
    double *eta = (double *)malloc(sizeof(double) * no);

    double *g = (double *)malloc(sizeof(double) * ni);
    double *xv = (double *)malloc(sizeof(double) * ni);

    // xv stores sum(v * x * x), xm store sum(v * x)
    double *xm = (double *)malloc(sizeof(double) * ni);



    int iz = warm;
    double rsqc = 0;

    double nulldev =
        fam->null_deviance(y, v, w, intr, eta, has_offset, offset, &aint, no);
    *nulldev_ptr = nulldev;

    if(warm) {
        X->compute_eta(eta, a, aint, has_offset, offset);
    }


    // Compute max_lambda here instead of using user defined lambdas
    double *lambdas = (double *)malloc(sizeof(double) * nlambda);
    if (flmin < 1.0) {
        // Note that, for plink matrix the mean is not removed when computing max_lambda
        // this is find since the residual is orthogonal to all-ones vector
        double max_lambda = X->max_grad(w, ju, vp);
        max_lambda /= fmax(alpha, 1e-3);
        lambdas = (double *)malloc(sizeof(double) * nlambda);
        lambdas[0] = max_lambda;
        double ratio = pow(flmin, 1.0 / (nlambda - 1));
        for (int i = 1; i < nlambda; ++i) {
            lambdas[i] = lambdas[i - 1] * ratio;
        }
    } else {
        for (int i = 0; i < nlambda; ++i) {
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
            double sumbuf[2];
            fam->get_workingset(eta, y, v, w, z, no, sumbuf);
            double rsum = sumbuf[0];
            double vsum = sumbuf[1];

            wls_base(alm0, almc, alpha, m, no, ni, X, z, w, intr, ju, vp, cl,
                     nx, thr, maxit, a, &aint, g, ia, iy, &iz, mm, &nino, &rsqc,
                     nlp, xv, xm, jerr, j, rsum, vsum);

            assert((*jerr) == 0);

            X->compute_eta(eta, a, aint, has_offset, offset);
            current_dev = fam->get_deviance(y, eta, v, no);

            double rel_diff =
                (prev_dev - current_dev) / fmax(fabs(prev_dev), 1.0);
            prev_dev = current_dev;
            if (fabs(rel_diff) < 1e-4) {
                break;
            }
        }
        double devratio = 1 - current_dev / nulldev;

        // Copy data to output
        for (int k = 0; k < nino; ++k) {
            ca[m * nx + k] = a[ia[k]];
        }
        fam->get_residual(y, eta, v, residuals+m*no, no);

        a0[m] = aint;
        nin[m] = nino;

        devratio_vec[m] = devratio;
        alm[m] = almc;
        (*lmu)++;
        if ((devratio > 0.999) || (nino > nx)) {
            break;
        }

        if((m > 0) &&  ((devratio - devratio_vec[m-1]) < devratio*1e-5)){
            break;
        }


        alm0 = almc;

    }


    free(xm);
    free(g);
    free(w);
    free(z);
    free(eta);
    free(xv);
    free(lambdas);

    return;
}