#include <cmath>
#include <cstdlib>

#include "glmnetMatrix.h"
void wls_base(double alm0, double almc, double alpha, int m, int no, int ni,
              MatrixGlmnet *X, double *r, const double *v, int intr,
              const int *ju, const double *vp, const double *cl, int nx,
              double thr, int maxit, double *__restrict a, double *aint,
              double *__restrict g, int *__restrict ia, int *__restrict iy,
              int *iz, int *__restrict mm, int *nino, double *rsqc, int *nlp, double *__restrict xv,
              int *jerr) {
    //double *__restrict xv = (double *)malloc(sizeof(double) * ni);
    double xmz = MatrixGlmnet::sumv(v, no);
    double ab = almc * alpha;
    double dem = almc * (1.0 - alpha);
    double tlam = alpha * (2.0 * almc - alm0);

    for (int j = 0; j < ni; ++j) {
        if (ju[j]) {
            g[j] = abs(X->dot_product(j, r));

        } else {
            continue;
        }

        if (iy[j]) {
            xv[j] = X->vx2(j, v);
        } else if (g[j] > tlam * vp[j]) {
            iy[j] = 1;
            xv[j] = X->vx2(j, v);
        }
    }

    bool jz = true;

    while (true) {
        if (!((*iz) && jz)) {
            (*nlp)++;
            double dlx = 0.0;
            for (int j = 0; j < ni; ++j) {
                if (!iy[j]) {
                    continue;
                }

                double gj = X->dot_product(j, r);
                double aj = a[j];
                double u = gj + aj * xv[j];
                double au = abs(u) - vp[j] * ab;
                if (au < 0.0) {
                    a[j] = 0.0;
                } else {
                    a[j] = fmax(cl[2 * j],
                                fmin(cl[2 * j + 1],
                                     copysign(au, u) / (xv[j] + vp[j] * dem)));
                }

                if (a[j] == aj) {
                    continue;
                }

                if (mm[j] == 0) {
                    (*nino)++;
                    if ((*nino) > nx) {
                        break;
                    }
                    mm[j] = (*nino);
                    ia[(*nino) - 1] = j;
                }
                double d = a[j] - aj;
                (*rsqc) += d * (2.0 * gj - d * xv[j]);
                X->update_res(j, d, v, r);
                dlx = fmax(xv[j] * d * d, dlx);
            }
            if ((*nino) > nx) {
                break;
            }
            if (intr) {
                double sumr = MatrixGlmnet::sumv(r, no);
                double d = sumr / xmz;
                (*aint) += d;
                (*rsqc) += d * (2.0 * sumr - d * xmz);

                dlx = fmax(dlx, xmz * d * d);

                for (int i = 0; i < no; ++i) {
                    r[i] -= d * v[i];
                }
            }

            // KKT checking here
            if (dlx < thr) {
                bool ixx = false;
                for (int j = 0; j < ni; ++j) {
                    if (iy[j] || (!ju[j])) {
                        continue;
                    }
                    g[j] = abs(X->dot_product(j, r));
                    if (g[j] > ab * vp[j]) {
                        iy[j] = 1;
                        xv[j] = X->vx2(j, v);
                        ixx = true;
                    }
                }

                if (ixx) {
                    continue;
                }
                break;
            }

            if ((*nlp) > maxit) {
                *jerr = -m;
                return;
            }
        }
        (*iz) = 1;

        while (true) {
            (*nlp)++;
            double dlx = 0.0;
            for (int l = 0; l < (*nino); ++l) {
                int k = ia[l];
                double gk = X->dot_product(k, r);
                double ak = a[k];
                double u = gk + ak * xv[k];
                double au = abs(u) - vp[k] * ab;
                if (au < 0.0) {
                    a[k] = 0.0;
                } else {
                    a[k] = fmax(cl[2 * k],
                                fmin(cl[2 * k + 1],
                                     copysign(au, u) / (xv[k] + vp[k] * dem)));
                }

                if (ak == a[k]) {
                    continue;
                }
                double d = a[k] - ak;
                (*rsqc) += d * (2.0 * gk - d * xv[k]);
                X->update_res(k, d, v, r);
                dlx = fmax(xv[k] * d * d, dlx);
            }

            if (intr) {
                double sumr = MatrixGlmnet::sumv(r, no);
                double d = sumr / xmz;
                (*aint) += d;
                (*rsqc) += d * (2.0 * sumr - d * xmz);

                dlx = fmax(dlx, xmz * d * d);

                for (int i = 0; i < no; ++i) {
                    r[i] -= d * v[i];
                }
            }

            if (dlx < thr) {
                break;
            }

            if ((*nlp) > maxit) {
                *jerr = -m;
                return;
            }
        }
        jz = false;
    }
    //free(xv);
    return;
}