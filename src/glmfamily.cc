#include "glmfamily.h"

#include <math.h>


GlmFamily::~GlmFamily() {}
Gaussian::Gaussian() {}

// Copy these vectors, though maybe only need to do it once?
void Gaussian::get_workingset(const double *eta, const double *y,
                              const double *v, double *w, double *z, int len, double* sumbuf) {
    sumbuf[0] = 0;
    sumbuf[1] = 0;
    for (int i = 0; i < len; ++i) {
        w[i] = v[i];
        z[i] = v[i] * (y[i] - eta[i]);
        sumbuf[0] += z[i];
        sumbuf[1] += w[i];
    }
}

double Gaussian::get_deviance(const double *y, const double *eta,
                              const double *v, int len) {
    double result = 0.0;
    for (int i = 0; i < len; ++i) {
        result += (y[i] - eta[i]) * (y[i] - eta[i]) * v[i];
    }
    return result;
}

void Gaussian::get_residual(const double *y, const double *eta, const double *v,
                            double *r, int len) {
    for (int i = 0; i < len; ++i) {
        r[i] = v[i] * (y[i] - eta[i]);
    }
}

double Gaussian::null_deviance(const double *y, const double *v, double *r,
                               int intr, double *eta, bool has_offset,
                               const double *offset, double *aint, int len) {
    if (!has_offset) {
        double weightysquare = 0;
        double weighty = 0;
        double weightsum = 0;
        for (int i = 0; i < len; ++i) {
            weightysquare += v[i] * y[i] * y[i];
            weighty += v[i] * y[i];
            weightsum += v[i];
        }
        if (!intr) {
            for (int i = 0; i < len; ++i) {
                eta[i] = 0;
                r[i] = v[i] * y[i];
            }
            return weightysquare;
        }

        double beta0 = weighty / weightsum;
        for (int i = 0; i < len; ++i) {
            eta[i] = beta0;
            r[i] = v[i] * (y[i] - beta0);
        }
        *aint = beta0;
        return (weightysquare - 2 * weighty * beta0 +
                weightsum * beta0 * beta0);
    }

    double weightysquare = 0;
    double weighty = 0;
    double weightsum = 0;
    for (int i = 0; i < len; ++i) {
        double y0 = y[i] - offset[i];
        weightysquare += v[i] * y0 * y0;
        weighty += v[i] * y0;
        weightsum += v[i];
    }
    if (!intr) {
        for (int i = 0; i < len; ++i) {
            eta[i] = offset[i];
            r[i] = v[i] * (y[i] - offset[i]);
        }
        return weightysquare;
    }
    double beta0 = weighty / weightsum;
    for (int i = 0; i < len; ++i) {
        eta[i] = offset[i] + beta0;
        r[i] = v[i] * (y[i] - beta0 - offset[i]);
    }
    return (weightysquare - 2 * weighty * beta0 + weightsum * beta0 * beta0);
}

Logistic::Logistic() {}
void Logistic::get_workingset(const double *eta, const double *y,
                              const double *v, double *w, double *z, int len, double* sumbuf) {
    sumbuf[0] = 0;
    sumbuf[1] = 0;
    for (int i = 0; i < len; ++i) {
        double p = 1 / (1 + exp(-eta[i]));
        w[i] = p * (1 - p) * v[i];
        z[i] = (y[i] - p) * v[i];
        sumbuf[0] += z[i];
        sumbuf[1] += w[i];
    }
}

double Logistic::get_deviance(const double *y, const double *eta,
                              const double *v, int len) {
    double result = 0;
    for (int i = 0; i < len; ++i) {
        result += v[i] * (log(1 + exp(eta[i])) - y[i] * eta[i]);
    }
    result *= 2;
    return result;
}

void Logistic::get_residual(const double *y, const double *eta, const double *v,
                            double *r, int len) {
    for (int i = 0; i < len; ++i) {
        double eeta = exp(eta[i]);
        r[i] = v[i] * (y[i] - eeta / (1 + eeta));
    }
}

double Logistic::null_deviance(const double *y, const double *v, double *r,
                               int intr, double *eta, bool has_offset,
                               const double *offset, double *aint, int len) {
    if (!has_offset) {
        double count0 = 0;
        double count1 = 0;
        for (int i = 0; i < len; ++i) {
            count0 += v[i] * y[i];
            count1 += v[i] * (1 - y[i]);
        }
        if (!intr) {
            for (int i = 0; i < len; ++i) {
                eta[i] = 0;
                r[i] = v[i] * (y[i] - 0.5);
            }
            return (2 * log(2) * (count0 + count1));
        }
        double p0inv = 1 + (count1 / count0);
        double p0 = 1/p0inv;
        double aint_val = -log(p0inv - 1);
        *aint = aint_val;
        for (int i = 0; i < len; ++i) {
            eta[i] = aint_val;
            r[i] = v[i] * (y[i] - p0);
        }
        return (2 * log(p0inv) * (count0 + count1));
    }

    throw "Offset has not been implemented for logstic regression\n";

    // double count0 = 0;
    // double count1 = 0;
    // for(int i = 0; i < len; ++i) {
    //     count0 += w[i] * y[i] * log(1 + exp(-offset[i]));
    //     count1 += w[i] * (1-y[i]) * log(1 + exp(offset[i]));
    // }
    // if(!intr){
    //     return (2*(count0 + count1));
    // }
}
