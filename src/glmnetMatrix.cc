#include "glmnetMatrix.h"

#include <math.h>

#include <cstdlib>


MatrixGlmnet::~MatrixGlmnet() {}

int MatrixGlmnet::get_ni() { return this->ni; }
int MatrixGlmnet::get_no() { return this->no; }

double MatrixGlmnet::max_grad(const double* r, const int* ju, const double* vp)
{
    double result = 0.0;
    for(int i = 0; i < ni; ++i){
        if((!ju[i]) || (vp[i] <= 0.0)){
            continue;
        }
        result =fmax(result, fabs(this->dot_product(i, r)));
    }
    return result;
}

double MatrixGlmnet::sumv(const double *v, int len) {
    double result = 0.0;
    for (int i = 0; i < len; ++i) {
        result += v[i];
    }
    return result;
}

DenseM::DenseM(int no, int ni, const double *x) {
    this->no = no;
    this->ni = ni;
    data = x;
}
DenseM::~DenseM() { data = nullptr; }

double DenseM::dot_product(int j, const double *v) {
    // return std::inner_product(data + j * no, data + (j + 1) * no, v, 0.0);

    double result = 0.0;
// If there's no auto vectorization then we can do #pragma clang loop vectorize(enable) interleave(enable)
    for (int i = 0; i < no; ++i) {
        result += data[j * no + i] * v[i];
    }
    return result;
}

double DenseM::vx2(int j, const double *v) {
    double result = 0.0;
    for (int i = 0; i < no; ++i) {
        result += data[j * no + i] * data[j * no + i] * v[i];
    }
    return result;
}

void DenseM::update_res(int j, double d, const double *v,
                        double *__restrict r) {
    for (int i = 0; i < no; ++i) {
        r[i] -= d * v[i] * data[j * no + i];
    }
}

void DenseM::compute_eta(double *__restrict eta, const double *a, double aint,
                         bool has_offset, const double *offset) {
    for (int i = 0; i < no; ++i) {
        eta[i] = aint;
    }
    for (int j = 0; j < ni; ++j) {
        double aj = a[j];
        for (int i = 0; i < no; ++i) {
            eta[i] += data[j * no + i] * aj;
        }
    }
    if (has_offset) {
        for (int i = 0; i < no; ++i) {
            eta[i] += offset[i];
        }
    }
}