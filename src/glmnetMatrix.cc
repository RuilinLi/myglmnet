#include "glmnetMatrix.h"
#include "R.h"
#include <math.h>

#include <cstdlib>
#include <iostream>
#include <numeric>

MatrixGlmnet::~MatrixGlmnet() {}

int MatrixGlmnet::get_ni() { return this->ni; }
int MatrixGlmnet::get_no() { return this->no; }

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
    return std::inner_product(data + j * no, data + (j + 1) * no, v, 0.0);
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

void DenseM::compute_eta(double *eta, const double *a, double aint,
                         bool has_offset, const double *offset) {
    Rprintf("Aint is %f\n", aint);
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
        Rprintf("should not reach here!\n");
        for (int i = 0; i < no; ++i) {
            eta[i] += offset[i];
        }
    }
}