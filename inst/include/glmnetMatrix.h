#ifndef GLMNET_MATRIX
#define GLMNET_MATRIX

#include <stdint.h>
class MatrixGlmnet {
   public:
    // Compute the inner product of X[,j] and v
    virtual double dot_product(int j, const double* v) = 0;

    // Compute the inner product of X[,j]^2 and v
    virtual double vx2(int j, const double* v) = 0;

    // Set r = r - d*v*x[,j]
    virtual void update_res(int j, double d, const double* v, double* r) = 0;

    // Set eta = X * a + aint + offset, offset is optional
    virtual void compute_eta(double* eta, const double* a, double aint,
                             bool has_offset, const double* offset) = 0;

    // Compute weighted mean and standard deviation of each variable,
    // ignore variables with ju[j] == 0
    // if variable is constant, set ju[j] = 0
    // virtual void get_xmxs(const double* v, const int* ju, double* xm,
    //                       double* xs) = 0;

    static double sumv(const double* v, int len);
    virtual ~MatrixGlmnet();

    int get_no();
    int get_ni();

    double max_grad(const double* r, const int* ju, const double* vp);

   protected:
    int no;  // Number of rows
    int ni;  // Number of variables
};

class DenseM : public MatrixGlmnet {
   public:
    DenseM(int no, int ni, const double* x);
    ~DenseM();

    double dot_product(int j, const double* v);

    double vx2(int j, const double* v);

    void update_res(int j, double d, const double* v, double* r);

    void compute_eta(double* eta, const double* a, double aint, bool has_offset,
                     const double* offset);

   private:
    const double* data;
};

class PlinkMatrix : public MatrixGlmnet {
   public:
    PlinkMatrix(int no, int ni, const uintptr_t* x, const double* xim);
    ~PlinkMatrix();

    double dot_product(int j, const double* v);

    double vx2(int j, const double* v);

    void update_res(int j, double d, const double* v, double* r);

    void compute_eta(double* eta, const double* a, double aint, bool has_offset,
                     const double* offset);
    
   private:
     const uintptr_t* data;
     uint32_t word_ct;
     const double *xim; // mean imputation
};

#endif