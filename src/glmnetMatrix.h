#ifndef GLMNET_MATRIX
#define GLMNET_MATRIX
#include <cstdlib>


class MatrixGlmnet {
   public:
    // Compute the inner product of X[,j] and v
    virtual double dot_product(int j, const double* v) = 0;
    // Compute the inner product of X[,i] and X[,j]
    virtual double column_product(int i, int j) = 0;
    // Compute the inner product of X[,j]^2 and v
    virtual double vx2(int j, const double* v) = 0;

    // Set r = r - d*v*x[,j]
    virtual void update_res(int j, double d, const double* v, double* r) = 0;

    static double sumv(const double* v, int len);
    virtual ~MatrixGlmnet();

    int get_no(){return this->no;}
    int get_ni(){return this->ni;}

   protected:
    int no;  // Number of rows
    int ni;  // Number of variables
};

class DenseM : public MatrixGlmnet {
   public:
    DenseM(int no, int ni, const double* x);
    ~DenseM();

    double dot_product(int j, const double* v);

    double column_product(int i, int j);

    double vx2(int j, const double* v);

    void update_res(int j, double d, const double* v, double* r);

   private:
    const double* data;
};


#endif