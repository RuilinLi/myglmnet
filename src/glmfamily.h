#ifndef GLM_FAMILY
#define GLM_FAMILY

class GlmFamily {
   public:
    // Compute the z and w in a IRLS, z here is actually a weighted working
    // residual
    virtual void get_workingset(const double *eta, const double *y,
                                const double *v, double *w, double *z,
                                int len) = 0;
    virtual double get_deviance(const double *y, const double *eta,
                                const double *v, int len) = 0;

    // Residual is vector such that r^TX is the gradient of the log-likelihood
    virtual void get_residual(const double *y, const double *eta,
                              const double *v, double *r, int len) = 0;

    // This function does three things:
    // 1. Compute the null deviance
    // 2. Get parameter estimate for the intercept, if needed
    // 3. Compute weighted residual of a model with intercept only
    virtual double null_deviance(const double *y, const double *v, double *r,
                                 int intr, double *eta, bool has_offset,
                                 const double *offset, double *aint,
                                 int len) = 0;
    virtual ~GlmFamily();
};

class Gaussian : public GlmFamily {
   public:
    Gaussian();

    // Do nothing
    void get_workingset(const double *eta, const double *y, const double *v,
                        double *w, double *z, int len);

    double get_deviance(const double *y, const double *eta, const double *v,
                        int len);

    void get_residual(const double *y, const double *eta, const double *v,
                      double *r, int len);

    double null_deviance(const double *y, const double *v, double *r, int intr,
                         double *eta, bool has_offset, const double *offset,
                         double *aint, int len);
    static const char *get_name();
};

class Logistic : public GlmFamily {
   public:
    Logistic();
    void get_workingset(const double *eta, const double *y, const double *v,
                        double *w, double *z, int len);

    double get_deviance(const double *y, const double *eta, const double *v,
                        int len);

    void get_residual(const double *y, const double *eta, const double *v,
                      double *r, int len);

    double null_deviance(const double *y, const double *v, double *r, int intr,
                         double *eta, bool has_offset, const double *offset,
                         double *aint, int len);
    static const char *get_name();
};

#endif