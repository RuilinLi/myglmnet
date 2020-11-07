#include <math.h>
#include "glmnetMatrix.h"

// Computation heavily adapted from
// https://github.com/chrchang/plink-ng/blob/master/2.0/pgenlib_ffi_support.cc
inline uintptr_t DivUp(uintptr_t val, uint32_t divisor) {
    return (val + divisor - 1) / divisor;
}
static const uintptr_t kMask5555 = (~static_cast<uintptr_t>(0)) / 3;

static const int32_t kBitsPerWord = 64;
static const int32_t kBitsPerWordD2 = kBitsPerWord / 2;
static const int32_t kCacheline = 64;
static const int32_t kBytesPerWord = kBitsPerWord / 8;
static const int32_t kWordsPerCacheline = kCacheline / kBytesPerWord;
static const int32_t kNypsPerCacheline = kCacheline * 4;

#ifdef _WIN64
inline uint32_t ctzw(unsigned long long ullii) {
    return __builtin_ctzll(ullii);
}
#else
inline uint32_t ctzw(unsigned long ulii) { return __builtin_ctzl(ulii); }
#endif

void GetWeightsByValueNoDosage(const double* weights, const uintptr_t* genoarr,
                               uint32_t sample_ct, double* buf) {
    const uint32_t word_ct = DivUp(sample_ct, kBitsPerWordD2);
    double result = 0.0;
    double result2 = 0.0;
    double miss_weight = 0.0;
    for (uint32_t widx = 0; widx != word_ct; ++widx) {
        const uintptr_t geno_word = genoarr[widx];
        if (!geno_word) {
            continue;
        }
        const double* cur_weights = &(weights[widx * kBitsPerWordD2]);
        uintptr_t geno_word1 = geno_word & kMask5555;
        uintptr_t geno_word2 = (geno_word >> 1) & kMask5555;
        uintptr_t geno_missing_word = geno_word1 & geno_word2;
        geno_word1 ^= geno_missing_word;
        while (geno_word1) {
            const uint32_t sample_idx_lowbits = ctzw(geno_word1) / 2;
            result += cur_weights[sample_idx_lowbits];
            geno_word1 &= geno_word1 - 1;
        }
        geno_word2 ^= geno_missing_word;
        while (geno_word2) {
            const uint32_t sample_idx_lowbits = ctzw(geno_word2) / 2;
            result2 += cur_weights[sample_idx_lowbits];
            geno_word2 &= geno_word2 - 1;
        }
        while (geno_missing_word) {
            const uint32_t sample_idx_lowbits = ctzw(geno_missing_word) / 2;
            miss_weight += cur_weights[sample_idx_lowbits];
            geno_missing_word &= geno_missing_word - 1;
        }
    }
    buf[0] = result;
    buf[1] = result2;
    buf[2] = miss_weight;
}

PlinkMatrix::PlinkMatrix(int no, int ni, const uintptr_t* x,
                         const double* xim) {
    this->no = no;
    this->ni = ni;
    data = x;
    this->xim = xim;
    const uint32_t cache_line_ct = DivUp(no, kNypsPerCacheline);
    word_ct = kWordsPerCacheline * cache_line_ct;
}
PlinkMatrix::~PlinkMatrix() {
    data = nullptr;
    xim = nullptr;
}

double PlinkMatrix::dot_product(int j, const double* v) {
    // assert((j > 0) && (j < ni));

    double buf[3];
    GetWeightsByValueNoDosage(v, &(data[j * word_ct]), no, buf);
    double result = buf[0] + 2 * buf[1] + buf[2] * xim[j];
    return result;
}

double PlinkMatrix::vx2(int j, const double* v) {
    double buf[3];
    GetWeightsByValueNoDosage(v, &(data[j * word_ct]), no, buf);
    double result = buf[0] + 4 * buf[1] + buf[2] * xim[j] * xim[j];
    return result;
}

void PlinkMatrix::update_res(int j, double d, const double* weights,
                             double* r) {
    const uintptr_t* genoarr = &(data[j * word_ct]);

    const uint32_t word_ct_local = DivUp(no, kBitsPerWordD2);
    for (uint32_t widx = 0; widx != word_ct_local; ++widx) {
        const uintptr_t geno_word = genoarr[widx];
        if (!geno_word) {
            continue;
        }
        const double* cur_weights = &(weights[widx * kBitsPerWordD2]);
        uintptr_t geno_word1 = geno_word & kMask5555;
        uintptr_t geno_word2 = (geno_word >> 1) & kMask5555;
        uintptr_t geno_missing_word = geno_word1 & geno_word2;
        geno_word1 ^= geno_missing_word;
        while (geno_word1) {
            const uint32_t sample_idx_lowbits = ctzw(geno_word1) / 2;
            r[widx * kBitsPerWordD2 + sample_idx_lowbits] -=
                d * cur_weights[sample_idx_lowbits];
            geno_word1 &= geno_word1 - 1;
        }
        geno_word2 ^= geno_missing_word;
        while (geno_word2) {
            const uint32_t sample_idx_lowbits = ctzw(geno_word2) / 2;
            r[widx * kBitsPerWordD2 + sample_idx_lowbits] -=
                2 * d * cur_weights[sample_idx_lowbits];
            geno_word2 &= geno_word2 - 1;
        }
        while (geno_missing_word) {
            const uint32_t sample_idx_lowbits = ctzw(geno_missing_word) / 2;
            r[widx * kBitsPerWordD2 + sample_idx_lowbits] -=
                xim[j] * d * cur_weights[sample_idx_lowbits];
            geno_missing_word &= geno_missing_word - 1;
        }
    }
    return;
}

void PlinkMatrix::compute_eta(double* eta, const double* weights, double aint,
                              bool has_offset, const double* offset) {
    for (int i = 0; i < no; ++i) {
        eta[i] = aint;
    }

    // This is only useful for small no, maybe I should just not use it at all?
    const uint32_t word_ct_local = DivUp(no, kBitsPerWordD2);
    for (int j = 0; j < ni; ++j) {
        const uintptr_t* genoarr = &(data[j * word_ct]);
        double ximpute = xim[j];
        double wj = weights[j];

        for (uint32_t widx = 0; widx != word_ct_local; ++widx) {
            const uintptr_t geno_word = genoarr[widx];
            if (!geno_word) {
                continue;
            }
            uintptr_t geno_word1 = geno_word & kMask5555;
            uintptr_t geno_word2 = (geno_word >> 1) & kMask5555;
            uintptr_t geno_missing_word = geno_word1 & geno_word2;
            geno_word1 ^= geno_missing_word;
            while (geno_word1) {
                const uint32_t sample_idx_lowbits = ctzw(geno_word1) / 2;
                eta[widx * kBitsPerWordD2 + sample_idx_lowbits] += wj;
                geno_word1 &= geno_word1 - 1;
            }
            geno_word2 ^= geno_missing_word;
            while (geno_word2) {
                const uint32_t sample_idx_lowbits = ctzw(geno_word2) / 2;
                eta[widx * kBitsPerWordD2 + sample_idx_lowbits] += 2 * wj;
                geno_word2 &= geno_word2 - 1;
            }
            while (geno_missing_word) {
                const uint32_t sample_idx_lowbits = ctzw(geno_missing_word) / 2;
                eta[widx * kBitsPerWordD2 + sample_idx_lowbits] += ximpute * wj;
                geno_missing_word &= geno_missing_word - 1;
            }
        }
    }
    if(has_offset) {
        for(int i = 0; i < no; ++i) {
            eta[i] += offset[i];
        }
    }
}