#include <time.h>
#include <stdio.h>
#include <windows.h>
#include <intrin.h>
#include <stddef.h>
#include "mat.h"

#define OPT_BLOCK(T) _Generic(((T*)0), \
    int*: 32,                          \
    float*: 32,                        \
    double*: 16,                       \
    default: 32\                       
)
#define VALLOC(size)                                                   \
    VirtualAlloc(NULL, size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE)
#define VFREE(ptr)                   \
    VirtualFree(ptr, 0, MEM_RELEASE)
#define FORMAT_MAT_ELEMENT(T) _Generic(((T*)0), \
    int*: " %5.1d",                             \
    float*: " %5.1f",                           \
    double*: " %5.1f",                          \
    default: ""                                 \
)


#define DEFINE_mat_create(T)                                  \
    Mat1D_##T* mat_create_1d_##T(size_t size) {               \
        Mat1D_##T* M = (Mat1D_##T*)VALLOC(sizeof(Mat1D_##T)); \
                                                              \
        M->size = size;                                       \
        M->data = (T*)VALLOC(size * sizeof(T));               \
                                                              \
        return M;                                             \
    }                                                         \
    Mat2D_##T* mat_create_2d_##T(size_t size) {               \
        Mat2D_##T* M = (Mat2D_##T*)VALLOC(sizeof(Mat2D_##T)); \
                                                              \
        M->size = size;                                       \
        M->data = (T**)VALLOC(size * sizeof(T*));             \
                                                              \
        for (size_t i = 0; i < size; i++)                     \
            M->data[i] = (T*)VALLOC(size * sizeof(T));        \
                                                              \
        return M;                                             \
    }

#define DEFINE_mat_fill(T)                           \
    void mat_fill_1d_##T(Mat1D_##T* M) {             \
        for (size_t i = 0; i < M->size; i++)         \
            M->data[i] = i % 10;                     \
    }                                                \
    void mat_fill_2d_##T(Mat2D_##T* M) {             \
        size_t size = M->size;                       \
                                                     \
        for (size_t i = 0; i < size; i++)            \
            for (size_t j = 0; j < size; j++)        \
                M->data[i][j] = (i * size + j) % 10; \
    }

#define DEFINE_mat_tran(T)                                                          \
    void mat_tran_##T(Mat2D_##T* A) {                                               \
        T* Ai;                                                                      \
        size_t size = A->size;                                                      \
                                                                                    \
        _Pragma("omp parallel for collapse(2)")                                     \
        for (size_t i = 0; i < size; i += OPT_BLOCK(T)) {                           \
            for (size_t j = i; j < size; j += OPT_BLOCK(T)) {                       \
                size_t i_max = (i + OPT_BLOCK(T) < size) ? i + OPT_BLOCK(T) : size; \
                size_t j_max = (j + OPT_BLOCK(T) < size) ? j + OPT_BLOCK(T) : size; \
                                                                                    \
                if (i == j) {                                                       \
                    for (size_t ii = i; ii < i_max; ++ii) {                         \
                        Ai = A->data[ii];                                           \
                                                                                    \
                        for (size_t jj = ii + 1; jj < j_max; ++jj) {                \
                            T tmp = Ai[jj];                                         \
                            Ai[jj] = A->data[jj][ii];                               \
                            A->data[jj][ii] = tmp;                                  \
                        }                                                           \
                    }                                                               \
                } else {                                                            \
                    for (size_t ii = i; ii < i_max; ++ii) {                         \
                        Ai = A->data[ii];                                           \
                                                                                    \
                        for (size_t jj = j; jj < j_max; ++jj) {                     \
                            T tmp = Ai[jj];                                         \
                            Ai[jj] = A->data[jj][ii];                               \
                            A->data[jj][ii] = tmp;                                  \
                        }                                                           \
                    }                                                               \
                }                                                                   \
            }                                                                       \
        }                                                                           \
    }

#define DEFINE_mat_print(T)                                   \
    void mat_print_1d_##T(Mat1D_##T* M) {                     \
        printf("[");                                          \
                                                              \
        for (size_t i = 0; i < M->size; i++)                  \
            printf(FORMAT_MAT_ELEMENT(T), M->data[i]);        \
                                                              \
        printf(" ]\n");                                       \
    }                                                         \
    void mat_print_2d_##T(Mat2D_##T* M) {                     \
        size_t size = M->size;                                \
                                                              \
        printf("[[");                                         \
                                                              \
        for (size_t i = 0; i < size; i++) {                   \
            for (size_t j = 0; j < size; j++)                 \
                printf(FORMAT_MAT_ELEMENT(T), M->data[i][j]); \
                                                              \
            if (i + 1 < size)                                 \
                printf("\n  ");                               \
        }                                                     \
                                                              \
        printf(" ]]\n");                                      \
    }

#define DEFINE_mat_delete(T)                 \
    void mat_delete_1d_##T(Mat1D_##T* M) {   \
        VFREE(M->data);                      \
        VFREE(M);                            \
    }                                        \
    void mat_delete_2d_##T(Mat2D_##T* M) {   \
        for (size_t i = 0; i < M->size; i++) \
            VFREE(M->data[i]);               \
                                             \
        VFREE(M->data);                      \
        VFREE(M);                            \
    }

#define DEFINE_mat_mul(T)                                              \
    void mat_mul_1d_2d_##T(Mat1D_##T* A, Mat2D_##T* B, Mat1D_##T* R) { \
        T result;                                                      \
        size_t size = A->size;                                         \
                                                                       \
        for (size_t i = 0; i < size; i++) {                            \
            result = (T)0;                                             \
                                                                       \
            for (size_t j = 0; j < size; j++)                          \
                result += A->data[j] * B->data[j][i];                  \
                                                                       \
            R->data[i] = result;                                       \
        }                                                              \
    }                                                                  \
    void mat_mul_2d_1d_##T(Mat2D_##T* A, Mat1D_##T* B, Mat1D_##T* R) { \
        T result;                                                      \
        size_t size = A->size;                                         \
                                                                       \
        for (size_t i = 0; i < size; i++) {                            \
            result = (T)0;                                             \
                                                                       \
            for (size_t j = 0; j < size; j++)                          \
                result += A->data[i][j] * B->data[j];                  \
                                                                       \
            R->data[i] = result;                                       \
        }                                                              \
    }                                                                  \

void mat_mul_opt_1d_2d_int(Mat1D_int* A, Mat2D_int* B, Mat1D_int* R) {
    mat_tran(B);
    mat_mul(B, A, R);
}

void mat_mul_opt_2d_1d_int(Mat2D_int* A, Mat1D_int* B, Mat1D_int* R) {
    int* Ai;
    size_t size = A->size;
    __m128i result;

    #pragma omp parallel for
    for (size_t i = 0; i < size; i++) {
        Ai = A->data[i];
        result = _mm_setzero_si128();

        for (size_t j = 0; j < size; j += 4) {
            __m128i a = _mm_load_si128((__m128i*)&Ai[j]);
            __m128i b = _mm_load_si128((__m128i*)&B->data[j]);
            result = _mm_add_epi32(result, _mm_mullo_epi32(a, b));
        }

        int tmp[4];
        _mm_store_si128((__m128i*)tmp, result);

        R->data[i] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
    }
}

void mat_mul_opt_1d_2d_float(Mat1D_float* A, Mat2D_float* B, Mat1D_float* R) {
    mat_tran(B);
    mat_mul(B, A, R);
    mat_tran(B);
}

void mat_mul_opt_2d_1d_float(Mat2D_float* A, Mat1D_float* B, Mat1D_float* R) { 
    float* Ai;
    size_t size = A->size;
    __m128 result;

    #pragma omp parallel for
    for (size_t i = 0; i < size; i++) {
        Ai = A->data[i];
        result = _mm_setzero_ps();

        for (size_t j = 0; j < size; j += 4) {
            __m128 a = _mm_load_ps(&Ai[j]);
            __m128 b = _mm_load_ps(&B->data[j]);
            result = _mm_add_ps(result, _mm_mul_ps(a, b));
        }

        result = _mm_hadd_ps(result, result);
        result = _mm_hadd_ps(result, result);

        R->data[i] = _mm_cvtss_f32(result);
    }
}

void mat_mul_opt_1d_2d_double(Mat1D_double* A, Mat2D_double* B, Mat1D_double* R) {
    mat_tran(B);
    mat_mul(B, A, R);
}

void mat_mul_opt_2d_1d_double(Mat2D_double* A, Mat1D_double* B, Mat1D_double* R) {
    double* Ai;
    size_t size = A->size;
    __m128d result;

    #pragma omp parallel for
    for (size_t i = 0; i < size; i++) {
        Ai = A->data[i];
        result = _mm_setzero_pd();

        for (size_t j = 0; j < size; j += 2) {
            __m128d a = _mm_load_pd(&Ai[j]);
            __m128d b = _mm_load_pd(&B->data[j]);
            result = _mm_add_pd(result, _mm_mul_pd(a, b));
        }

        R->data[i] = result[0] + result[1];
    }
}


APPLY_FOR_TYPES(DEFINE_mat_create)
APPLY_FOR_TYPES(DEFINE_mat_fill)
APPLY_FOR_TYPES(DEFINE_mat_print)
APPLY_FOR_TYPES(DEFINE_mat_delete)
APPLY_FOR_TYPES(DEFINE_mat_mul)
APPLY_FOR_TYPES(DEFINE_mat_tran)
