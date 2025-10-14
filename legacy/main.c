#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <windows.h>
#include "mat.h"


#define TIME_TEST(function_call) ({               \
    clock_t start = clock();                      \
    function_call;                                \
    ((double)(clock() - start)) / CLOCKS_PER_SEC; \
})
#define TEST(T, function_call) ({ \
    V = mat_create_1d_##T(N);  \
    S = mat_create_2d_##T(N);  \
    mat_fill(V);               \
    mat_fill(S);               \
    function_call;             \
    mat_delete(V);             \
    mat_delete(S);             \
})


#define DEFINE_test(T)                           \
void test_##T(int N) {           \
    Mat1D_##T* V, *R = mat_create_1d_##T(N);     \
    Mat2D_##T* S;                                \
                                                 \
    TEST(T,                                      \
        printf(                                  \
            "[%s     N=%5d]  V * S  -  %fs.\n",  \
        #T, N, TIME_TEST(mat_mul(V, S, R)));     \
    );                                           \
                                                 \
    TEST(T,                                      \
        printf(                                  \
            "[%s     N=%5d]  S * V  -  %fs.\n",  \
        #T, N, TIME_TEST(mat_mul(S, V, R)));     \
    );                                           \
                                                 \
    TEST(T,                                      \
        printf(                                  \
            "[%s opt N=%5d]  V * S  -  %fs.\n",  \
        #T, N, TIME_TEST(mat_mul_opt(V, S, R))); \
    );                                           \
                                                 \
    TEST(T,                                      \
        printf(                                  \
            "[%s opt N=%5d]  S * V  -  %fs.\n",  \
        #T, N, TIME_TEST(mat_mul_opt(S, V, R))); \
    );                                           \
                                                 \
    mat_delete(R);                               \
}
#define DEFINE_test_print(T)                                                              \
void test_print_##T(int N) {                                                              \
    Mat1D_##T* V, *R = mat_create_1d_##T(N);                                              \
    Mat2D_##T* S;                                                                         \
                                                                                          \
    TEST(T,                                                                               \
        printf("V:\n");                                                                   \
        mat_print(V);                                                                     \
                                                                                          \
        printf("S:\n");                                                                   \
        mat_print(S);                                                                     \
    );                                                                                    \
                                                                                          \
    TEST(T,                                                                               \
        printf(                                                                           \
            "[%s     N=%5d]  V * S  -  %fs. R:\n", #T, N, TIME_TEST(mat_mul(V, S, R))     \
        );                                                                                \
        mat_print(R);                                                                     \
    );                                                                                    \
                                                                                          \
    TEST(T,                                                                               \
        printf(                                                                           \
            "[%s     N=%5d]  S * V  -  %fs. R:\n", #T, N, TIME_TEST(mat_mul(S, V, R))     \
        );                                                                                \
        mat_print(R);                                                                     \
    );                                                                                    \
                                                                                          \
    TEST(T,                                                                               \
        printf(                                                                           \
            "[%s opt N=%5d]  V * S  -  %fs. R:\n", #T, N, TIME_TEST(mat_mul_opt(V, S, R)) \
        );                                                                                \
        mat_print(R);                                                                     \
    );                                                                                    \
                                                                                          \
    TEST(T,                                                                               \
        printf(                                                                           \
            "[%s opt N=%5d]  S * V  -  %fs. R:\n", #T, N, TIME_TEST(mat_mul_opt(S, V, R)) \
        );                                                                                \
        mat_print(R);                                                                     \
    );                                                                                    \
                                                                                          \
    mat_delete(R);                                                                        \
}

#define DEFINE_start(T)    \
void start_##T(size_t N) { \
    if (N <= 0)            \
        return;            \
    else if (N <= 16)      \
        test_print_##T(N); \
    else                   \
        test_##T(N);       \
}


APPLY_FOR_TYPES(DEFINE_test)
APPLY_FOR_TYPES(DEFINE_test_print)
APPLY_FOR_TYPES(DEFINE_start)


int main() {
    omp_set_num_threads(32);
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);

    size_t N;
    
    printf("Enter N: ");
    scanf("%d", &N);

    start_int(N);
    start_float(N);
    start_double(N);

    return 0;
}
