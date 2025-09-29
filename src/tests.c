#include <time.h>
#include <stdio.h>

#include "mat.h"

// #define PRINT_MATRIX
#define SHORT_OUTPUT

static int n = 1000;


double get_time(double start_time)
{
    return ((double)(clock() - start_time)) / CLOCKS_PER_SEC;
}


void tests_int()
{
    Mat1Di* result_mat;
    clock_t start_time;

#ifdef SHORT_OUTPUT
    printf("%d ", n);
#else
    printf("\n=== Замір швидкості для типу int ===\n\n");
#endif


    // створення квадратної матриці розміром n на n
    start_time = clock();

    Mat2Di* square_mat = mat_2di(n);
    mat_fill(square_mat);

#ifdef SHORT_OUTPUT
    printf("%f ", get_time(start_time));
#else
    printf("1. створення квадратної матриці розміром %d на %d - %fс:\n", n, n, get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(square_mat);
    #endif
#endif
    //


    // створення вектору розміром n
    start_time = clock();

    Mat1Di* vector_mat = mat_1di(n);
    mat_fill(vector_mat);

#ifdef SHORT_OUTPUT
    printf("%f ", get_time(start_time));
#else
    printf("2. створення вектору розміром %d - %fс\n", n, get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(vector_mat);
    #endif
#endif
    //


    // перемноження квадратної матриці на вектор
    result_mat = mat_1di(n);
    start_time = clock();

    mat_mul(square_mat, vector_mat, result_mat);

#ifdef SHORT_OUTPUT
    printf("%f ", get_time(start_time));
#else
    printf("3. перемноження квадратної матриці на вектор - %fс\n", get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(result_mat);
    #endif
#endif
    mat_delete(result_mat);
    //


    // перемноження квадратної матриці на вектор
    result_mat = mat_1di(n);
    start_time = clock();

    mat_mul(vector_mat, square_mat, result_mat);

#ifdef SHORT_OUTPUT
    printf("%f\n", get_time(start_time));
#else
    printf("4. перемноження вектора на квадратну матрицю - %fс\n", get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(result_mat);
    #endif
#endif
    mat_delete(result_mat);
    //
    

    mat_delete(square_mat);
    mat_delete(vector_mat);
}


void tests_float()
{
    Mat1Df* result_mat;
    clock_t start_time;

#ifdef SHORT_OUTPUT
    printf("%d ", n);
#else
    printf("\n=== Замір швидкості для типу float ===\n\n");
#endif


    // створення квадратної матриці розміром n на n
    start_time = clock();

    Mat2Df* square_mat = mat_2df(n);
    mat_fill(square_mat);

#ifdef SHORT_OUTPUT
    printf("%f ", get_time(start_time));
#else
    printf("1. створення квадратної матриці розміром %d на %d - %fс:\n", n, n, get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(square_mat);
    #endif
#endif
    //


    // створення вектору розміром n
    start_time = clock();

    Mat1Df* vector_mat = mat_1df(n);
    mat_fill(vector_mat);

#ifdef SHORT_OUTPUT
    printf("%f ", get_time(start_time));
#else
    printf("2. створення вектору розміром %d - %fс\n", n, get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(vector_mat);
    #endif
#endif
    //


    // перемноження квадратної матриці на вектор
    result_mat = mat_1df(n);
    start_time = clock();

    mat_mul(square_mat, vector_mat, result_mat);

#ifdef SHORT_OUTPUT
    printf("%f ", get_time(start_time));
#else
    printf("3. перемноження квадратної матриці на вектор - %fс\n", get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(result_mat);
    #endif
#endif
    mat_delete(result_mat);
    //


    // перемноження квадратної матриці на вектор
    result_mat = mat_1df(n);
    start_time = clock();

    mat_mul(vector_mat, square_mat, result_mat);

#ifdef SHORT_OUTPUT
    printf("%f\n", get_time(start_time));
#else
    printf("4. перемноження вектора на квадратну матрицю - %fс\n", get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(result_mat);
    #endif
#endif
    mat_delete(result_mat);
    //
    

    mat_delete(square_mat);
    mat_delete(vector_mat);
}


void tests_double()
{
    Mat1Dd* result_mat;
    clock_t start_time;

#ifdef SHORT_OUTPUT
    printf("%d ", n);
#else
    printf("\n=== Замір швидкості для типу double ===\n\n");
#endif


    // створення квадратної матриці розміром n на n
    start_time = clock();

    Mat2Dd* square_mat = mat_2dd(n);
    mat_fill(square_mat);

#ifdef SHORT_OUTPUT
    printf("%f ", get_time(start_time));
#else
    printf("1. створення квадратної матриці розміром %d на %d - %fс:\n", n, n, get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(square_mat);
    #endif
#endif
    //


    // створення вектору розміром n
    start_time = clock();

    Mat1Dd* vector_mat = mat_1dd(n);
    mat_fill(vector_mat);

#ifdef SHORT_OUTPUT
    printf("%f ", get_time(start_time));
#else
    printf("2. створення вектору розміром %d - %fс\n", n, get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(vector_mat);
    #endif
#endif
    //


    // перемноження квадратної матриці на вектор
    result_mat = mat_1dd(n);
    start_time = clock();

    mat_mul(square_mat, vector_mat, result_mat);

#ifdef SHORT_OUTPUT
    printf("%f ", get_time(start_time));
#else
    printf("3. перемноження квадратної матриці на вектор - %fс\n", get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(result_mat);
    #endif
#endif
    mat_delete(result_mat);
    //


    // перемноження квадратної матриці на вектор
    result_mat = mat_1dd(n);
    start_time = clock();

    mat_mul(vector_mat, square_mat, result_mat);

#ifdef SHORT_OUTPUT
    printf("%f \n", get_time(start_time));
#else
    printf("4. перемноження вектора на квадратну матрицю - %fс\n", get_time(start_time));

    #ifdef PRINT_MATRIX
        mat_print(result_mat);
    #endif
#endif
    mat_delete(result_mat);
    //
    

    mat_delete(square_mat);
    mat_delete(vector_mat);
}


int main() {
    for (int i = 1; i < 47; i++)
    {
        n = i * 1000;

        tests_int();
        tests_float();
        tests_double();
    }

    return 0;
}
