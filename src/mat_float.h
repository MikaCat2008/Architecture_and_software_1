#include <stdlib.h>


typedef struct
{
    size_t size;
    float* data;
} Mat1Df;

typedef struct
{
    size_t size;
    float* data;
} Mat2Df;


Mat1Df* mat_1df(size_t size);

Mat2Df* mat_2df(size_t size);

void mat_fill_1df(Mat1Df* mat);

void mat_fill_2df(Mat2Df* mat);

void mat_mul_1df_2df(Mat1Df* mat1d, Mat2Df* mat2d, Mat1Df* result_mat);

void mat_mul_2df_1df(Mat2Df* mat2d, Mat1Df* mat1d, Mat1Df* result_mat);

void mat_print_1df(Mat1Df* mat);

void mat_print_2df(Mat2Df* mat);

void mat_delete_1df(Mat1Df* mat);

void mat_delete_2df(Mat2Df* mat);
