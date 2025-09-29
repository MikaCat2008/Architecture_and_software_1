#include <stdlib.h>


typedef struct
{
    size_t size;
    double* data;
} Mat1Dd;

typedef struct
{
    size_t size;
    double* data;
} Mat2Dd;


Mat1Dd* mat_1dd(size_t size);

Mat2Dd* mat_2dd(size_t size);

void mat_fill_1dd(Mat1Dd* mat);

void mat_fill_2dd(Mat2Dd* mat);

void mat_mul_1dd_2dd(Mat1Dd* mat1d, Mat2Dd* mat2d, Mat1Dd* result_mat);

void mat_mul_2dd_1dd(Mat2Dd* mat2d, Mat1Dd* mat1d, Mat1Dd* result_mat);

void mat_print_1dd(Mat1Dd* mat);

void mat_print_2dd(Mat2Dd* mat);

void mat_delete_1dd(Mat1Dd* mat);

void mat_delete_2dd(Mat2Dd* mat);
