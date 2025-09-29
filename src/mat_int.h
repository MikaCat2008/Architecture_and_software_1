#include <stdlib.h>


typedef struct
{
    size_t size;
    int* data;
} Mat1Di;

typedef struct
{
    size_t size;
    int* data;
} Mat2Di;


Mat1Di* mat_1di(size_t size);

Mat2Di* mat_2di(size_t size);

void mat_fill_1di(Mat1Di* mat);

void mat_fill_2di(Mat2Di* mat);

void mat_mul_1di_2di(Mat1Di* mat1d, Mat2Di* mat2d, Mat1Di* result_mat);

void mat_mul_2di_1di(Mat2Di* mat2d, Mat1Di* mat1d, Mat1Di* result_mat);

void mat_print_1di(Mat1Di* mat);

void mat_print_2di(Mat2Di* mat);

void mat_delete_1di(Mat1Di* mat);

void mat_delete_2di(Mat2Di* mat);
