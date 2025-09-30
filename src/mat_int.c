#include <stdio.h>
#include <stdlib.h>
#include "mat_int.h"


Mat1Di* mat_1di(size_t size)
{
    Mat1Di* mat = (Mat1Di*)malloc(sizeof(Mat1Di));
    mat->size = size;
    mat->data = _aligned_malloc(size * sizeof(int), 16);

    return mat;
}

Mat2Di* mat_2di(size_t size)
{
    Mat2Di* mat = (Mat2Di*)malloc(sizeof(Mat2Di));

    mat->size = size;
    mat->data = _aligned_malloc(size * size * sizeof(int), 16);

    return mat;
}

void mat_fill_1di(Mat1Di* mat) 
{
    size_t size = mat->size;
    int* data = mat->data;
    
    for (int i = 0; i < size; i++)
        data[i] = i % 10;
}

void mat_fill_2di(Mat2Di* mat) 
{
    size_t size = mat->size;
    int* data = mat->data;

    for (int i = 0; i < size * size; i++)
        data[i] = i % 10;
}

void mat_mul_1di_2di(Mat1Di* mat1d, Mat2Di* mat2d, Mat1Di* result_mat) 
{
    int jd;
    size_t size = mat1d->size; 
    int result;

    for (int i = 0; i < size; i++)
    {
        jd = i;
        result = 0;

        for (int j = 0; j < size; j++) 
        {
            result += mat1d->data[j] * mat2d->data[jd];
            jd += size;
        }

        result_mat->data[i] = result;
    }
}

void mat_mul_2di_1di(Mat2Di* mat2d, Mat1Di* mat1d, Mat1Di* result_mat) 
{
    int id; 
    size_t size = mat1d->size;
    int result;

    for (int i = 0; i < size; i++)
    {
        id = i * size;
        result = 0;

        for (int j = 0; j < size; j++)
            result += mat2d->data[id + j] * mat1d->data[j];

        result_mat->data[i] = result;
    }
}

void mat_print_1di(Mat1Di* mat)
{
    size_t size = mat->size;
    int* data = mat->data;

    printf("[");

    for (int i = 0; i < size; i++)
        printf(" %5d", data[i]);

    printf(" ]\n");
}

void mat_print_2di(Mat2Di* mat)
{
    size_t size = mat->size;
    int* data = mat->data;

    printf("[[");
    
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) 
            printf(" %5d", data[i * size + j]);

        if (i + 1 < size)
            printf("\n  ");
    }
    
    printf(" ]]\n");
}

void mat_delete_1di(Mat1Di* mat) 
{
    _aligned_free(mat->data);
    free(mat);
}

void mat_delete_2di(Mat2Di* mat) 
{
    _aligned_free(mat->data);
    free(mat);
}
