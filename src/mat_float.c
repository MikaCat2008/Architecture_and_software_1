#include <stdio.h>
#include <stdlib.h>
#include "mat_float.h"


Mat1Df* mat_1df(size_t size)
{
    Mat1Df* mat = (Mat1Df*)malloc(sizeof(Mat1Df));
    mat->size = size;
    mat->data = _aligned_malloc(size * sizeof(float), 16);

    return mat;
}

Mat2Df* mat_2df(size_t size)
{
    Mat2Df* mat = (Mat2Df*)malloc(sizeof(Mat2Df));

    mat->size = size;
    mat->data = _aligned_malloc(size * size * sizeof(float), 16);

    return mat;
}

void mat_fill_1df(Mat1Df* mat) 
{
    size_t size = mat->size;
    float* data = mat->data;
    
    for (int i = 0; i < size; i++)
        data[i] = i % 10;
}

void mat_fill_2df(Mat2Df* mat) 
{
    size_t size = mat->size;
    float* data = mat->data;

    for (int i = 0; i < size * size; i++)
        data[i] = i % 10;
}

void mat_mul_1df_2df(Mat1Df* mat1d, Mat2Df* mat2d, Mat1Df* result_mat) 
{
    int jd;
    size_t size = mat1d->size; 
    float result;

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

void mat_mul_2df_1df(Mat2Df* mat2d, Mat1Df* mat1d, Mat1Df* result_mat) 
{
    int id; 
    size_t size = mat1d->size;
    float result;

    for (int i = 0; i < size; i++)
    {
        id = i * size;
        result = 0;

        for (int j = 0; j < size; j++)
            result += mat2d->data[id + j] * mat1d->data[j];

        result_mat->data[i] = result;
    }
}

void mat_print_1df(Mat1Df* mat)
{
    size_t size = mat->size;
    float* data = mat->data;

    printf("[");

    for (int i = 0; i < size; i++)
        printf(" %5.1f", data[i]);

    printf(" ]\n");
}

void mat_print_2df(Mat2Df* mat)
{
    size_t size = mat->size;
    float* data = mat->data;

    printf("[[");
    
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) 
            printf(" %5.1f", data[i * size + j]);

        if (i + 1 < size)
            printf("\n  ");
    }
    
    printf(" ]]\n");
}

void mat_delete_1df(Mat1Df* mat) 
{
    _aligned_free(mat->data);
    free(mat);
}

void mat_delete_2df(Mat2Df* mat) 
{
    _aligned_free(mat->data);
    free(mat);
}
