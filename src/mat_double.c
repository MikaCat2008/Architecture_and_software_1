#include <stdio.h>
#include <stdlib.h>
#include "mat_double.h"


Mat1Dd* mat_1dd(size_t size)
{
    Mat1Dd* mat = (Mat1Dd*)malloc(sizeof(Mat1Dd));
    mat->size = size;
    mat->data = _aligned_malloc(size * sizeof(double), 16);

    return mat;
}

Mat2Dd* mat_2dd(size_t size)
{
    Mat2Dd* mat = (Mat2Dd*)malloc(sizeof(Mat2Dd));

    mat->size = size;
    mat->data = _aligned_malloc(size * size * sizeof(double), 16);

    return mat;
}

void mat_fill_1dd(Mat1Dd* mat) 
{
    size_t size = mat->size;
    double* data = mat->data;
    
    for (int i = 0; i < size; i++)
        data[i] = i % 10;
}

void mat_fill_2dd(Mat2Dd* mat) 
{
    size_t size = mat->size;
    double* data = mat->data;

    for (int i = 0; i < size * size; i++)
        data[i] = i % 10;
}

void mat_mul_1dd_2dd(Mat1Dd* mat1d, Mat2Dd* mat2d, Mat1Dd* result_mat) 
{
    int jd;
    size_t size = mat1d->size; 
    double result;

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

void mat_mul_2dd_1dd(Mat2Dd* mat2d, Mat1Dd* mat1d, Mat1Dd* result_mat) 
{
    int id; 
    size_t size = mat1d->size;
    double result;

    for (int i = 0; i < size; i++)
    {
        id = i * size;
        result = 0;

        for (int j = 0; j < size; j++)
            result += mat2d->data[id + j] * mat1d->data[j];

        result_mat->data[i] = result;
    }
}

void mat_print_1dd(Mat1Dd* mat)
{
    size_t size = mat->size;
    double* data = mat->data;

    printf("[");

    for (int i = 0; i < size; i++)
        printf(" %5.1f", data[i]);

    printf(" ]\n");
}

void mat_print_2dd(Mat2Dd* mat)
{
    size_t size = mat->size;
    double* data = mat->data;

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

void mat_delete_1dd(Mat1Dd* mat) 
{
    _aligned_free(mat->data);
    free(mat);
}

void mat_delete_2dd(Mat2Dd* mat) 
{
    _aligned_free(mat->data);
    free(mat);
}
