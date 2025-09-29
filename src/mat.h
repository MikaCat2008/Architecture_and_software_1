#include "mat_int.h"
#include "mat_float.h"
#include "mat_double.h"


#define mat_fill(mat) _Generic((mat), \
    Mat1Di*: mat_fill_1di, \
    Mat2Di*: mat_fill_2di, \
    Mat1Df*: mat_fill_1df, \
    Mat2Df*: mat_fill_2df, \
    Mat1Dd*: mat_fill_1dd, \
    Mat2Dd*: mat_fill_2dd, \
    default: NULL \
)(mat)

#define mat_mul(mat_a, mat_b, mat_r) _Generic((mat_a), \
    Mat1Di*: _Generic((mat_b), \
        Mat2Di*: mat_mul_1di_2di, \
        default: NULL \
    ), \
    Mat2Di*: _Generic((mat_b), \
        Mat1Di*: mat_mul_2di_1di, \
        default: NULL \
    ), \
    Mat1Df*: _Generic((mat_b), \
        Mat2Df*: mat_mul_1df_2df, \
        default: NULL \
    ), \
    Mat2Df*: _Generic((mat_b), \
        Mat1Df*: mat_mul_2df_1df, \
        default: NULL \
    ), \
    Mat1Dd*: _Generic((mat_b), \
        Mat2Dd*: mat_mul_1dd_2dd, \
        default: NULL \
    ), \
    Mat2Dd*: _Generic((mat_b), \
        Mat1Dd*: mat_mul_2dd_1dd, \
        default: NULL \
    ), \
    default: NULL \
)(mat_a, mat_b, mat_r)

#define mat_print(mat) _Generic((mat), \
    Mat1Di*: mat_print_1di, \
    Mat2Di*: mat_print_2di, \
    Mat1Df*: mat_print_1df, \
    Mat2Df*: mat_print_2df, \
    Mat1Dd*: mat_print_1dd, \
    Mat2Dd*: mat_print_2dd, \
    default: NULL \
)(mat)

#define mat_delete(mat) _Generic((mat), \
    Mat1Di*: mat_delete_1di, \
    Mat2Di*: mat_delete_2di, \
    Mat1Df*: mat_delete_1df, \
    Mat2Df*: mat_delete_2df, \
    Mat1Dd*: mat_delete_1dd, \
    Mat2Dd*: mat_delete_2dd, \
    default: NULL \
)(mat)
