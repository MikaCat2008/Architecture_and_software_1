#define GENDEF_mat_func(T, f)                        \
    Mat1D_##T*: f##_1d_##T,                          \
    Mat2D_##T*: f##_2d_##T 
#define APPLY_FOR_TYPES(f) f(int) f(float) f(double)
#define APPLY_FOR_TYPES_VA_ARGS(f, ...)              \
    f(int, __VA_ARGS__),                             \
    f(float, __VA_ARGS__),                           \
    f(double, __VA_ARGS__)
#define APPLY_FOR_MAT_FUNCTION(M, f) _Generic((M),   \
    APPLY_FOR_TYPES_VA_ARGS(GENDEF_mat_func, f),     \
    default: NULL                                    \
)(M)


#define GENDEF_Mat1D(T)           \
    typedef struct Mat1D_##T {    \
        size_t size;              \
        T* data;                  \
    } Mat1D_##T;          
#define GENDEF_Mat2D(T)           \
    typedef struct Mat2D_##T {    \
        size_t size;              \
        T** data;                 \
    } Mat2D_##T;
#define GENDEF_MAT_TYPES()        \
    APPLY_FOR_TYPES(GENDEF_Mat1D) \
    APPLY_FOR_TYPES(GENDEF_Mat2D)


#define mat_fill(M) APPLY_FOR_MAT_FUNCTION(M, mat_fill)
#define mat_tran(M) _Generic((M),  \
    Mat2D_int*: mat_tran_int,       \
    Mat2D_float*: mat_tran_float,   \
    Mat2D_double*: mat_tran_double, \
    default: NULL                  \
)(M)
#define mat_print(M) APPLY_FOR_MAT_FUNCTION(M, mat_print)
#define mat_delete(M) APPLY_FOR_MAT_FUNCTION(M, mat_delete)


#define GENDEF_mat_mul(T, m, B)                              \
    Mat2D_##T*: _Generic((B),                                \
        Mat1D_##T*: m##_2d_1d_##T,                           \
        default: NULL                                        \
    ),                                                       \
    Mat1D_##T*: _Generic((B),                                \
        Mat2D_##T*: m##_1d_2d_##T,                           \
        default: NULL                                        \
    )                                                                      
#define mat_mul(A, B, R) _Generic((A),                       \
    APPLY_FOR_TYPES_VA_ARGS(GENDEF_mat_mul, mat_mul, B),     \
    default: NULL                                            \
)(A, B, R)
#define mat_mul_opt(A, B, R) _Generic((A),                   \
    APPLY_FOR_TYPES_VA_ARGS(GENDEF_mat_mul, mat_mul_opt, B), \
    default: NULL                                            \
)(A, B, R)


#define DEFINE_mat_create(T)                   \
    Mat1D_##T* mat_create_1d_##T(size_t size); \
    Mat2D_##T* mat_create_2d_##T(size_t size);

#define DEFINE_mat_fill(T)                \
    void mat_fill_1d_##T(Mat1D_##T* mat); \
    void mat_fill_2d_##T(Mat2D_##T* mat);

#define DEFINE_mat_tran(T)             \
    void mat_tran_##T(Mat2D_##T* mat);

#define DEFINE_mat_print(T)                \
    void mat_print_1d_##T(Mat1D_##T* mat); \
    void mat_print_2d_##T(Mat2D_##T* mat);

#define DEFINE_mat_delete(T)                \
    void mat_delete_1d_##T(Mat1D_##T* mat); \
    void mat_delete_2d_##T(Mat2D_##T* mat);

#define DEFINE_mat_mul(T)                                                         \
    void mat_mul_1d_2d_##T(Mat1D_##T* mat_a, Mat2D_##T* mat_b, Mat1D_##T* mat_r); \
    void mat_mul_2d_1d_##T(Mat2D_##T* mat_a, Mat1D_##T* mat_b, Mat1D_##T* mat_r);

#define DEFINE_mat_mul_opt(T)                                                         \
    void mat_mul_opt_1d_2d_##T(Mat1D_##T* mat_a, Mat2D_##T* mat_b, Mat1D_##T* mat_r); \
    void mat_mul_opt_2d_1d_##T(Mat2D_##T* mat_a, Mat1D_##T* mat_b, Mat1D_##T* mat_r);


GENDEF_MAT_TYPES()
APPLY_FOR_TYPES(DEFINE_mat_create)
APPLY_FOR_TYPES(DEFINE_mat_fill)
APPLY_FOR_TYPES(DEFINE_mat_tran)
APPLY_FOR_TYPES(DEFINE_mat_print)
APPLY_FOR_TYPES(DEFINE_mat_delete)
APPLY_FOR_TYPES(DEFINE_mat_mul)
APPLY_FOR_TYPES(DEFINE_mat_mul_opt)
