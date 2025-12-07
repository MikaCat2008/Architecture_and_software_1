#include <time.h>
#include <stdio.h>
#include <typeinfo>

inline double get_time(clock_t start) {
    return ((double)(clock() - start)) / CLOCKS_PER_SEC;
}

template <typename T>
void test(size_t N) {
    const char *type_name = typeid(T).name(), 
               *print_format = type_name[0] == 'i' ? " %5.1d" : " %5.1f";

    T* V = new T[N];
    for (size_t i = 0; i < N; i++)
        V[i] = (T)(i % 10);

    T* R = new T[N];
    for (size_t i = 0; i < N; i++)
        R[i] = (T)0;
        
    T** S = new T*[N];
    for (size_t i = 0; i < N; i++) {
        S[i] = new T[N];

        for (size_t j = 0; j < N; j++)
            S[i][j] = (T)((i * N + j) % 10);
    }

    if (N <= 16) {
        printf("V:\n");

        for (size_t i = 0; i < N; i++)
            printf(print_format, V[i]);
        
        printf("\nS:\n");

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++)
                printf(print_format, S[i][j]);

            printf("\n");
        }
    }

    clock_t start = clock();
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            R[i] += S[i][j] * V[j];
    printf("[%s     N=%5d]  S * V  -  %fs.\n", type_name, N, get_time(start));

    if (N <= 16) 
        for (size_t i = 0; i < N; i++)
            printf(print_format, R[i]);
    for (size_t i = 0; i < N; i++)
        R[i] = 0;

    start = clock();
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            R[i] += V[j] * S[j][i];
    printf("\n[%s     N=%5d]  V * S  -  %fs.\n", type_name, N, get_time(start));

    if (N <= 16) 
        for (size_t i = 0; i < N; i++)
            printf(print_format, R[i]);
    for (int i = 0; i < N; i++)
        R[i] = 0;

    clock_t start = clock();
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            R[j] += V[i] * S[i][j];
    printf("\n[%s opt N=%5d]  V * S  -  %fs.\n", type_name, N, get_time(start));

    if (N <= 16) 
        for (size_t i = 0; i < N; i++)
            printf(print_format, R[i]);

    for (size_t i = 0; i < N; i++)
        delete[] S[i];
    delete[] V, R, S;
}

int main() {
    size_t N;

    printf("Enter N: ");
    scanf("%d", &N);

    test<int>(N);
    test<float>(N);
    test<double>(N);
    
    return 0;
}