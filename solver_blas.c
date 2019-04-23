/*
 * Tema 2 ASC
 * 2019 Spring
 * Catalin Olaru / Vlad Spoiala
 */
#include <cblas.h>

#include "utils.h"
/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
	double *first_multiplication = malloc(N * N * sizeof(double));
	double *second_multiplication = malloc(N * N * sizeof(double));
	double *final_result = malloc(N * N * sizeof(double));
	double alpha = 1.0;
	double beta = 0.0;
	double *identity_matrix = malloc(N * N * sizeof(double));
	int i, j;

/* Compute identity matrix so we can use it later on. */
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			identity_matrix[i * N + j] = j == i ? 1.0 : 0.0;

/* Compute Transpose(A) * B in first_multiplication. */
	cblas_dgemm(CblasRowMajor,
				CblasTrans,
				CblasNoTrans,
				N,
				N,
				N,
				alpha,
				A,
				N,
				B,
				N,
				beta,
				first_multiplication,
				N);
/* Compute Transpose(B) * A in second_multiplication. */
	cblas_dgemm(CblasRowMajor,
				CblasTrans,
				CblasNoTrans,
				N,
				N,
				N,
				alpha,
				B,
				N,
				A,
				N,
				beta,
				second_multiplication,
				N);
/* Add the previous 2 results in second multiplication. */
	cblas_dgemm(CblasRowMajor,
				CblasNoTrans,
				CblasNoTrans,
				N,
				N,
				N,
				alpha,
				first_multiplication,
				N,
				identity_matrix,
				N,
				alpha,
				second_multiplication,
				N);

/* Apply zerotr(second_multiplication) as everything involving
 * CblasLower or CblasUpper didn't work for some reason. :(
 */
	for (i = 0; i < N; ++i)
		for (j = 0; j < i; ++j)
			second_multiplication[i * N + j] = 0;

/* Multiply second_multiplication by itself so we get the
 * result in final_result. */
	cblas_dgemm(CblasRowMajor,
				CblasNoTrans,
				CblasNoTrans,
				N,
				N,
				N,
				alpha,
				second_multiplication,
				N,
				second_multiplication,
				N,
				beta,
				final_result,
				N);

/* Free everything aside from the return value. */
	free(first_multiplication);
	free(second_multiplication);
	free(identity_matrix);

	return final_result;
}
