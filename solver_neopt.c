/*
 * Tema 2 ASC
 * 2019 Spring
 * Catalin Olaru / Vlad Spoiala
 */
#include "utils.h"

/* Returns a newly allocated matrix representing the transpose
 * of matrix A, of size N x N.
 */
static double* get_transpose(double *A, int N) {
	int i, j;
	double *result = malloc(N * N * sizeof(double));

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			result[j * N + i] = A[i * N + j];

	return result;
}

/* Multiply 2 matrices of size N x N the silliy way.
 * Returns a newly allocated matrix.
 */
static double* multiply(double *A, double *B, int N) {
	int i, j, k;
	double *result = malloc(N * N * sizeof(double));

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j) {
			result[i * N + j] = 0.0;
			for (k = 0; k < N; ++k) {
				double temp = A[i * N + k] * B[k * N + j];
				result[i * N + j] += temp;
			}
		}

	return result;
}

/* Add 2 matrices of size N x N. Returns newly allocated matrix. */
static double* add(double *A, double *B, int N) {
	int i, j;
	double *result = malloc(N * N * sizeof(double));

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			result[i * N + j] = A[i * N + j] + B[i * N + j];

	return result;
}

/* Set elements below first diagonal to 0. */
static double* zerotr(double *A, int N) {
	int i, j;
	double *result = malloc(N * N * sizeof(double));

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			result[i * N + j] = j < i ? 0 : A[i * N + j];

	return result;
}

double* my_solver(int N, double *A, double* B) {
	/* Get transposes of factors. */
	double *transpose_A = get_transpose(A, N);
	double *transpose_B = get_transpose(B, N);

	/* Compute the tranpose(A) * B and transpose(B) * A. */
	double *first_multiplication = multiply(transpose_A, B, N);
	double *second_multiplication = multiply(transpose_B, A, N);

	/* Add previous 2 results together. */
	double *addition = add(first_multiplication, second_multiplication, N);
	/* Set bottom half to 0. */
	double *upper_triangle = zerotr(addition, N);

	/* Raise to the power of 2 previous result. */
	double *final_result = multiply(upper_triangle, upper_triangle, N);

	/* Free everything. */
	free(transpose_A);
	free(transpose_B);
	free(first_multiplication);
	free(second_multiplication);
	free(addition);;
	free(upper_triangle);

	return final_result;
}
