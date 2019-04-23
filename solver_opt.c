/*
 * Tema 2 ASC
 * 2019 Spring
 * Catalin Olaru / Vlad Spoiala
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 */

/* Block size used to compute multiplication faster.
 * MIN returns the smallest of 2 numbers.
 */
#define block_size 40
#define MIN(a, b) (((a) < (b) ? (a) : (b)))

/* Returns a newly allocated matrix representing the transpose
 * of the input of size N x N.
 */
static double* get_transpose(double *A, int N) {
	register int i, j;
	double *result = malloc(N * N * sizeof(double));

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			*(result + j * N + i) = *(A + i * N + j);

	return result;
}

/* Returns a newly allocated matrix computing
 * A * B (both of size N x N).
 */
static double* multiply(double *A, double *B, int N) {
	register int i, j, k;
	double *result = calloc(N * N, sizeof(double));

	for (i = 0; i < N; ++i)
		for (k = 0; k < N; ++k) {
			for (j = 0; j < N; ++j)
				*(result + i * N + j) += *(A + i * N + k) * *(B + k * N + j);
		}

	return result;
}

/* Adds 2 matrices into the first one. */
static void add(double *A, double *B, int N) {
	register int i, j;

	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			*(A + i * N + j) += *(B + i * N + j);
}

/* Set lower bottom half of matrix to 0. */
static void zerotr(double *A, int N) {
	register int i, j;

	for (i = 0; i < N; ++i)
		for (j = 0; j < i; ++j)
			*(A + i * N + j) = 0;
}

double* my_solver(int N, double *A, double* B) {
	double *transpose_A = get_transpose(A, N);
	double *transpose_B = get_transpose(B, N);

	double *first_multiplication = multiply(transpose_A, B, N);
	double *second_multiplication = multiply(transpose_B, A, N);

	add(first_multiplication, second_multiplication, N);
	zerotr(first_multiplication, N);

	double *final_result = multiply(first_multiplication, first_multiplication, N);

	free(transpose_A);
	free(transpose_B);
	free(first_multiplication);
	free(second_multiplication);

	return final_result;
}

