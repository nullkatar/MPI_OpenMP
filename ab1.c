/**
 * MPI matrix multiplication.
 * Author: Nikolay I. Khokhlov <k_h@inbox.ru>, 2018.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

typedef struct {
	int size; /* Число процессов. */
	int rank; /* Глобальный номер процесса. */
	int sk, fk;	

	int n;
	int *bA;
	int *x;
} grid_t;

void multiply_matrix(int n, int sk, int fk, int* bA, int *x);
void init(grid_t *g, int n);
void print_matrix(int n, int *A);
void print_vector(int n, int *x);
void decomposition(int n, int p, int k, int *sk, int *fk);

int main(int argc, char **argv)
{
	FILE *out = fopen("ab1.dat", "a");
	MPI_Init(&argc, &argv);
	grid_t g;
	init(&g, 6);
	int i, sk, fk;
	int *counts, *displs;
	counts = (int *) malloc(sizeof(int)*g.size);
	displs = (int *) malloc(sizeof(int)*g.size);
	for (i = 0; i < g.size; i++) {
		decomposition(g.n, g.size, i, &sk, &fk);
		counts[i] = fk - sk;
		displs[i] = sk;
	}
	double time = MPI_Wtime();
	for (i = 0; i < 3; i++) {
		multiply_matrix(g.n, g.sk, g.fk, g.bA, g.x);
		MPI_Allgatherv(g.x+g.sk, g.fk-g.sk, MPI_INT, g.x, counts, displs, MPI_INT, MPI_COMM_WORLD);
	}
	time = MPI_Wtime() - time;
	if (g.rank==0) {	
		//print_matrix(g.n, g.bA);
		//print_vector(g.n, g.x);
		fprintf(out, "%d\t%lf\n", g.size, time);
	}
	MPI_Finalize();
	fclose(out);
	free(counts);
	free(displs);
	return 0;
}


void multiply_matrix(int n, int sk, int fk, int *bA, int *x)
{
	int i, j;
	int *x_new = (int *) calloc(sizeof(int), n);
	for (i = 0; i < fk - sk; i++)
		for (j = 0; j < n; j++)
			x_new[i+sk] += bA[i*n+j]*x[j];
	for (i = 0; i < n; i++) x[i] = x_new[i];
	free(x_new);
}

void print_matrix(int n, int *A)
{
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%d ", A[i * n + j]);
		printf("\n");
	}
	printf("\n");
}

void print_vector(int n, int *x)
{
	int i;
	for (i = 0; i < n; i++) printf("%d\n", x[i]);
	printf("\n");
}

void init(grid_t *g, int n) {
	MPI_Comm_size(MPI_COMM_WORLD, &g->size);
	MPI_Comm_rank(MPI_COMM_WORLD, &g->rank);	
	g->n = n;
	decomposition(n, g->size, g->rank, &g->sk, &g->fk);
	g->bA = (int*)malloc(sizeof(int) * (g->fk - g->sk) * n);
	g->x = (int*)malloc(sizeof(int) * n);
	int i, j;
	for (i = 0; i < n; i++) g->x[i] = 1;
	for (i = 0; i < (g->fk - g->sk); i++) {
		for (j = 0; j < n; j++) {
			g->bA[i * n + j] = (i + g->sk) * n + j;//(i+g->sk+j) % 2;
		}
	}
}

void decomposition(int n, int p, int k, int *sk, int *fk) {
	*sk = k * (n / p);
	*fk = *sk + n / p;
	if (k == p - 1) {
		*fk = n;
	}
}
