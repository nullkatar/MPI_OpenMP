/**
 * MPI matrix multiplication.
 * Author: Nikolay I. Khokhlov <k_h@inbox.ru>, 2018.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define N 4

typedef struct {
	int size; /* Число процессов. */
	int rank; /* Глобальный номер процесса. */
	int q;
	MPI_Comm grid_comm; /**/
	MPI_Comm row_comm; /**/
	MPI_Comm col_comm; /**/
	int row;    // Coordinates of our processes in the row_communicator
	int col;    // Same thing for column
	
	int n;
	int ln;
	int *bA;
	int *bB;
	int *bC;
	
	/*
	 * Gather/Scatter
	 */
	MPI_Datatype blocktype;
	int *counts;
	int *disps;
} grid_t;

void init_grid(grid_t *g);
void multiply_matrix(int n, int* bA, int* bB, int* bC);
void init_matrix(grid_t *g, int n);
void gather_matrix(grid_t *g, int n, int *A);
void fox(grid_t *g);
void cannon(grid_t *g);
void print_matrix(int n, int *A);

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	grid_t g;
	init_grid(&g);
	init_matrix(&g, N);
	double time = MPI_Wtime();
	fox(&g);
	//cannon(&g);
	//print_matrix(g.ln, g.bC);
	int *C;
	if (g.rank == 0) {
		C = (int*)malloc(sizeof(int) * g.n * g.n);
	}
	gather_matrix(&g, g.n, C);
	if (g.rank==0) {	
		print_matrix(g.n, C);
	}
	MPI_Finalize();
	return 0;
}

void fox(grid_t *g)
{
	int s, k, i, j;
	int *tmpA = (int *)malloc(sizeof(int)*g->ln*g->ln);
	for (s = 0; s < g->q; s++) {
		k = (g->row + s) % g->q;
		if (k == g->col) {
			MPI_Bcast(g->bA, g->ln*g->ln, MPI_INT, k, g->row_comm);
			multiply_matrix(g->ln, g->bA, g->bB, g->bC);
		}
		else {
			MPI_Bcast(tmpA, g->ln*g->ln, MPI_INT, k, g->row_comm);
			multiply_matrix(g->ln, tmpA, g->bB, g->bC);
		}
		if (s != g->q-1) ring_shift(g->bB, g->ln*g->ln, g->col_comm, -1);
	}
	free(tmpA);	
}

void ring_shift(int *data, size_t count, MPI_Comm ring, int delta)
{
	int src, dst;
	MPI_Cart_shift(ring, 0, delta, &src, &dst);
	MPI_Sendrecv_replace(data, count, MPI_INT, dst, 0, src, 0, ring, MPI_STATUS_IGNORE);
}

void cannon(grid_t *g)
{
	int s;
	ring_shift(g->bA, g->ln*g->ln, g->row_comm, -g->row);
	ring_shift(g->bB, g->ln*g->ln, g->col_comm, -g->col);
	for (s = 0; s < g->q; s++) {
		multiply_matrix(g->ln, g->bA, g->bB, g->bC);
		ring_shift(g->bA, g->ln*g->ln, g->row_comm, 1);	
		ring_shift(g->bB, g->ln*g->ln, g->col_comm, 1);
	}
}


void multiply_matrix(int n, int* bA, int* bB, int* bC)
{
	int i, j, k;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				bC[i*n+j] += bA[i*n+k]*bB[k*n+j];
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

void init_matrix(grid_t *g, int n)
{
	g->n  = n;
	g->ln = n / g->q;
	g->bA = (int*)malloc(sizeof(int) * g->ln * g->ln);
	g->bB = (int*)malloc(sizeof(int) * g->ln * g->ln);
	g->bC = (int*)calloc(g->ln * g->ln, sizeof(int));
	int *A;
	int *B;
	int i, j;
	if (g->rank == 0) {
		A = (int*)malloc(sizeof(int) * n * n);
		B = (int*)malloc(sizeof(int) * n * n);
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				A[i * n + j] = i * n + j;
				B[i * n + j] = i * n + j;
			}
		}
	}

	MPI_Datatype blocktype2;
	MPI_Type_vector(g->ln, g->ln, g->n, MPI_INT, &blocktype2);
	MPI_Type_create_resized(blocktype2, 0, sizeof(int), &g->blocktype);
	MPI_Type_commit(&g->blocktype);

	g->disps = (int*)malloc(sizeof(int) * g->size);
	g->counts = (int*)malloc(sizeof(int) * g->size);
	for (i=0; i < g->q; i++) {
		for (j=0; j < g->q; j++) {
			g->disps[i*g->q+j] = i*n*g->ln+j*g->ln;
			g->counts[i*g->q+j] = 1;
		}
	}
	MPI_Scatterv(A, g->counts, g->disps, g->blocktype, g->bA, g->ln * g->ln, MPI_INT, 0, g->grid_comm);
	MPI_Scatterv(B, g->counts, g->disps, g->blocktype, g->bB, g->ln * g->ln, MPI_INT, 0, g->grid_comm);
}

void gather_matrix(grid_t *g, int n, int *A)
{
	MPI_Gatherv(g->bC, g->ln * g->ln, MPI_INT, A, g->counts, g->disps, g->blocktype, 0, g->grid_comm);
}

void init_grid(grid_t *g) {
	int rank;
	int dims[2];
	int period[2];
	int coords[2];
	int free_coords[2];
	MPI_Comm_size(MPI_COMM_WORLD, &g->size);
	g->q = dims[0] = dims[1] = (int)sqrt(g->size);
	period[0] = period[1] = 1;

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 0, &g->grid_comm);
	MPI_Comm_rank(g->grid_comm, &g->rank);
	MPI_Cart_coords(g->grid_comm, g->rank, 2, coords);
	g->row = coords[0];
	g->col = coords[1];

	free_coords[0] = 0; 
	free_coords[1] = 1;
	MPI_Cart_sub(g->grid_comm, free_coords, &g->row_comm);

	free_coords[0] = 1; 
	free_coords[1] = 0;
	MPI_Cart_sub(g->grid_comm, free_coords, &g->col_comm);
}
