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
	int size; /* Число процессов. */ //Запускаем такое, что size - квадрат целого числа
	int rank, rankx, ranky;
	int sx, fx, sy, fy;
	int sizex, sizey; 	
	MPI_Comm cart_comm;
	MPI_Comm line_comm;
	MPI_Comm col_comm;
	int n;
	int *bA;
	int *x;
} grid_t;

void multiply_matrix(int n, int sx, int fx, int sy, int fy, int* bA, int *x);
void init(grid_t *g, int n);
void print_matrix(int n, int m, int *A);
void print_vector(int n, int *x);
void decomposition(int n, int p, int k, int *sk, int *fk);

int main(int argc, char **argv)
{
	FILE *out = fopen("ab3.dat", "a");
	MPI_Init(&argc, &argv);
	grid_t g;
	init(&g, 4);
	int *x_new = (int *)malloc(sizeof(int)*g.n); 
	int *tmp = (int *)malloc(sizeof(int)*g.n);
	int i, j, sx, fx;
	int *counts, *displs;
	counts = (int *) malloc(sizeof(int)*g.sizex);
	displs = (int *) malloc(sizeof(int)*g.sizex);
	for (i = 0; i < g.sizex; i++) {
		decomposition(g.n, g.sizex, i, &sx, &fx);
		counts[i] = fx - sx;
		displs[i] = sx;
	}
	MPI_Comm_split(MPI_COMM_WORLD, g.rankx, g.ranky, &g.line_comm);
	MPI_Comm_split(MPI_COMM_WORLD, g.ranky, g.rankx, &g.col_comm);
	double time = MPI_Wtime();
	for (i = 0; i < 3; i++) {
		multiply_matrix(g.n, g.sx, g.fx, g.sy, g.fy, g.bA, g.x);
		MPI_Reduce(g.x+displs[g.rankx], x_new+g.sx, counts[g.rankx], MPI_INT, MPI_SUM, g.rankx, g.line_comm);
		MPI_Bcast(x_new+displs[g.ranky], counts[g.ranky], MPI_INT, g.ranky, g.col_comm);
		tmp = x_new;
		x_new = g.x;
		g.x = tmp;
	}
	MPI_Gatherv(g.x+g.sy, g.fy-g.sy, MPI_INT, g.x, counts, displs, MPI_INT, 0, g.line_comm); //Собираем результат
	time = MPI_Wtime() - time;
	if (g.rank==0) {	
		print_vector(g.n, g.x);
		fprintf(out, "%d\t%lf\n", g.size, time);
	}
	free(tmp);
	free(x_new);
	free(counts);
	free(displs);
	MPI_Finalize();
	fclose(out);
	return 0;
}


void multiply_matrix(int n, int sx, int fx, int sy, int fy, int *bA, int *x)
{
	int i, j;
	int *x_new = (int *) calloc(sizeof(int), n);
	for (i = 0; i < fx-sx; i++)
		for (j = 0; j < fy-sy; j++)
			x_new[i+sx] += bA[i*(fy-sy)+j]*x[j+sy];
	for (i = 0; i < n; i++) x[i] = x_new[i];
	free(x_new);
}

void print_matrix(int n, int m, int *A)
{
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++)
			printf("%d ", A[i * m + j]);
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
	int dims[2] = {0, 0};
	MPI_Dims_create(g->size, 2, dims);
	int periods[2] = {0, 0};
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &(g->cart_comm));
	g->sizex = dims[0];
	g->sizey = dims[1];
	int coords[2]; 
	MPI_Cart_coords(g->cart_comm, g->rank, 2, coords);
	g->rankx = coords[0];
	g->ranky = coords[1];
	g->n = n;
	decomposition(g->n, g->sizey, g->ranky, &g->sy, &g->fy);
	decomposition(g->n, g->sizex, g->rankx, &g->sx, &g->fx);	
	g->bA = (int*)malloc(sizeof(int) * (g->fx - g->sx) * (g->fy - g->sy));
	g->x = (int*)malloc(sizeof(int) * n);
	int i, j;
	for (i = 0; i < n; i++) g->x[i] = 1;
	for (i = 0; i < g->fx-g->sx; i++) {
		for (j = 0; j < g->fy-g->sy; j++) {
			g->bA[i * (g->fy-g->sy) + j] = (i+g->sx) * n + j+g->sy;//(i+g->sk+j) % 2;
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
