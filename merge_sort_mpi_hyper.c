/*
 * Николай Игоревич Хохлов, k_h@inbox.ru, 2011-2014.
 * Реализация алгоритма сортировки слиянием.
 * gcc -O2 merge_sort.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

void merge(int *a, int *b, int *c, int na, int nb);
void merge_sort(int *a, int na);
void print(int *a, int na);
int check_sort(int *a, int n);
double timer();

int main (int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
        int size, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int i, n, *a;
	double t;
	int *b;
	int *nk;

	int d;
	d = ceil(log2(size));
	int mask = 0;
	int partner;
	int j;
	if (argc < 2) {
		printf("Usage: %s num_elements.\n", argv[0]);
		return 1;
	}
	n = atoi(argv[1]);
	a = (int*)malloc(sizeof(int) * n);
	srand(time(NULL));

	for (i = 0; i < n; i++) a[i] = rand() % 1000;

	if (rank == 0) {
		if (n < 101) print(a, n);
	}

	nk = (int*)malloc(sizeof(int) * size);
        for(i = 0; i < size-1; i++) {
                nk[i] = n / size;
        }

        nk[size-1] = n - ((size-1) * (n/size));

	t = MPI_Wtime();
	merge_sort(a + rank * (n/size), nk[rank]);
	for(i = 0; i < d; i++) {
		if ((rank & mask) == 0) {
			partner = (rank ^ (1 << i));
			if (partner < size) {
				if ((rank & (1 << i)) != 0) {
					MPI_Send(a + rank * (n/size), nk[rank], MPI_INT, partner, 0, MPI_COMM_WORLD);
				}
				else {
					MPI_Recv(a + i * (n/size), nk[i], MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					b = (int*)malloc(sizeof(int) * (nk[rank] + nk[i]));
					merge(a + rank * (n/size), a + i * (n/size), b, nk[rank], nk[i]);
					nk[rank] = nk[rank] + nk[i];
						for(j = 0; j < nk[rank]; j++) {
							a[j+rank*(n/size)] = b[j];
						}
					free(b);
				}
			}
		}
	mask = (mask ^ (1 << i));
	}

	t = MPI_Wtime() - t;
	if (rank == 0) {
		if (n < 101) print(a, n);
		printf("%d %f\n", size, t);
	}
	free(a);
	MPI_Finalize();
	return 0;
}

double timer()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + 1e-6 * (double)ts.tv_usec;
}

int check_sort(int *a, int n)
{
	int i;
	for (i = 0; i < n - 1; i++) {
		if (a[i] > a[i+1]) {
			return 0;
		}
	}
	return 1;
}

void print(int *a, int na)
{
	int i;
	for (i = 0; i < na; i++) printf("%d ", a[i]);
	printf("\n");
}

/*
 * Процедура слияния массивов a и b в массив c.
 */
void merge(int *a, int *b, int *c, int na, int nb)
{
	/* Реализовать процедуру слияния. */
	int i = 0;
	int j = 0;
	while((i < na) && (j < nb)) {
		if (a[i] <= b[j]) {
			c[i+j] = a[i];
			i += 1;
		}
		else {
			c[i+j] = b[j];
			j += 1;
		}
	}
	if (i < na) {
		while(i < na) {
			c[i+j] = a[i];
			i += 1;
		}
	}
	else {
		while(j < nb) {
                        c[i+j] = b[j];
                        j += 1;
		}
	}
}

/*
 * Процедура сортировки слиянием.
 */
void merge_sort(int *a, int na)
{
	/* Реализовать процедуру сортировки, используя процедуру слияния. */
	if (na < 2) {
		return;
	}
	int i;
	int *b;
	b = (int*)malloc(sizeof(int) * na);
	merge_sort(a, (na/2));
	merge_sort(a + (na/2), na - (na/2));
	merge(a, a + (na/2), b, (na/2), na - (na/2));
	for(i = 0; i < na; i++) {
		a[i] = b[i];
	}
	free(b);
}
