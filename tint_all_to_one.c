#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define A 2.0

double f(double x)
{
	return sqrt(4.0 - x * x);
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		printf("Usage: %s num_intervals.\n", argv[0]);
		return 1;
	}
	MPI_Init(&argc, &argv);
	int k, p;
	int i;
	int j;
	int N = atoi(argv[1]);
	double h = A / N;
	double S = 0.0;
	int d = ceil(log2(p));
        int mask;
        int partner;
        double Sk = S;
        int twoi;
        double t;
	MPI_Comm_rank(MPI_COMM_WORLD, &k);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	t = MPI_Wtime();
	for (j = 0; j < 100000; j++) {
	mask = 0;
	for (i = k; i < N; i+=p) {
		S += h * (f(h * i) + f(h * (i + 1))) / 2.0;
	}
	for (i = 0; i < d; i++) {
		twoi = (1 << i);
		if ((k & mask) == 0) {
			partner = (k ^ twoi);
			if (partner < p) {
			if ((k & twoi) != 0) {
				MPI_Send(&Sk, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD);
			}
			else {
				MPI_Recv(&S, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				Sk += S;
			}
			}
		}
		mask = (mask ^ twoi);
	}
	}
//	MPI_Reduce(&S, &Sk, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	t = MPI_Wtime() - t;
	if (k == 0) {
		printf("%d %f\n", p, t);
	}
	MPI_Finalize();
	return 0;
}
