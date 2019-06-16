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
        int N = atoi(argv[1]);
        double h = A / N;
        double S = 0.0;

        MPI_Comm_rank(MPI_COMM_WORLD, &k);
        MPI_Comm_size(MPI_COMM_WORLD, &p);

        double d = ceil(log2(p));
        int mask = 0;
        int partner = 0;
        int kapow = 0;
	double t = MPI_Wtime();
        for (i = k; i < N; i+=p) {
                S += h * (f(h * i) + f(h * (i + 1))) / 2.0;
        }

        double sum = S;
        for ( i = 0; i < d; i++) {
                kapow = pow(2, i);
                if ((k & mask) == 0) {
                        partner = k ^ kapow;
                        if (partner < p) {
                                if ((k & kapow) != 0){
                                        MPI_Send(&sum, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD);}
                                else{
                                        MPI_Recv(&S, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                        sum += S;}
                                }
                        }
        mask = mask ^ kapow;
        }

        //double Sk;
        //MPI_Reduce(&S, &Sk, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	t = MPI_Wtime() - t;
        if (k == 0) {
                printf("%f\n", sum);
        }
        MPI_Finalize();
        return 0;
}

