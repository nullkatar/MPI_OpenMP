
#include <mpi.h>
#include <stdio.h>

#define X 1000
#define Y 1000
#define C 1
#define H 0.1
#define DT (H * H / 4 / C)
#define STEPS 10000

typedef float real;
#define REAL MPI_FLOAT

int rank, size;
int from, to;

int rows_count(int my_rank);
void bounds(real *u);
void gather(real *u);

void init(real *u);
void save(real *u, const char *path);
void step(real *u, real *u1);

int main(int argc, char *argv[])
{
	real **u[2];
	for (i = 0; i < 2; i++) {
		u[i] = (real *)malloc(X*Y*sizeof(real));
	}

	int i;
	int z = 0;
	double loop_time;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	from = 0;
	for (i = 0; i < rank; i++) {
		from += rows_count(i);
	}
	to = from + rows_count(rank);
	printf("%d: from = %d, to = %d\n", rank, from, to);
	
	init((real*)u[z]);
	init((real*)u[1-z]);
	
	gather((real*)u[z]);
	if (rank == 0) {
		save((real*)u[z], "in.dat");
	}
	
	loop_time = MPI_Wtime();

	for (i = 0; i < STEPS; i++) {
		bounds((real*)u[z]);
		step((real*)u[z], (real*)u[1-z]);
		z = 1 - z;
	}

	loop_time = MPI_Wtime() - loop_time;
	
	gather((real*)u[z]);
	
	if (rank == 0) {
		printf("%d\t%f\n", size, loop_time);

	}
	MPI_Finalize();

	for (i = 0; i < 2; i++) {
		free(u[i]);
	}
	return 0;
}

void gather(real *u)
{
	int i, c;
	real *ru = u;
	if (rank == 0) {
		ru += Y * rows_count(0);
		for (i = 1; i < size; i++) {
			c = rows_count(i);
			MPI_Recv(ru, c * Y, REAL, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			ru += c * Y;
		}
	} else {
		MPI_Send(u + from * Y, Y * (to-from), REAL, 0, 0, MPI_COMM_WORLD);
	}
}

void bounds(real *u)
{
	int left  = rank - 1;
	int right = rank + 1;
	
	if (rank % 2) {
		if (right < size) {
			MPI_Sendrecv(u+(to-1)*Y, Y, REAL, right, 0, u+to*Y,       Y, REAL, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if (left >= 0) {
			MPI_Sendrecv(u+from*Y,   Y, REAL, left,  0, u+(from-1)*Y, Y, REAL, left,  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		if (left >= 0) {
			MPI_Sendrecv(u+from*Y,   Y, REAL, left,  0, u+(from-1)*Y, Y, REAL, left,  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if (right < size) {
			MPI_Sendrecv(u+(to-1)*Y, Y, REAL, right, 0, u+to*Y,       Y, REAL, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
}


int rows_count(int my_rank)
{
	return X / size + (X % size > my_rank);
}

void init(real *u)
{
	int i, j;
	for (i = from; i < to; i++) {
	for (j = 0; j < Y; j++) {
		if (i > X / 2 - X / 20 && j > Y / 2 - Y / 20 && i < X / 2 + X / 20 && j < Y / 2 + Y / 20) {
			*(u+j+i*Y) = 1.0;
		} else {
			*(u+j+i*Y) = 0.0;
		}
	}
	}
}

void save(real *u, const char *path)
{
	int i, j;
	FILE *fp;
	
	fp = fopen(path, "w");
	for (j = Y-1; j >= 0; j--) {
	for (i = 0; i <= X-1; i++) {
		fprintf(fp, "%8.3f", *(u+i*Y+j));
		if (i != X-1) {
			fprintf(fp, " ");
		} else {
			fprintf(fp, "\n");
		}
	}
	}
	fclose(fp);
	//printf("%s\n", path);

}

void step(real *u, real *u1)
{
	int i, j;
	
	for (i = (from==0)+from; i < to && i < X-1; i++) {
	for (j = 1; j < Y-1; j++) {
		*(u1+i*Y+j) = *(u+i*Y+j) +
		C * DT / (H * H) * (
		*(u+(i+1)*Y+j) + 
		*(u+(i-1)*Y+j) - 
		*(u+i*Y+j) * 4 +
		*(u+i*Y+j+1) + 
		*(u+i*Y+j-1)
		);
	}
	}
}
