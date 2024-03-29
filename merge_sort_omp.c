/*
 * Николай Игоревич Хохлов, k_h@inbox.ru, 2014.
 * Реализация алгоритма сортировки слиянием.
 * Распараллелить используя OpenMP и алгоритм гиперкуб.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

void merge(int *a, int *b, int *c, int na, int nb);
void merge_sort(int *a, int na);
void print(int *a, int na);
int check_sort(int *a, int n);
double timer();
void dec(const int p, const int k, const int n, int *sk, int *fk);

int pow2(int i)
{
	return 1 << i;
}

int main (int argc, char *argv[])
{
	int i, n, *a;
	double t;
	if (argc < 2) {
		printf("Usage: %s num_elements.\n", argv[0]);
		return 1;
	}
	n = atoi(argv[1]);
	a = (int*)malloc(sizeof(int) * n);
	srand(time(NULL));

	for (i = 0; i < n; i++) a[i] = rand() % 1000;

	if (n < 101) print(a, n);

	t = timer();

	int p;
	#pragma omp parallel
	{
		p = omp_get_num_threads();
		int k = omp_get_thread_num();
		int sk, fk;
		dec(p, k, n, &sk, &fk);
		printf("#%d: sk = %d, fk = %d\n", k, sk, fk);
		merge_sort(a + sk, fk - sk);
		#pragma omp barrier
		// if (k == 0)
		#pragma omp master
		{
			int i;
			for (i = 1; i < p; i++) {
				int si, fi;
				dec(p, i, n, &si, &fi);
				int *c = (int*)malloc(sizeof(int) * (fi - si + fk - sk));
				merge(a + sk, a + si, c, fk - sk, fi - si);
				memcpy(a + sk, c, sizeof(int) * (fi - si + fk - sk));
				fk = fi;
				free(c);
			}
		}
	}

	t = timer() - t;

	if (n < 101) print(a, n);
	printf("Time: %f sec, sorted: %d\n", t, check_sort(a, n));
	free(a);
	return 0;
}

void dec(const int p, const int k, const int n, int *sk, int *fk)
{
	*sk = k * (n / p);
	*fk = *sk + n / p;
	if (k == p - 1) *fk = n;
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
	int i = 0, j = 0;
	while (i < na && j < nb) {
		if (a[i] <= b[j]) {
			c[i + j] = a[i];
			i++;
		} else {
			c[i + j] = b[j];
			j++;
		}
	}
	if (i < na) {
		memcpy(c + i + j, a + i, (na - i) * sizeof(int));
	} else {
		memcpy(c + i + j, b + j, (nb - j) * sizeof(int));
	}
}

/*
 * Процедура сортировки слиянием.
 */
void merge_sort(int *a, int na)
{
	if(na < 2) return;
	merge_sort(a, na / 2);
	merge_sort(a + na / 2, na - na / 2);

	int *b = (int*)malloc(sizeof(int) * na);
	
	merge(a, a + na / 2, b, na / 2, na - na / 2);
	
	memcpy(a, b, sizeof(int) * na);
	
	free(b);
}
