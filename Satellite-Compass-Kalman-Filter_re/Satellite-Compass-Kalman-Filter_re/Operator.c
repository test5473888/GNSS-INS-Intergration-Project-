/*
matrix operate : mul,inverse(gaussj),Transpose
generate ramdom varible
and two test function
*/
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
//#include<fstream>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

//-------------------------
//matrix multiplication
void  MatrixMul(double **m1, int i1, int j1, double **m2, int i2, int j2, double **p)
{
	int i, j, k;
	for (i = 0; i<i1; i++)
	{
		for (j = 0; j<j2; j++)
		{
			p[i][j] = 0.0;
			for (k = 0; k<j1; k++)
			{
				p[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
};

void  MatrixMul_S(double **m1, int i1, int j1, double *m2, int j2, double *p)
{
	/*if(i2!=j1)
	{
	return NULL;
	} */
	int i, j;
	for (i = 0; i<i1; i++)
	{
		p[i] = 0.0;
		for (j = 0; j<j2; j++)
		{
			p[i] += m1[i][j] * m2[j];
		}
	}
};

void MatPlus(double **m1,double **m2,int a, int b) {
	for (int i = 0; i < a; i++)
		for (int j = 0; j < b; j++)
			m1[i][j] = m1[i][j] + m2[i][j];
};

void MatMinus(double **m1, double **m2, int a, int b) {
	for (int i = 0; i < a; i++)
		for (int j = 0; j < b; j++)
			m1[i][j] = m1[i][j] - m2[i][j];
};


void  SWAP(double *a, double *b)
{
	double *temp;
	temp = a;
	a = b;
	b = temp;
};


//-------------------------
//inverse matrix
void  gaussj(double **a, int n)
{
	int *indxc, *indxr, *ipiv;
	int i, icol, irow, j, k, l, ll;
	double dum, big, pivinv;	 
	
	indxc = (int*)malloc(n * sizeof(int));		//indxc = new int[n];
	indxr = (int*)malloc(n * sizeof(int));		//indxr = new int[n];
	ipiv = (int*)malloc(n * sizeof(int));		//ipiv = new int[n]
	
	irow = 0;
	icol = 0;

	for (j = 0; j<n; j++)
		ipiv[j] = 0;

	for (i = 0; i<n; i++) {
		big = 0.0;
		for (j = 0; j<n; j++)
			if (ipiv[j] != 1)
				for (k = 0; k<n; k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big = (float)fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l = 0; l<n; l++)
				SWAP(&a[irow][l], &a[icol][l]);
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0)
			printf("gaussj: Singular Matrix");
		pivinv = (float)1.0 / a[icol][icol];
		a[icol][icol] = 1.0;
		for (l = 0; l<n; l++)
			a[icol][l] *= pivinv;
		for (ll = 0; ll<n; ll++)
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l = 0; l<n; l++)
					a[ll][l] -= a[icol][l] * dum;
			}
	}

	for (l = n - 1; l >= 0; l--) {
		if (indxr[l] != indxc[l])
			for (k = 0; k<n; k++)
				SWAP(&a[k][indxr[l]], &a[k][indxc[l]]);
	}

	free(indxc);		//delete[]indxc;
	free(indxr);		//delete[]indxr;
	free(ipiv);			//delete[]ipiv;
};

//-------------------------
//Transpose
void mat_trans_N(int **A, int** rst, int m, int n) {

	int tempt;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
		{
			tempt = A[j][i];
			rst[i][j] = tempt;
		}
}

void mat_trans(double **A, double** rst, int m, int n) {

	double tempt;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
		{
			tempt = A[j][i];
			rst[i][j] = tempt;
		}
}

double** mat_eye(int e) {
	double **temp = (double**)malloc(e * sizeof(double*) );
	for (int i=0; i < e; i++) {
		temp[i] = (double*)calloc(e, sizeof(double));
		for (int j=0; j < e; j++)
			if(i==j) temp[i][j] = 1;
	}
	return temp;
		
}

//-------------------------
//generate r.v.

double generaterv(double mean, double std)
{
	double u = rand() / (double)RAND_MAX;
	double v = rand() / (double)RAND_MAX;
	double x = sqrt(-2 * log(u)) * cos(2 * M_PI * v) * std + mean;
	return x;
}

//-------------------------
//Textbook Example (2)
void test_matopr() {

	void mat_trans(int **A, int** rst, int m, int n);	
	void mat_trans(double **A, double** rst, int m, int n);
	void  MatrixMul(double **m1, int i1, int j1, double **m2, int i2, int j2, double **p);
	void gaussj(double **a, int n);

	double **A, **N, **At, **N_inv;
	double **L, **x, *l;
	int sq_msize = 3, sq_nsize = 2;

	A = (double**)malloc(sq_msize * sizeof(double*));
	N = (double**)malloc(sq_nsize * sizeof(double*));
	At = (double**)malloc(sq_nsize * sizeof(double*));
	N_inv = (double**)malloc(sq_nsize * sizeof(double*));
	L = (double**)malloc(sq_msize * sizeof(double*));
	l = (double*)malloc(sq_msize * sizeof(double*));
	x = (double**)malloc(sq_nsize * sizeof(double*));	

	double m1[2] = { 1, 1 };	double l1[1] = { 393.65 };
	double m2[2] = { 1, 0 };	double l2[1] = { 190.4 };
	double m3[2] = { 0, 1 };	double l3[1] = { 203.16 };

	A[0] = m1, A[1] = m2, A[2] = m3;
	L[0] = l1, L[1] = l2, L[2] = l3;

	At[0] = (double*)malloc(sq_msize * sizeof(double)), At[1] = (double*)malloc(sq_msize * sizeof(double));
	N[0] = (double*)malloc(sq_nsize * sizeof(double)), N[1] = (double*)malloc(sq_nsize * sizeof(double));
	N_inv[0] = (double*)malloc( sizeof(double)), N_inv[1] = (double*)malloc(sizeof(double));
	x[0] = (double*)malloc(sq_nsize * sizeof(double)), x[1] = (double*)malloc(sq_nsize * sizeof(double));
	

	//Âà¸m¯x°}
	mat_trans(A, At, 3, 2);
	MatrixMul(At, 2, 3, A,3, 2,N);

	//°f¯x°}
	gaussj(N, 2);
	MatrixMul(At,2, 3,  L, 3, 1,x);
	MatrixMul(N, 2, 2,  x, 2, 1, N_inv);
	printf("%lf\n%lf", N_inv[0][0], N_inv[1][0]);
	
		free(At[0]), free(At[1]);
		free(N[0]), free(N[1]);
		free(N_inv[0]), free(N_inv[1]);
		free(x[0]), free(x[1]);		
		free(x), free(L), free(N_inv), free(N), free(At), free(A);
}

void test_rv() {
	double generaterv(double mean, double std);
	FILE *infile;
	fopen_s(&infile, ".\\data\\rbdata.csv", "a");	

		for (int j = 0; j < 9999; j++) 			
			fprintf_s(infile, "%lf,%lf,%lf\n",generaterv(0, 1), generaterv(0, 1), generaterv(0, 1));	
		
}