#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include"Kalman.h"
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif


double generaterv(double mean, double std);
//matrix math
void  MatrixMul(double **m1, int i1, int j1, double **m2, int i2, int j2, double **p);
void  MatrixMul_S(double **m1, int i1, int j1, double *m2, int j2, double *p);
void mat_trans(double **A, double** rst, int m, int n);
void mat_trans_N(int **A, int** rst, int m, int n);
double** mat_eye(int e);
void MatPlus(double **m1, double **m2, int a, int b);
void MatMinus(double **m1, double **m2, int a, int b);
//string work
size_t getline(char* line, size_t n, FILE *fp);



int main()
{
	FILE *fgyro, *frtk,*fdeg;
	struct Kalman filter;
	struct Kalman ini_KF(double *x, double**p);
	void KF_set_process(struct Kalman *KF, double** q, double**f, double**g, double**u);
	void KF_set_measure(struct Kalman *KF, double **h, double *z, double *r);
	void  KF_pridict(struct Kalman* kf);
	void  KF_update(struct Kalman* kf);
	char **line;


	fgyro = fopen(".\\data\\2019-3-6_acc_gyro-s.txt", "r");
	frtk = fopen(".\\data\\Degree_065-s.txt", "r");
	fdeg = fopen(".\\data\\Out_065.txt", "a+");
	double degree, z0, *X, *z,*u;
	double **P, **F, **Q, **H, **R,**G;
	char buff[10];
	char *rtk_line, *gyro_line, *gy_temp;
	int i,j = 0, n = 2;	
	size_t sn = 67;
	double rtk_deg, f_temp, gz;	
	//initialize memory
	X = calloc(2, sizeof(double));
	z = calloc(2, sizeof(double));
	u = calloc(2, sizeof(double));
	P = malloc(2 * sizeof(double*));
	F = malloc(2 * sizeof(double*));
	Q = malloc(2 * sizeof(double*));
	G = malloc(2 * sizeof(double*));
	filter.z = calloc(2, sizeof(double));
	for (i = 0; i < n; i++) P[i] = calloc(n, sizeof(double)), P[i][i] = 0.01;
	for (int i = 0; i<n; i++) Q[i] = calloc(n, sizeof(double)), G[i] = calloc(n, sizeof(double));
	rtk_line = malloc(sn * sizeof(char));
	gyro_line = malloc(sn * sizeof(char));
	gy_temp = malloc(sn * sizeof(char));
	
	F = mat_eye(n);
	H = mat_eye(n);
	R = mat_eye(n);
	G = mat_eye(n);

	F[1][1] = 0;
	Q[0][0] = 1.1*1.1;
	Q[1][1] = 1.1*1.1;	
	R[0][0] = 0.02*0.02;
	R[1][1] = 0.02*0.02;
	G[0][0] = 0;
	G[0][1] = 1;


	getline(rtk_line, sn, frtk);						//get first data
	getline(gyro_line, sn, fgyro);
	getline(gyro_line, sn, fgyro);
	sscanf(rtk_line, "%lf,%lf", &rtk_deg, &f_temp);
	sscanf(gyro_line, "%s %s %s %s %s %s %s %lf %lf %lf",
		gy_temp, gy_temp, gy_temp, gy_temp, gy_temp, gy_temp,
		gy_temp, &f_temp, &f_temp, &gz);

	z0 = rtk_deg;
	X[0] = rtk_deg;	X[1] = 0;
	u[0] = 0;	u[1] = gz;

	filter = ini_KF(X, P);
	KF_set_process(&filter, Q, F,G,u);


	//getline	
	while (j != -1) {
		if (j = getline(rtk_line, sn, frtk))
			sscanf(rtk_line, "%lf,%lf", &rtk_deg, &f_temp);
		if (j = getline(gyro_line, sn, fgyro))
			sscanf(gyro_line, "%s %s %s %s %s %s %s %lf %lf %lf",
				gy_temp, gy_temp, gy_temp, gy_temp, gy_temp, gy_temp,
				gy_temp, &f_temp, &f_temp, &gz);

		z[0] = rtk_deg;
		z[1] = rtk_deg - z0;
		z0 = rtk_deg;
		u[1] = gz;
		KF_set_process(&filter, Q, F, G, u);
		//if (z[1] > 30 || z[1] < -30) filter.Q[0][0] = 1000;
		KF_pridict(&filter);
		KF_set_measure(&filter, H, z, R);
		KF_update(&filter);
		//printf("%.4lf\n", filter.X[0]);
		fprintf(fdeg,"%.4lf\n", filter.X[0]);
		//fprintf(fdeg, "%.4lf,%.4lf\n", rtk_deg, gz * 180 / M_PI);
	}


	fclose(fgyro);
	fclose(frtk);
	fclose(fdeg);

	//release memory
	for (i = 0; i < n; i++) free(P[i]), free(Q[i]);
	free(X);	free(z);
	free(P);	free(F);
	free(Q);	free(filter.z);

	free(rtk_line);
	return 0;
}




struct Kalman ini_KF(double *x, double**p) {
	struct Kalman temp;
	temp.X = x;
	temp.P = p;
	return temp;
}

void KF_set_process(struct Kalman *KF, double** q, double**f, double**g, double**u) {
	KF->Q = q;		KF->F = f;		KF->G = g;		KF->u = u;
}

void  KF_pridict(struct Kalman* kf) {
	double **xt, **FX, **x, **FT, **FP,**Gu,**u,**FQ,**FQFT;
	int n = 2;

	x = (double**)malloc(sizeof(double*));
	xt = (double**)malloc(n * sizeof(double*));
	u = (double**)malloc(n * sizeof(double*));
	FX = (double**)malloc(n * sizeof(double*));
	FT = (double**)malloc(n * sizeof(double*));
	FP = (double**)malloc(n * sizeof(double*));
	FQ = (double**)malloc(n * sizeof(double*));
	FQFT = (double**)malloc(n * sizeof(double*));
	Gu = (double**)malloc(n * sizeof(double*));
	x[0] = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++)
		xt[i] = (double*)malloc(sizeof(double)),
		u[i] = (double*)malloc(sizeof(double)),
		FX[i] = (double*)malloc(sizeof(double)),
		FT[i] = (double*)malloc(n * sizeof(double)),
		Gu[i] = (double*)malloc(n * sizeof(double)),
		FQ[i] = (double*)malloc(n * sizeof(double)),
		FQFT[i] = (double*)malloc(n * sizeof(double)),
		FP[i] = (double*)malloc(n * sizeof(double));

	for (int i = 0; i < n; i++) x[0][i] = kf->X[i], u[i][0] = kf->u[i];

	mat_trans(x, xt, 1, n);
	mat_trans(kf->F, FT, n, n);
	MatrixMul(kf->F, n, n, xt, n, 1, FX);					//FX
	MatrixMul(kf->G, n, n, u, n, 1, Gu);					//gu
	MatPlus(FX, Gu, n, 1);									//x = FX+Gu

	MatrixMul(kf->F, n, n, kf->Q, n, n, FQ);					//FQ
	MatrixMul(FQ, n, n, FT, n, n, FQFT);						//FQ*FT
	for (int i = 0; i < n; i++)  kf->X[i] = FX[i][0];
	MatrixMul(kf->F, n, n, kf->P, n, n, FP);
	MatrixMul(FP, n, n, FT, n, n, kf->P);			// = FP*FT + FQ*FT
	MatPlus(kf->P, FQFT, n, n);
	//for (int i = 0; i < 3; i++)  kf->F[i][i] = kf->F[i][i] + 1;		//FPFT + Q	
	for (int i = 0; i < n; i++)
		free(xt[i]), free(FX[i]), free(FT[i]), free(FP[i]);
	free(x[0]);
	free(x), free(xt), free(FX), free(FT), free(FP);

}

void KF_set_measure(struct Kalman *KF, double **h, double *z, double *r) {
	KF->H = h;	KF->z = z;	KF->R = r;
}

void  KF_update(struct Kalman* kf) {
	double **K, **PHT, **x, **HT, **S, **HX, **Kv;
	double **KH, **Pu;
	double **v;
	int n = 2;
	x = (double**)malloc(n * sizeof(double*));
	K = (double**)malloc(n * sizeof(double*));
	PHT = (double**)malloc(n * sizeof(double*));
	HT = (double**)malloc(n * sizeof(double*));
	S = (double**)malloc(sizeof(double*));
	HX = (double**)malloc(sizeof(double*)),
		Kv = (double**)malloc(n * sizeof(double*));
	KH = (double**)malloc(n * sizeof(double*));
	Pu = (double**)malloc(n * sizeof(double*));
	v = (double**)malloc(n * sizeof(double*));

	for (int i = 0; i < n; i++)
		x[i] = (double*)malloc(sizeof(double)), x[i][0] = kf->X[i],
		K[i] = (double*)malloc(n * sizeof(double)),
		PHT[i] = (double*)malloc(n * sizeof(double)),
		HT[i] = (double*)malloc(n * sizeof(double)),
		Kv[i] = (double*)malloc(n * sizeof(double)),
		KH[i] = (double*)malloc(n * sizeof(double)),
		Pu[i] = (double*)malloc(n * sizeof(double)),
		S[i] = (double*)malloc(n * sizeof(double)),
		HX[i] = (double*)malloc(sizeof(double)),
		v[i] = (double*)calloc(1, sizeof(double)), v[i][0] = kf->z[i];


	mat_trans(kf->H, HT, n, n);
	MatrixMul(kf->P, n, n, HT, n, n, PHT);
	MatrixMul(kf->H, n, n, PHT, n, n, S);
	MatPlus(S, kf->R, n, n);	// S = S + R
	gaussj(S, n);
	MatrixMul(PHT, n, n, S, n, n, K);	// get K

	MatrixMul(kf->H, n, n, x, n, 1, HX);// X updata	
	MatMinus(v, HX, n, 1);
	MatrixMul(K, n, n, v, n, 1, Kv);
	MatPlus(x, Kv, n, 1);

	for (int i = 0; i < n; i++)	kf->X[i] = x[i][0];

	MatrixMul(K, n, n, kf->H, n, n, KH);		//// P updata
	double ** I = mat_eye(n);
	MatMinus(I, KH, n, n);
	MatrixMul(I, n, n, kf->P, n, n, Pu);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			kf->P[i][j] = Pu[i][j];


	// release memory
	/*
	for (int i = 0; i < n; i++)
	x[i] = (double*)malloc(sizeof(double)), x[i][0] = kf->X[i],
	K[i] = (double*)malloc(n * sizeof(double)),
	PHT[i] = (double*)malloc(n * sizeof(double)),
	HT[i] = (double*)malloc(n * sizeof(double)),
	Kv[i] = (double*)malloc(n * sizeof(double)),
	KH[i] = (double*)malloc(n * sizeof(double)),
	Pu[i] = (double*)malloc(n * sizeof(double)),
	S[i] = (double*)malloc(n * sizeof(double)),
	HX[i] = (double*)malloc(sizeof(double)),
	v[i] = (double*)calloc(1, sizeof(double)), v[i][0] = kf->z[i];
	*/
}


size_t getline(char* line, size_t n, FILE *fp) {
	char *buf = line;
	char c;
	size_t i = 0;
	while ((c = fgetc(fp)) != '\n')
	{
		if (c == EOF)
			return -1;
		if (i < n - 2)
		{
			buf[i] = c;
			i++;
		}
		else
		{
			n = n + 10;
			buf = realloc(buf, n);
			buf[i] = c;
			i++;
		}
	}
	buf[i] = '\n';
	buf[i + 1] = '\0';
	i = i + 2;
	return i;
}