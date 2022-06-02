#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "solver_RK4.h"



int main(int argc, char *argv[]) {
	double gamma = 0;
	double alpha = 0.0;
	double dt = 0.0;
	double ku = 0;
	double m[3];
	double H[3];
	double loops = 0;
	int plots = 0;
	init(m, H, &dt, &alpha, &gamma, &ku, &loops, &plots);

	H[0] = sin(0.01)*1.5e4; H[1] = 0; H[2] = -cos(0.01)*1.5e4;

	int t = 0;
	for (t = 0; t <= loops; t++) {
		RK4(m, H, dt, alpha, gamma, ku);
		if (t % plots == 0) {
			if (strcmp(argv[1], "3") == 0) {
				printout3(m, t);
			} else if (strcmp(argv[1], "2") == 0) {
				printout2(m, t, dt);
			}
		}
	}
	return 0;
}

// ---------------------------------------
//	初期設定
// ---------------------------------------
int init(double* m, double* H, double* dt, double* alpha, double* gamma, double* ku, double* loops, int* plots) {
	FILE *fp;
	const int n = 256;
	char fname[] = "init.data";
	char line[n];
	char str[16];
	
	fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("init.data not open!\n");
		exit(1);
	}

	// 中身の取得
	while(fgets(line, n, fp) != NULL) {
		switch (line[0]){
		case '#':
			break;
		case 'm':
			sscanf(line, "%s %lf %lf %lf", str, &m[0], &m[1], &m[2]);
			break;
		case 'H':
			sscanf(line, "%s %lf %lf %lf", str, &H[0], &H[1], &H[2]);
			break;
		case 'd':
			sscanf(line, "%s %lf", str, dt);
			break;
		case 'a':
			sscanf(line, "%s %lf", str, alpha);
			break;
		case 'g':
			sscanf(line, "%s %lf", str, gamma);
			break;
		case 'l':
			sscanf(line, "%s %lf", str, loops);
			break;
		case 'k':
			sscanf(line, "%s %lf", str, ku);
			break;
		case 'p':
			sscanf(line, "%s %d", str, plots);
			break;
		default:
			break;
		}
	}

	/* デバッグ用
	printf("mx: %lf, my: %lf, mz: %lf\n", m[0], m[1], m[2]);
	printf("Hx: %lf, Hy: %lf, Hz: %lf\n", H[0], H[1], H[2]);
	printf("dt: %lf\n", *dt);
	printf("alpha: %lf\n", *alpha);
	printf("gamma: %lf\n", *gamma);
	printf("loops: %lf\n", *loops);
	*/
	return 0;
}


// ---------------------------------------
//	ルンゲクッタ法
// ---------------------------------------
int RK4(double* m, double* H, double dt, double alpha, double gamma, double ku) {
	double k1[3], k2[3], k3[3], k4[4];
	double m0[3];
	int i = 0;
	for (i = 0; i < 3; i++) {
		m0[i] = m[i];
	}
	Euler(m, H, k1, dt, alpha, gamma, ku);
	vadd(m, m0, k1, 0.5);
	Euler(m, H, k2, dt, alpha, gamma, ku);
	vadd(m, m0, k2, 0.5);
	Euler(m, H, k3, dt, alpha, gamma, ku);
	vadd(m, m0, k3, 1.0);
	Euler(m, H, k4, dt, alpha, gamma, ku);
	vadd4(m, m0, k1, k2, k3, k4);
}

int Euler(double* m, double* H, double* k, double dt, double alpha, double gamma, double ku) {
	double mx = m[0], my = m[1], mz = m[2];
	double Hx, Hy, Hz;
	Heff(m, H, &Hx, &Hy, &Hz, ku);
	llg(k, mx, my, mz, Hx, Hy, Hz, dt, alpha, gamma);
	return 0;
}

int llg(double* k, double mx, double my, double mz, double Hx, double Hy, double Hz, double dt, double alpha, double gamma) {
	double c = mx*Hx + my*Hy + mz*Hz;
	k[0] = dt * (-fabs(gamma)) * (my*Hz - mz*Hy + alpha*(c*mx - Hx)) / (1 + alpha*alpha);
	k[1] = dt * (-fabs(gamma)) * (mz*Hx - mx*Hz + alpha*(c*my - Hy)) / (1 + alpha*alpha);
	k[2] = dt * (-fabs(gamma)) * (mx*Hy - my*Hx + alpha*(c*mz - Hz)) / (1 + alpha*alpha);
	return 0;
}

// ---------------------------------------
//	H に関する計算
// ---------------------------------------
int Heff(double* m, double* H, double* Hx, double* Hy, double* Hz, double ku) {
	Hext(H, Hx, Hy, Hz);
	Hanis(m, H, Hx, Hy, Hz, ku);
}

int Hext(double* H, double* Hx, double* Hy, double* Hz) {
	*Hx = H[0];
	*Hy = H[1];
	*Hz = H[2];
}

int Hanis(double* m, double* H, double* Hx, double* Hy, double* Hz, double ku) {
	*Hz += 2*ku*m[2] / sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
}

// ---------------------------------------
//	vadd 
// ---------------------------------------
int vadd(double* m, double* m0, double* k, double r) {
	int i = 0;
	for(i = 0; i < 3; i++) {
		m[i] = m0[i] + k[i]*r;
	}

	double M = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
	for(i = 0; i < 3; i++) {
		m[i] = m[i]/M;
	}
	return 0;
}

int vadd4(double* m, double* m0, double* k1, double* k2, double* k3, double* k4) {
	int i = 0;
	for (i = 0; i < 3; i++) {
		m[i] = m0[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6;
	}

	double mabs = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
	for(i = 0; i < 3; i++) {
		m[i] = m[i]/mabs;
	}
	return 0;
}


// ---------------------------------------
//	出力用関数
// ---------------------------------------
int printout2(double* m, int t, double dt) {
	int i = 0;
	double ms = (double)t*dt*1.0e6;
	printf("%lf %.6lf", ms, m[2]);
	printf("\n");
	return 0;
}

int printout3(double* m, int t) {
	int i = 0;
	for (i = 0; i < 3; i++) {
		printf("%.6lf ", m[i]);
	}
	printf("\n");
	return 0;
}