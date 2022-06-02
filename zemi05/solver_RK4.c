#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "solver_RK4.h"

double mx[1024];
double my[1024];
double mz[1024];
double Hx[1024];
double Hy[1024];
double Hz[1024];
double dt = 0;
double alpha = 0;
double Gamma = 0;
double A = 0;
double ku = 0;
double M = 0;
double lw = 0;
double dx = 0;
double phi = 0;
double interval = 0;
double loops = 0;
int plots = 0;

int main(int argc, char *argv[]) {
	init();

	lw = M_PI * sqrt(A/ku);
	dx = lw/interval;
	phi = M_PI/2;

	// 初期条件設定
	int i = 0;
	for (i = 0; i < interval*2; i++) {
		double x = ((double)i-interval)*dx + dx/2;
		double theta = 2*atan2(exp(M_PI*x/lw), 1);
		mx[i] = sin(theta)*cos(phi);
		my[i] = sin(theta)*sin(phi);
		mz[i] = cos(theta);
		printf("%.8lf %lf\n", x, theta);
//		printf("%lf %lf %lf\n", mx[i], my[i], mz[i]);
	}
	printf("---------------------------------\n");

	int t = 0;
	for (t = 0; t <= loops; t++) {
		if (t == loops) {
			printf("timeout\n");
		}
		RK4();
		if (judge_break()) {
			break;
		}
	}

	for (i = 0; i < interval*2; i++) {
		double x = ((double)i-interval)*dx + dx/2;
		double theta = acos(mz[i]);
		printf("%.8lf %lf\n", x, theta);
//		printf("%lf %lf %lf\n", mx[i], my[i], mz[i]);
	}
	printf("---------------------------------\n");

	double grad = (acos(mz[(int)interval]) - acos(mz[(int)interval-1]))/dx;
	double simlw = M_PI/grad;
	printf("   lw = %.12lf\n", lw);
	printf("simlw = %.12lf\n", simlw);
	return 0;
}

// ---------------------------------------
// 収束判定
// ---------------------------------------
int judge_break() {
	static bool IsFirst = true;
	double avgTmp = 0;
	static double avgFirst = 0;
	static double avgLast = 0;
	static int count = 0;
	int i;
	for (i = 0; i < interval*2; i++) {
		avgTmp += mx[i]*Hx[i] + my[i]*Hy[i] + mz[i]*Hz[i];
	}
	avgTmp = avgTmp/(interval*2);
	count += 1;

	if(IsFirst) {
		avgFirst = avgTmp;
		avgLast = avgTmp;
		IsFirst = false;
		return 0;
	} else {
		if (fabs((avgLast-avgTmp)*1.0e12) <= fabs(avgFirst-avgTmp)) {
			return 1;
		} else {
			avgLast = avgTmp;
			if (count%1000 == 0) {
//				printf("First: %lf  tmp: %lf\n", avgFirst, avgTmp);
			}
			return 0;
		}
	}
}

// ---------------------------------------
//	初期設定
// ---------------------------------------
int init() {
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
		case 'd':
			sscanf(line, "%s %lf", str, &dt);
			break;
		case 'a':
			sscanf(line, "%s %lf", str, &alpha);
			break;
		case 'g':
			sscanf(line, "%s %lf", str, &Gamma);
			break;
		case 'A':
			sscanf(line, "%s %lf", str, &A);
			break;
		case 'k':
			sscanf(line, "%s %lf", str, &ku);
			break;
		case 'M':
			sscanf(line, "%s %lf", str, &M);
			break;
		case 'i':
			sscanf(line, "%s %lf", str, &interval);
			break;
		case 'l':
			sscanf(line, "%s %lf", str, &loops);
			break;
		case 'p':
			sscanf(line, "%s %d", str, &plots);
			break;
		default:
			break;
		}
	}

/*
// デバッグ用
	printf("dt: %lf\n", dt);
	printf("alpha: %lf\n", alpha);
	printf("Gamma: %lf\n", Gamma);
	printf("A    : %lf\n", A);
	printf("ku   : %lf\n", ku);
	printf("M    : %lf\n", M);
	printf("inter: %lf\n", interval);
	printf("loops: %lf\n", loops);
	printf("plots: %d\n", plots);
*/
	return 0;
}


// ---------------------------------------
//	ルンゲクッタ法
// ---------------------------------------
int RK4() {
	double k1x[1024], k1y[1024], k1z[1024];
	double k2x[1024], k2y[1024], k2z[1024];
	double k3x[1024], k3y[1024], k3z[1024];
	double k4x[1024], k4y[1024], k4z[1024];
	double mx0[1024], my0[1024], mz0[1024];

	int i = 0;
	for(i = 0; i < interval*2; i++) {
		mx0[i] = mx[i];
		my0[i] = my[i];
		mz0[i] = mz[i];
	}


	Euler(k1x, k1y, k1z);
	vadd(mx0, my0, mz0, k1x, k1y, k1z, 0.5);
	Euler(k2x, k2y, k2z);
	vadd(mx0, my0, mz0, k2x, k2y, k2z, 0.5);
	Euler(k3x, k3y, k3z);
	vadd(mx0, my0, mz0, k3x, k3y, k3z, 1.0);
	Euler(k4x, k4y, k4z);
	vadd4(mx0, my0, mz0, k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z);
}

int Euler(double* kx, double* ky, double* kz) {
	double Hx_ef[1024], Hy_ef[1024], Hz_ef[1024];
	Heff(Hx_ef, Hy_ef, Hz_ef);
	llg(kx, ky, kz, Hx_ef, Hy_ef, Hz_ef);
	return 0;
}

int llg(double* kx, double* ky, double* kz, double* Hx_ef, double* Hy_ef, double* Hz_ef) {
	int i = 0;
	for (i = 0; i < interval*2; i++) {
		double c = mx[i]*Hx_ef[i] + my[i]*Hy_ef[i] + mz[i]*Hz_ef[i];
		kx[i] = dt * (-fabs(Gamma)) * (my[i]*Hz_ef[i] - mz[i]*Hy_ef[i] + alpha*(c*mx[i] - Hx_ef[i])) / (1 + alpha*alpha);
		ky[i] = dt * (-fabs(Gamma)) * (mz[i]*Hx_ef[i] - mx[i]*Hz_ef[i] + alpha*(c*my[i] - Hy_ef[i])) / (1 + alpha*alpha);
		kz[i] = dt * (-fabs(Gamma)) * (mx[i]*Hy_ef[i] - my[i]*Hx_ef[i] + alpha*(c*mz[i] - Hz_ef[i])) / (1 + alpha*alpha);
	}
	return 0;
}

// ---------------------------------------
//	H に関する計算
// ---------------------------------------
int Heff(double* Hx_ef, double* Hy_ef, double* Hz_ef) {
	Hext(Hx_ef, Hy_ef, Hz_ef);
	HK(Hx_ef, Hy_ef, Hz_ef);
	HA(Hx_ef, Hy_ef, Hz_ef);
	HD(Hx_ef, Hy_ef, Hz_ef);
	// Hに保存
	int i = 0;
	for (i = 0; i < interval*2; i++) {
		Hx[i] = Hx_ef[i];
		Hy[i] = Hy_ef[i];
		Hz[i] = Hz_ef[i];
	}
}

int Hext(double* Hx_ef, double* Hy_ef, double* Hz_ef) {
	int i = 0;
	for (i = 0; i < interval*2; i++) {
		Hx_ef[i] = 0;
		Hy_ef[i] = 0;
		Hz_ef[i] = 0;
	}
}

int HK(double* Hx_ef, double* Hy_ef, double* Hz_ef) {
	int i = 0;
	for (i = 0; i < interval*2; i++) {
		Hz_ef[i] += 2 * ku * mz[i] / M;
	}
}

int HA(double* Hx_ef, double* Hy_ef, double* Hz_ef) {
	int i = 0;
	for (i = 0; i < interval*2; i++) {
		double mxp, mxm, myp, mym, mzp, mzm;
		if (i == 0) {
			mxp = mx[i+1];
			mxm = 0;
			myp = my[i+1];
			mym = 0;
			mzp = mz[i+1];
			mzm = 1;
		} else if(i == interval*2 -1) {
			mxp = 0;
			mxm = mx[i-1];
			myp = 0;
			mym = my[i-1];
			mzp = -1;
			mzm = mz[i-1];
		} else {
			mxp = mx[i+1];
			mxm = mx[i-1];
			myp = my[i+1];
			mym = my[i-1];
			mzp = mz[i+1];
			mzm = mz[i-1];
		}
		Hx_ef[i] += 2*A*(mxp - 2*mx[i] + mxm)/(M*dx*dx);
		Hy_ef[i] += 2*A*(myp - 2*my[i] + mym)/(M*dx*dx);
		Hz_ef[i] += 2*A*(mzp - 2*mz[i] + mzm)/(M*dx*dx);
	}
}

int HD(double* Hx_ef, double* Hy_ef, double* Hz_ef) {
	int i = 0;
	for (i = 0; i < interval*2; i++) {
		Hx_ef[i] += - 4*M_PI*mx[i];
	}
}


// ---------------------------------------
//	vadd 
// ---------------------------------------
int vadd(double* mx0, double* my0, double* mz0, double* kx, double* ky, double* kz, double r) {
	int i = 0;
	for (i = 0; i < interval*2; i++) {
		mx[i] = mx0[i] + kx[i]*r;
		my[i] = my0[i] + ky[i]*r;
		mz[i] = mz0[i] + kz[i]*r;

		double absm = sqrt(mx[i]*mx[i] + my[i]*my[i] + mz[i]*mz[i]);
		mx[i] = mx[i]/absm;
		my[i] = my[i]/absm;
		mz[i] = mz[i]/absm;
	}
	return 0;
}

int vadd4(double* mx0, double* my0, double* mz0, double* k1x, double* k1y, double* k1z, double* k2x, double* k2y, double* k2z, double* k3x, double* k3y, double* k3z, double* k4x, double* k4y, double* k4z) {
	int i = 0;
	for (i = 0; i < interval*2; i++) {
		mx[i] = mx0[i] + (k1x[i] + 2*k2x[i] + 2*k3x[i] + k4x[i])/6;
		my[i] = my0[i] + (k1y[i] + 2*k2y[i] + 2*k3y[i] + k4y[i])/6;
		mz[i] = mz0[i] + (k1z[i] + 2*k2z[i] + 2*k3z[i] + k4z[i])/6;

		double absm = sqrt(mx[i]*mx[i] + my[i]*my[i] + mz[i]*mz[i]);
		mx[i] = mx[i]/absm;
		my[i] = my[i]/absm;
		mz[i] = mz[i]/absm;
	}
	return 0;
}


// ---------------------------------------
//	出力用関数
// ---------------------------------------
