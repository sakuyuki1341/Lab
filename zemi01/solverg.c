#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int printout(int n, double* T, double saveT[][11], int t, int issave) {
	int i = 0, j = 0;
	double dx = (double)1.0 / (double)n;
	if (issave) {
		for (i = 0; i <= n; i++) {
			saveT[i][t] = T[i];
		}
	} else {
		for (i = 0; i <= n; i++) {
			printf("%.4lf ", dx*i);
			for (j = 0; j <= 10; j++) {
				printf("%.4lf ", saveT[i][j]);
			}
			printf("\n");
		}
	}
}

// ΔTの算出
int TT(double* T, double* k, double dx, double dt, int n) {
	int i = 0;
	k[0] = 0; k[n] = 0;
	for (i = 1; i < n; i++) {
		k[i] = dt*(T[i-1] - 2*T[i] + T[i+1])/(dx*dx);
	}
}

// ΔTからTの算出
int vadd(double* T, double* k, int n) {
	int i = 0;
	for (i = 1; i < n; i++) {
		T[i] = T[i] + k[i];
	}
}

int main() {
	// nの取得
	int n = 0;
	scanf("%d", &n);
	// dtの取得
	double dt = 0.0;
	scanf("%lf", &dt);
	// dxの計算
	double dx;
	dx = (double)1.0 / (double)n;
	// T, ΔTの用意
	double T[n+1];
	double DT[n+1];
	double saveT[n+1][11];

	// 初期値設定
	int i = 0;
	for (i = 0; i <= n/2; i++) {
		T[i] = 2*dx*i;
		T[n-i] = T[i];
	}
	printout(n, T, saveT, 0, 1);

	// 計算のメインループ. 時間は10*dtまで計算する.
	int t = 0;
	for (t = 1; t <= 10; t++) {
		TT(T, DT, dx, dt, n);
		vadd(T, DT, n);
		printout(n, T, saveT, t, 1);
	}

	printout(n, T, saveT, 0, 0);
	return 0;
}
