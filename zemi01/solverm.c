#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int printout(int n, double* T, int t) {
	// 最初に呼ばれたときのみ実行
	static bool IsFirst = true;
	int i;
	if(IsFirst) {
		// 表最初の行を出力
		printf("time & ");
		for(i = 0; i <= n; i++){
			printf("T%d & ", i);
		}
		printf("\n");
		IsFirst = false;
	}

	// Tの中身を出力
	printf("%d*dt & ", t);
	for(i = 0; i <= n; i++) {
		printf("%.4lf & ", T[i]);
	}
	printf("\n");
}

// ΔTの算出
int calcDT(double* T, double* DT, double dx, double dt, int n) {
	int i = 0;
	DT[0] = 0; DT[n] = 0;
	for (i = 1; i < n; i++) {
		DT[i] = dt*(T[i-1] - 2*T[i] + T[i+1])/(dx*dx);
	}
}

// ΔTからTの算出
int calcT(double* T, double* DT, int n) {
	int i = 0;
	for (i = 1; i < n; i++) {
		T[i] = T[i] + DT[i];
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

	// 初期値設定
	int i = 0;
	for (i = 0; i <= n/2; i++) {
		T[i] = 2*dx*i;
		T[n-i] = T[i];
	}
	printout(n, T, 0);

	// 計算のメインループ. 時間は10*dtまで計算する.
	int t = 0;
	for (t = 1; t <= 10; t++) {
		calcDT(T, DT, dx, dt, n);
		calcT(T, DT, n);
		printout(n, T, t);
	}

	return 0;
}
