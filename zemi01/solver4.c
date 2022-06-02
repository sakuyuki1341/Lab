#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

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
	int n = 0;
	double dt = 0.0;
	double dx = 0.0;

	int i = 0;
	for (n = 6; n <= 48; n++) {
		for (dt = 0.0001; dt <= 0.05; dt += 0.0001) {
			int isbreak = 0;
			double T[n+1];
			double DT[n+1];
			double T0[n+1];
			double saveT[n+1][11];
			dx = (double)1.0 / (double)n;

			for (i = 0; i <= n/2; i++) {
				T[i] = 2*dx*i;
				T[n-i] = T[i];
			}

			int t = 0;
			for (t = 1; t <= 300; t++) {
				for (i = 0; i <= n; i++) {
					T0[i] = T[i];
				}
				TT(T, DT, dx, dt, n);
				vadd(T, DT, n);
				for (i = 0; i <= n; i++) {
					if (T0[i] < T[i]) {
						isbreak = 1;
					}
				}
			}

			if (isbreak) {
				printf("%d %.4lf\n", n, dt-0.0001);
				break;
			}
		}
	}

	return 0;
}
