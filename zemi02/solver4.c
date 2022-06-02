#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int TT(double* T, double* k, double dx, double dt, int n) {
    int i = 0;
    k[0] = 0; k[n] = 0;
    for (i = 1; i < n; i++) {
        k[i] = dt*(T[i-1] - 2*T[i] + T[i+1])/(dx*dx);
    }
}

int vadd(double* T, double* T0, double* k, int n, double r) {
    int i = 0;
    for(i = 1; i < n; i++) {
        T[i] = T0[i] + k[i]*r;
    }
}

int vadd4(double* T, double* T0, double* k1, double* k2, double* k3, double* k4, int n) {
    int i = 0;
    for (i = 1; i < n; i++) {
        T[i] = T0[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6;
    }
}

int main() {
    int n = 0; 
    double dt = 0.0;
    double dx = 0.0;
    

    int i = 0;
	int j = 0;
	for (n = 6; n <= 48; n++) {
		for (dt = 0.0001; dt <= 0.05; dt+=0.0001) {
			double isbreak = 0;
			double T[n+1];
			double T0[n+1];
			double k1[n+1];
			double k2[n+1];
			double k3[n+1];
			double k4[n+1];
			dx = (double)1.0 / (double)n;

			for(i = 0; i <= n/2; i++) {
				T[i] = 2*dx*i;
				T[n-i] = T[i];
			}

			int t = 0;
			for (t = 1; t <= 300; t++) {
				int l = 0;
				for (l = 0; l <= n; l++) {
					T0[l] = T[l];
					if (T0[l] < 0) {
						isbreak = 1;
					}
				}
				TT(T, k1, dx, dt, n);
				vadd(T, T0, k1, n, 0.5);
				TT(T, k2, dx, dt, n);
				vadd(T, T0, k2, n, 0.5);
				TT(T, k3, dx, dt, n);
				vadd(T, T0, k3, n, 1.0);
				TT(T, k4, dx, dt, n);
				vadd4(T, T0, k1, k2, k3, k4, n);
			}

			if (isbreak) {
				printf("%d %.4lf\n", n, dt-0.0001);
				break;
			}
		}
	}

	return 0;
}
