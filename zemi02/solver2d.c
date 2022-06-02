#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int printout(int n, double T[][11], double saveT[][11][11], int t, int issave) {
	int i = 0, j = 0, k = 0;
	double dx = (double)1.0 / (double)n;
	if (issave) {
		for (i = 0; i <= n; i++) {
			for (j = 0; j <= n; j++) {
				saveT[i][j][t] = T[i][j];
			}
		}
	} else {
		for (i = 0; i <= n; i++) {
			for (j = 0; j <= n; j++) {
				//printf("test");
				printf("%.4lf %.4lf ", dx*i, dx*j);
				for (k = 0; k <= 10; k++) {
					//printf("a");
					printf("%.4lf ", saveT[i][j][k]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}
}

int TT(double T[][11], double k[][11], double dx, double dt, int n) {
	int i = 0, j = 0;
	for (i = 0; i <= n; i++) {
		k[0][i] = 0; k[n][i] = 0; k[i][0] = 0; k[i][n] = 0; 
	}
	for (i = 1; i < n; i++) {
		for (j = 1; j < n; j++) {
			k[i][j] = dt*(T[i-1][j] - 2*T[i][j] + T[i+1][j] + T[i][j-1] - 2*T[i][j] + T[i][j+1])/(dx*dx);
		}
	}
}

int vadd(double T[][11], double T0[][11], double k[][11], int n, double r) {
    int i = 0, j = 0;
    for(i = 1; i < n; i++) {
		for (j = 1; j < n; j++) {
			T[i][j] = T0[i][j] + k[i][j]*r;
		}
    }
}

int vadd4(double T[][11], double T0[][11], double k1[][11], double k2[][11], double k3[][11], double k4[][11], int n) {
    int i = 0, j = 0;
    for (i = 1; i < n; i++) {
		for (j = 1; j < n; j++) {
			T[i][j] = T0[i][j] + (k1[i][j] + 2*k2[i][j] + 2*k3[i][j] + k4[i][j])/6;
		}
    }
}

int main() {
	int n = 10;

	double dt = 0.001;

	double dx = 0.0;
	dx = 1.0 / (double)n;

	double T[n+1][n+1];
	double T0[n+1][n+1];
	double k1[n+1][n+1];
	double k2[n+1][n+1];
	double k3[n+1][n+1];
	double k4[n+1][n+1];
	double saveT[n+1][n+1][11];

	int i = 0, j = 0;
	for(i = 0; i <= n/2; i++) {
		for (j = 0; j <= n/2; j++) {
			T[i][j] = 4*dx*dx*i*j;
			T[n-i][n-j] = T[i][j];
			T[i][n-j] = T[i][j];
			T[n-i][j] = T[i][j];
		}
	}
    printout(n, T, saveT, 0, 1);

    int t = 0;
    for (t = 1; t <= 100; t++) {
		for (i = 0; i <= n; i++) {
			for (j = 0; j <= n; j++) {
				T0[i][j] = T[i][j];
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
		if (t % 10 == 0) {
			printout(n, T, saveT, t/10, 1);
		}
	}

	printout(n, T, saveT, 0, 0);
	return 0;
}
