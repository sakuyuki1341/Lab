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
    scanf("%d", &n);
 
    double dt = 0.0;
    scanf("%lf", &dt);

    double dx = 0.0;
    dx = (double)1.0 / (double)n;
    
    double T[n+1];
    double T0[n+1];
    double k1[n+1];
    double k2[n+1];
    double k3[n+1];
    double k4[n+1];
	double saveT[n+1][11];
    
    int i = 0;
    for(i = 0; i <= n/2; i++) {
        T[i] = 2*dx*i;
        T[n-i] = T[i];
    }
    printout(n, T, saveT, 0, 1);

    int t = 0;
    for (t = 1; t <= 100; t++) {
        int l = 0;
		for (l = 0; l <= n; l++) {
			T0[l] = T[l];
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
