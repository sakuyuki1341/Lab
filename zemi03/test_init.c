#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

int init(double* m, double* H, double* dt, double* alpha, double* gamma) {
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
		default:
			break;
		}
	}

	printf("mx: %lf, my: %lf, mz: %lf\n", m[0], m[1], m[2]);
	printf("Hx: %lf, Hy: %lf, Hz: %lf\n", H[0], H[1], H[2]);
	printf("dt: %lf\n", *dt);
	printf("alpha: %lf\n", *alpha);
	printf("gamma: %lf\n", *gamma);
	return 0;
}

int main() {
	double gamma = 0;
	double alpha = 0.0;
	double dt = 0.0;
	double m[3];
	double H[3];
	init(m, H, &dt, &alpha, &gamma);
	
	printf("---------------------------------\n");
	printf("mx: %lf, my: %lf, mz: %lf\n", m[0], m[1], m[2]);
	printf("Hx: %lf, Hy: %lf, Hz: %lf\n", H[0], H[1], H[2]);
	printf("dt: %lf\n", dt);
	printf("alpha: %lf\n", alpha);
	printf("gamma: %lf\n", gamma);
}