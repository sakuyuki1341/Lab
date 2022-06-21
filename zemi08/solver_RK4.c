#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "solver_RK4.h"

// moment[ny][nx]
M_moment moment[1024][1024];

double dt = 0;
double dx = 80e-8;
double dy = 80e-8;
double dz = 0;
double alpha = 0;
double Gamma = 0;
double A = 0;
double ku = 0;
double M = 0;
double lw = 0;
double phi = 0;
double loops = 0;
int plots = 0;
int nx = 48;
int ny = 16;
int nz = 0;
double dnx = 48.0;
double dny = 16.0;
double dnz = 0;

// テスト用関数
//  in: (オプション文字の数, オプション名一文で)
// out: なし
// オプション名一覧
//	・c: 材料定数の出力
//	・m: 磁気モーメントの出力
//	・mx:座標系とmxの出力
//	・my:座標系とmyの出力
//	・mz:座標系とmzの出力
//	・H: 外部磁界の出力
//	・e: 実効磁界の出力
//	・1: k1の出力
//	・2: k2の出力
//	・3: k3の出力
//	・4: k4の出力
//	・q: qxx,qzz,qxzの出力
//	・qx: qxxの出力
//	・qz: qzzの出力
//	・q_: qxzの出力
//	・g: gnuplot用のベクトル出力(x-z平面)
//	・t: 計算時間の出力
int tester(int argc, char* argv) {
	int i, j, k;
	for (i = 0; i < argc; i++) {
		switch (argv[i]) {
		case 'c':
			printf("      dt = %.8e\n", dt);
			printf("      dx = %.8e\n", dx);
			printf("      dz = %.8e\n", dz);
			printf("   alpha = %.8e\n", alpha);
			printf("   Gamma = %.8e\n", Gamma);
			printf("       A = %.8e\n", A);
			printf("      ku = %.8e\n", ku);
			printf("       M = %.8e\n", M);
			printf("   loops = %.8e\n", loops);
			printf("   plots = %d\n", plots);
			printf("      nx = %d\n", nx);
			printf("      ny = %d\n", ny);
			printf("      nz = %d\n", nz);
			printf("-----------------------------------------\n");
			break;

		case 'm':
				printf("mx my mz:\n");
				for (j = 0; j < ny; j++) {
					for (k = 0; k < nx; k++) {
						printf("%.6e %.6e %.6e\n", moment[j][k].m[0], moment[j][k].m[1], moment[j][k].m[2]);
					}
				}

			printf("-----------------------------------------\n");
			break;

		case 'g':
			printf("#x(nm) y(nm) vx vz:\n");
			for (j = 0; j < ny; j++) {
				for (k = 0; k < nx; k++) {
					printf("%.6e %.6e %.6e %.6e\n", moment[j][k].x*10000000, moment[j][k].y*10000000, moment[j][k].m[0]*10, moment[j][k].m[1]*10);
				}
			}
			printf("#-----------------------------------------\n");
			break;

		default:
			break;
		}
	}
}

int main(int argc, char *argv[]) {
	// 各位置での(z,x)の値
	int i, j, l;
	for (i = 0; i < ny; i++) {
		double tmp_y = -(dny/2 - (1+(double)i))*dy - dy/2;
		for (j = 0; j < nx/2; j++) {
			moment[i][j].x = -(dnx/2 - (1+(double)j))*dx - dx/2;
			moment[i][nx-j-1].x = -moment[i][j].x;
			moment[i][j].y = tmp_y;
			moment[i][nx-j-1].y = tmp_y;
		}
	}

	for (i = 0; i  < ny; i++) {
		for (j = 0; j < nx; j++) {
			moment[i][j].m[0] = -moment[i][j].y;
			moment[i][j].m[1] = moment[i][j].x;
			moment[i][j].m[2] = 0;
			double absm = sqrt(pow(moment[i][j].m[0], 2.0) + pow(moment[i][j].m[1], 2.0) + pow(moment[i][j].m[2], 2.0));
			for (l = 0; l < 3; l++) {
				moment[i][j].m[l] = moment[i][j].m[l]/absm;
			}
		}
	}

	tester(1,"m");
	return 0;
}