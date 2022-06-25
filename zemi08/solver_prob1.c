#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "solver_RK4.h"
#include "qcalc.h"

// moment[ny][nx]
M_moment moment[1024][1024];

double dt = 0;
double dx = 0;
double dy = 0;
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
int nx = 0;
int ny = 0;
int nz = 0;
double dnx = 0;
double dny = 0;
double dnz = 0;

double qxx[2047][2047];
double qyy[2047][2047];
double qzz[2047][2047];
double qxy[2047][2047];
double qxz[2047][2047];
double qyz[2047][2047];

clock_t st, ed;

int main(int argc, char *argv[]) {
	st = clock();

	init();

	int i, j, k, l;
	// 初期条件計算
	for (i = 0; i  < ny; i++) {
		for (j = 0; j < nx; j++) {
			moment[i][j].m[0] = -moment[i][j].y;
			moment[i][j].m[1] = moment[i][j].x;
			moment[i][j].m[2] = 0.3;
			double absm = sqrt(pow(moment[i][j].m[0], 2.0) + pow(moment[i][j].m[1], 2.0) + pow(moment[i][j].m[2], 2.0));
			for (l = 0; l < 3; l++) {
				moment[i][j].m[l] = moment[i][j].m[l]/absm;
			}
		}
	}
	tester(1,"q");
	return 0;


	int io_x, io_y, is_x, is_y;
	int ib_x = nx-1;
	int ib_y = ny-1;
	io_y = ny/2; is_x = nx/2; is_y = ny/2;
	printf("mx my mz\n%.6e %.6e %.6e\n", moment[is_y][is_x].m[0], moment[is_y][is_x].m[1], moment[is_y][is_x].m[2]);
	for (io_x = 0; io_x < nx; io_x++){
		moment[io_y][io_x].H_ef[0] = qxx[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[0]  +
		  qxy[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[1] + qxz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[2];
		moment[io_y][io_x].H_ef[1] = qxy[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[0]  +
		  qyy[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[1] + qyz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[2];
		moment[io_y][io_x].H_ef[2] = qxz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[0]  +
		  qyz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[1] + qzz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[2];
	}
	tester(1, "cc");
	printf("x Hx Hy Hz\n");
	for (i = 0; i < nx; i++) {
		printf("%.6e %.6e %.6e %.6e\n", moment[ny/2][i].x, moment[ny/2][i].H_ef[0], moment[ny/2][i].H_ef[1], moment[ny/2][i].H_ef[2]);
	}

	return 0;
}

// ---------------------------------------
// 収束判定
// ---------------------------------------
int judge_break1() {
	static bool IsFirst = true;
	double avgTmp = 0;
	static double avgFirst = 0;
	static double avgLast = 0;
	static int count = 0;
	int i, j;
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			// 可読性のため一時的に別変数へ
			double* m = moment[i][j].m;
			double* H_ef = moment[i][j].H_ef;
			avgTmp += m[0]*H_ef[0] + m[1]*H_ef[1] + m[2]*H_ef[2];
		}
	}
	avgTmp = avgTmp/(nx*ny);
	count += 1;

	if(IsFirst) {
		avgFirst = avgTmp;
		avgLast = avgTmp;
		IsFirst = false;
		return 0;
	} else {
		if (fabs((avgLast-avgTmp)*1.0e6) <= fabs(avgFirst-avgTmp)) {
			return 1;
		} else {
			avgLast = avgTmp;
			// 収束状況を見る用
			if (count%1000 == 0) {
				//printf("First: %lf  tmp: %lf\n", avgFirst, avgTmp);
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
	char fname[] = "init.data";
	char line[256];
	char str[16];
	
	// 設定用ファイルを開ける
	fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("init.data not open!\n");
		exit(1);
	}

	// 中身の取得
	while(fgets(line, 256, fp) != NULL) {
		switch (line[0]){
		case '#':
			break;
		case 't':
			sscanf(line, "%s %lf", str, &dt);
			break;
		case 'd':
			if (line[1] == 'x') {
				sscanf(line, "%s %lf", str, &dx);
			}else if (line[1] == 'y') {
				sscanf(line, "%s %lf", str, &dy);
			}else if (line[1] == 'z') {
				sscanf(line, "%s %lf", str, &dz);
			}else{
				;
			}
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
		case 'l':
			sscanf(line, "%s %lf", str, &loops);
			break;
		case 'p':
			sscanf(line, "%s %d", str, &plots);
			break;
		case 'n':
			if (line[1] == 'x') {
				sscanf(line, "%s %d", str, &nx);
				dnx = (double)nx;
			}else if (line[1] == 'y') {
				sscanf(line, "%s %d", str, &ny);
				dny = (double)ny;
			}else if (line[1] == 'z') {
				sscanf(line, "%s %d", str, &nz);
				dnz = (double)nz;
			}else{
				;
			}
			break;
		default:
			break;
		}
	}

	// 各位置での(z,y,x)の値を計算
	int i, j, k, l;
	for (i = 0; i < dny/2; i++) {
		double tmp_y = -(dny/2 - (1+(double)i))*dy - dy/2;
		for (j = 0; j < dnx/2; j++) {
			moment[i][j].x = -(dnx/2 - (1+(double)j))*dx - dx/2;
			moment[ny-i-1][j].x = moment[i][j].x;
			moment[i][nx-j-1].x = -moment[i][j].x;
			moment[ny-i-1][nx-j-1].x = -moment[i][j].x;

			moment[i][j].y = tmp_y;
			moment[i][nx-j-1].y = tmp_y;
			moment[ny-i-1][j].y = -tmp_y;
			moment[ny-i-1][nx-j-1].y = -tmp_y;

			moment[i][j].z = 0;
			moment[i][nx-j-1].z = 0;
			moment[ny-i-1][j].z = 0;
			moment[ny-i-1][nx-j-1].z = 0;
		}
	}

	// 対称性を考えず全て計算
	double stdx = moment[0][0].x - moment[0][nx-1].x;
	double stdy = moment[0][0].y - moment[ny-1][0].y;
	printf("#x:%.6e  y:%.6e\n", moment[0][nx-1].x, moment[ny-1][0].y);
	for (k = 0; k < 2*ny-1; k++) {
		for (l = 0; l < 2*nx-1; l++) {
			Q_Calc(dx, dy, dz, dx, dy, dz, stdx+l*dx, stdy+k*dy, 0, 
					&qxx[k][l], &qyy[k][l], &qzz[k][l],
					&qxy[k][l], &qxz[k][l], &qyz[k][l],
					M);
		}
	}

	// 静磁界係数の計算部
	// for (k = 0; k < ny; k++) {
	// 	for (l = 0; l < nx; l++) {
	// 		Q_Calc(dx, dy, dz, dx, dy, dz, moment[k][l].x-moment[0][0].x, moment[k][l].y-moment[0][0].y, moment[k][l].z, 
	// 				&qxx[k+(ny-1)][l+(nx-1)], &qyy[k+(ny-1)][l+(nx-1)], &qzz[k+(ny-1)][l+(nx-1)],
	// 				&qxy[k+(ny-1)][l+(nx-1)], &qxz[k+(ny-1)][l+(nx-1)], &qyz[k+(ny-1)][l+(nx-1)],
	// 				M);
	// 		qxx[k+(ny-1)][-l+(nx-1)] = qxx[k+(ny-1)][l+(nx-1)];
	// 		qxx[-k+(ny-1)][l+(nx-1)] = qxx[k+(ny-1)][l+(nx-1)];
	// 		qxx[-k+(ny-1)][-l+(nx-1)] = qxx[k+(ny-1)][l+(nx-1)];

	// 		qyy[k+(ny-1)][-l+(nx-1)] = qyy[k+(ny-1)][l+(nx-1)];
	// 		qyy[-k+(ny-1)][l+(nx-1)] = qyy[k+(ny-1)][l+(nx-1)];
	// 		qyy[-k+(ny-1)][-l+(nx-1)] = qyy[k+(ny-1)][l+(nx-1)];

	// 		qzz[k+(ny-1)][-l+(nx-1)] = qzz[k+(ny-1)][l+(nx-1)];
	// 		qzz[-k+(ny-1)][l+(nx-1)] = qzz[k+(ny-1)][l+(nx-1)];
	// 		qzz[-k+(ny-1)][-l+(nx-1)] = qzz[k+(ny-1)][l+(nx-1)];

	// 		qxy[k+(ny-1)][-l+(nx-1)] = qxy[k+(ny-1)][l+(nx-1)];
	// 		qxy[-k+(ny-1)][l+(nx-1)] = qxy[k+(ny-1)][l+(nx-1)];
	// 		qxy[-k+(ny-1)][-l+(nx-1)] = qxy[k+(ny-1)][l+(nx-1)];

	// 		qxz[k+(ny-1)][-l+(nx-1)] = qxz[k+(ny-1)][l+(nx-1)];
	// 		qxz[-k+(ny-1)][l+(nx-1)] = qxz[k+(ny-1)][l+(nx-1)];
	// 		qxz[-k+(ny-1)][-l+(nx-1)] = qxz[k+(ny-1)][l+(nx-1)];

	// 		qyz[k+(ny-1)][-l+(nx-1)] = qyz[k+(ny-1)][l+(nx-1)];
	// 		qyz[-k+(ny-1)][l+(nx-1)] = qyz[k+(ny-1)][l+(nx-1)];
	// 		qyz[-k+(ny-1)][-l+(nx-1)] = qyz[k+(ny-1)][l+(nx-1)];
	// 	}
	// }

	return 0;
}

// ---------------------------------------
//	ルンゲクッタ法
// ---------------------------------------
int RK4() {
	int i, j, k;
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			for (k = 0; k < 3; k++) {
				moment[i][j].m0[k] = moment[i][j].m[k];
			}
		}
	}


	Euler(1);
	vadd(1, 0.5);
	Euler(2);
	vadd(2, 0.5);
	Euler(3);
	vadd(3, 1.0);
	Euler(4);
	vadd4();
	return 0;
}

int Euler(int target) {
	Heff();
	llg(target);
	return 0;
}

int llg(int target) {
	int i, j;
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			double* k = k_sub(i, j, target);
			double* m = moment[i][j].m;
			double* H_ef = moment[i][j].H_ef;
			double c = m[0]*H_ef[0] + m[1]*H_ef[1] + m[2]*H_ef[2];
			k[0] = dt * (-fabs(Gamma)) * (m[1]*H_ef[2] - m[2]*H_ef[1] + alpha*(c*m[0] - H_ef[0])) / (1 + alpha*alpha);
			k[1] = dt * (-fabs(Gamma)) * (m[2]*H_ef[0] - m[0]*H_ef[2] + alpha*(c*m[1] - H_ef[1])) / (1 + alpha*alpha);
			k[2] = dt * (-fabs(Gamma)) * (m[0]*H_ef[1] - m[1]*H_ef[0] + alpha*(c*m[2] - H_ef[2])) / (1 + alpha*alpha);
		}
	}
	return 0;
}

double* k_sub(int i, int j, int target) {
	switch (target) {
	case 1:
		return moment[i][j].k1;
		break;
	case 2:
		return moment[i][j].k2;
		break;
	case 3:
		return moment[i][j].k3;
		break;
	case 4:
		return moment[i][j].k4;
		break;
	default:
		printf("llg_sub: error\n");
		exit(1);
		break;
	}
}


// ---------------------------------------
//	H に関する計算
// ---------------------------------------
int Heff() {
	Hext();
	HK();
	HA();
	HD();
	return 0;
}

int Hext() {
	int i, j;
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			moment[i][j].H_ef[0] = 0;
			moment[i][j].H_ef[1] = 0;
			moment[i][j].H_ef[2] = 0;
		}
	}
	return 0;
}

int HK() {
	int i, j;
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			moment[i][j].H_ef[2] += 2 * ku * moment[i][j].m[2] / M;
		}
	}
	return 0;
}

int HA() {
	int i, j, k;
	double mp, mm;
	for (i = 0; i < ny; i++) {
		for (j = 0; j <nx; j++) {
			for (k = 0; k < 3; k++) {
				HA_sub_x(i, j, k, &mp, &mm);
				moment[i][j].H_ef[k] += 2*A*(mp - 2*moment[i][j].m[k] + mm) / (M*dx*dx);
				HA_sub_y(i, j, k, &mp, &mm);
				moment[i][j].H_ef[k] += 2*A*(mp - 2*moment[i][j].m[k] + mm) / (M*dy*dy);
			}
		}
	}
	return 0;
}

int HA_sub_x(int i, int j, int k, double* mp, double* mm) {
	if (j == 0) {
		*mp = moment[i][j+1].m[k];
		*mm = moment[i][j].m[k];
	} else if (j == nx-1) {
		*mp = moment[i][j].m[k];
		*mm = moment[i][j-1].m[k];
	} else {
		*mp = moment[i][j+1].m[k];
		*mm = moment[i][j-1].m[k];
	}
	return 0;
}

int HA_sub_y(int i, int j, int k, double* mp, double* mm) {
	if (i == 0) {
		*mp = moment[i+1][j].m[k];
		*mm = moment[i][j].m[k];

	} else if (i == ny-1) {
		*mp = moment[i][j].m[k];
		*mm = moment[i-1][j].m[k];
	} else {
		*mp = moment[i+1][j].m[k];
		*mm = moment[i-1][j].m[k];
	}
	return 0;
}

int HD() {
	int io_x, io_y, is_x, is_y;
	int ib_x = nx-1;
	int ib_y = ny-1;
	for (io_y = 0; io_y < ny; io_y++) {
		for (io_x = 0; io_x < nx; io_x++) {
			for (is_y = 0; is_y < nz; is_y++) {
				for (is_x = 0; is_x < nx; is_x++) {
					moment[io_y][io_x].H_ef[0] += qxx[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[0]  +
					  qxy[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[1] + qxz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[2];
					moment[io_y][io_x].H_ef[1] += qxy[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[0]  +
					  qyy[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[1] + qyz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[2];
					moment[io_y][io_x].H_ef[2] += qxz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[0]  +
					  qyz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[1] + qzz[is_y-io_y+ib_y][is_x-io_x+ib_x] * moment[is_y][is_x].m[2];
				}
			}
		}
	}
	return 0;
}


// ---------------------------------------
//	vadd 
// ---------------------------------------
int vadd(int target, double r) {
	int i, j, l;
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			double* k = k_sub(i, j, target);
			for (l = 0; l < 3; l++) {
				moment[i][j].m[l] = moment[i][j].m0[l] + k[l]*r;
			}
			// 正規化
			double absm = sqrt(pow(moment[i][j].m[0], 2.0) + pow(moment[i][j].m[1], 2.0) + pow(moment[i][j].m[2], 2.0));
			for (l = 0; l < 3; l++) {
				moment[i][j].m[l] = moment[i][j].m[l]/absm;
			}
		}
	}
	return 0;
}

int vadd4() {
	int i, j, l;
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			double* k1 = moment[i][j].k1;
			double* k2 = moment[i][j].k2;
			double* k3 = moment[i][j].k3;
			double* k4 = moment[i][j].k4;

			for (l = 0; l < 3; l++) {
				moment[i][j].m[l] = moment[i][j].m0[l] + (k1[l] + 2*k2[l] + 2*k3[l] +k4[l])/6;
			}
			// 正規化
			double absm = sqrt(pow(moment[i][j].m[0], 2.0) + pow(moment[i][j].m[1], 2.0) + pow(moment[i][j].m[2], 2.0));
			for (l = 0; l < 3; l++) {
				moment[i][j].m[l] = moment[i][j].m[l]/absm;
			}
		}
	}
	return 0;
}

// ---------------------------------------
//	出力用関数
// ---------------------------------------

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
			if (argv[i+1] == 'c') {
				printf("#cell_x cell_y x y\n");
				for (j = 0; j < ny; j++) {
					for (k = 0; k <nx; k++) {
						printf("%d %d %.6e %.6e\n", k, j, moment[j][k].x, moment[j][k].y);
					}
				}
			}else{
				printf("      dt = %.8e\n", dt);
				printf("      dx = %.8e\n", dx);
				printf("      dy = %.8e\n", dy);
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
			}
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
			printf("#x(nm) y(nm) vx vy:\n");
			for (j = 0; j < ny; j++) {
				for (k = 0; k < nx; k++) {
					printf("%.6e %.6e %.6e %.6e\n", moment[j][k].x*100000000, moment[j][k].y*100000000, moment[j][k].m[0]*10, moment[j][k].m[1]*10);
				}
			}
			printf("#-----------------------------------------\n");
			break;

		// ここでqすべてを処理する。
		case 'q':
			if (argv[i+1] == 'x') {
				printf("k l qxx[k][l]:\n");
				for (j = 0; j < 2*ny-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e\n", k-nx+1, j-ny+1, qxx[j][k]);
					}
				}
				i += 1;
			}else if (argv[i+1] == 'y') {
				printf("k l qyy[k][l]:\n");
				for (j = 0; j < 2*ny-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e\n", k-nx+1, j-ny+1, qyy[j][k]);
					}
				}
				i += 1;
			}else if (argv[i+1] == 'z') {
				printf("k l qzz[k][l]:\n");
				for (j = 0; j < 2*ny-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e\n", k-nx+1, j-ny+1, qzz[j][k]);
					}
				}
				i += 1;
			}else if (argv[i+1] == '1') {
				printf("k l qxy[k][l]:\n");
				for (j = 0; j < 2*ny-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e\n", k-nx+1, j-ny+1, qxy[j][k]);
					}
				}
				i += 1;
			}else if (argv[i+1] == '2') {
				printf("k l qxz[k][l]:\n");
				for (j = 0; j < 2*ny-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e\n", k-nx+1, j-ny+1, qxz[j][k]);
					}
				}
				i += 1;
			}else if (argv[i+1] == '3') {
				printf("k l qyz[k][l]:\n");
				for (j = 0; j < 2*ny-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e\n", k-nx+1, j-ny+1, qyz[j][k]);
					}
				}
				i += 1;
			}else{
				printf("k l qxx qyy qzz qxy qxz qyz:\n");
				for (j = 0; j < 2*ny-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e %.6e %.6e %.6e %.6e %.6e\n", k-nx+1, j-ny+1, qxx[j][k], qyy[j][k], qzz[j][k], qxy[j][k], qxz[j][k], qyz[j][k]);
					}
				}
			}
			printf("-----------------------------------------\n");
			break;

		case 't':
			ed = clock();
			printf("#nz t(s):\n");
			printf("%d %.6e\n", nz, (double)(ed-st)/CLOCKS_PER_SEC);
			break;

		default:
			break;
		}
	}
}
