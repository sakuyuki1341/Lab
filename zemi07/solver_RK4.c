#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "solver_RK4.h"

// moment[nz][nx]
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
double qxz[2047][2047];
double qzz[2047][2047];

clock_t st, ed;

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("the option is fault\n");
		printf("use: ./a.exe [Bloch-a, Bloch-b, Neel-a, Neel-b]\n");
		return 0;
	}

	st = clock();

	init();

	int i, j, k, l, t;
		// 初期条件設定
	if (strcmp(argv[1], "Bloch-a") == 0) {
		for (i = 0; i < nz; i++) {
			for (j = 0; j < nx/2; j++) {
				//myの設定
				moment[i][j].m[1] = -(dnx/2-(1+(double)j))*(2/dnx) - 1/dnx;
				moment[i][nx-j-1].m[1] = -moment[i][j].m[1];
				//mz,mxの設定
				//　膜の表面の場合
				if (i == 0) {
					moment[i][j].m[0] = sqrt(1-moment[i][j].m[1]*moment[i][j].m[1]);
					moment[i][nx-j-1].m[0] = moment[i][j].m[0];
				}else if (i == nz-1) {
					moment[i][j].m[0] = -sqrt(1-moment[i][j].m[1]*moment[i][j].m[1]);
					moment[i][nx-j-1].m[0] = moment[i][j].m[0];
				}else{
					moment[i][j].m[2] = sqrt(1-moment[i][j].m[1]*moment[i][j].m[1]);
					moment[i][nx-j-1].m[2] = moment[i][j].m[2];
				}
			}
		}
	}else if (strcmp(argv[1], "Bloch-b") == 0) {
		for (i = 0; i < nz; i++) {
			for (j = 0; j < nx/2; j++) {
				// myの設定
				moment[i][j].m[1] = -(dnx/2-(1+(double)j))*(2/dnx) - 1/dnx;
				moment[i][nx-j-1].m[1] = -moment[i][j].m[1];
				//mz,mxの設定
				// 膜の表面の場合
				if (i == 0) {
					moment[i][j].m[0] = sqrt(1-pow(moment[i][j].m[1], 2.0));
					moment[i][nx-j-1].m[0] = -moment[i][j].m[0];
				}else if (i == nz-1) {
					moment[i][j].m[0] = -sqrt(1-pow(moment[i][j].m[1], 2.0));
					moment[i][nx-j-1].m[0] = -moment[i][j].m[0];
				}else{
					moment[i][j].m[2] = sqrt(1-moment[i][j].m[1]*moment[i][j].m[1]);
					moment[i][nx-j-1].m[2] = moment[i][j].m[2];
				}
			}
		}
	}else if (strcmp(argv[1], "Neel-a") == 0) {
		for (i = 0; i < nz; i++) {
			for (j = 0; j < nx/2; j++) {
				// myの設定
				moment[i][j].m[1] = -(dnx/2-(1+(double)j))*(2/dnx) - 1/dnx;
				moment[i][nx-j-1].m[1] = -moment[i][j].m[1];
				moment[i][j].m[0] = sqrt(1-moment[i][j].m[1]*moment[i][j].m[1]);
				moment[i][nx-j-1].m[0] = moment[i][j].m[0];
			}
		}
	}else if (strcmp(argv[1], "Neel-b") == 0) {
		for (i = 0; i < nz; i++) {
			for (j = 0; j < nx/2; j++) {
				// myの設定
				moment[i][j].m[1] = -(dnx/2-(1+(double)j))*(2/dnx) - 1/dnx;
				moment[i][nx-j-1].m[1] = -moment[i][j].m[1];
				if (i == 0 || i == nz-1) {
					moment[i][j].m[0] = sqrt(1-moment[i][j].m[1]*moment[i][j].m[1]);
					moment[i][nx-j-1].m[0] = moment[i][j].m[0];
				}else{
					moment[i][j].m[0] = 0.9;
					moment[i][j].m[2] = 0.5;
					// 正規化
					double absm = sqrt(pow(moment[i][j].m[0], 2.0) + pow(moment[i][j].m[1], 2.0) + pow(moment[i][j].m[2], 2.0));
					for (l = 0; l < 3; l++) {
						moment[i][j].m[l] = moment[i][j].m[l]/absm;
					}
				}
			}
		}
	}else{
		printf("the option is fault\n");
		printf("use: ./a.exe [Bloch-a, Bloch-b, Neel]\n");
		return 0;
	}

	// 平衡状態計算
	for (t = 0; t < loops; t++) {
		RK4();
		// 収束判定
		if (judge_break1()) {
			break;
		}
		if (t == loops) {
			printf("timeout\n");
		}	
		if (t%plots == 0) {
			//tester(1,"m");
		}
	}
	tester(1,"g");
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
	for (i = 0; i < nz; i++) {
		for (j = 0; j < nx; j++) {
			// 可読性のため一時的に別変数へ
			double* m = moment[i][j].m;
			double* H_ef = moment[i][j].H_ef;
			avgTmp += m[0]*H_ef[0] + m[1]*H_ef[1] + m[2]*H_ef[2];
		}
	}
	avgTmp = avgTmp/(nx*nz);
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

	// qxx, qzzの計算
	double k, l;
	for (k = 0; k < nz; k++) {
		for (l = 0; l < nx; l++) {
			qxx[(int)k+(nz-1)][(int)l+(nx-1)] = -2*M*(atan((k+0.5)*dz/((l+0.5)*dx)) - atan((k+0.5)*dz/((l-0.5)*dx)) - 
			  atan((k-0.5)*dz/((l+0.5)*dx)) + atan((k-0.5)*dz/((l-0.5)*dx)));
			qxx[(int)k+(nz-1)][-(int)l+(nx-1)] = qxx[(int)k+(nz-1)][(int)l+(nx-1)];
			qxx[-(int)k+(nz-1)][(int)l+(nx-1)] = qxx[(int)k+(nz-1)][(int)l+(nx-1)];
			qxx[-(int)k+(nz-1)][-(int)l+(nx-1)] = qxx[(int)k+(nz-1)][(int)l+(nx-1)];

			qzz[(int)k+(nz-1)][(int)l+(nx-1)] = -2*M*(atan((l+0.5)*dx/((k+0.5)*dz)) - atan((l-0.5)*dx/((k+0.5)*dz)) -
			  atan((l+0.5)*dx/((k-0.5)*dz)) + atan((l-0.5)*dx/((k-0.5)*dz)));
			qzz[(int)k+(nz-1)][-(int)l+(nx-1)] = qzz[(int)k+(nz-1)][(int)l+(nx-1)];
			qzz[-(int)k+(nz-1)][(int)l+(nx-1)] = qzz[(int)k+(nz-1)][(int)l+(nx-1)];
			qzz[-(int)k+(nz-1)][-(int)l+(nx-1)] = qzz[(int)k+(nz-1)][(int)l+(nx-1)];
		}
	}

	// qxzの計算
	for (k = 0; k < 2*nz-1; k++) {
		for (l = 0; l < 2*nx-1; l++) {
			double km = k-(nz-1);
			double lm = l-(nx-1);
			qxz[(int)k][(int)l] = -M*log(pow(((lm+0.5)*dx),2.0)+pow(((km+0.5)*dz),2.0)) + M*log(pow(((lm-0.5)*dx),2.0)+pow(((km+0.5)*dz),2.0)) +
			  M*log(pow(((lm+0.5)*dx),2.0)+pow(((km-0.5)*dz),2.0)) - M*log(pow(((lm-0.5)*dx),2.0)+pow(((km-0.5)*dz),2.0));
		}
	}

	// 各位置での(z,x)の値
	int i, j;
	for (i = 0; i < nz; i++) {
		double tmp_z = -(dnz/2 - (1+(double)i))*dz - dz/2;
		for (j = 0; j < nx/2; j++) {
			moment[i][j].x = -(dnx/2 - (1+(double)j))*dx - dx/2;
			moment[i][nx-j-1].x = -moment[i][j].x;
			moment[i][j].z = tmp_z;
			moment[i][nx-j-1].z = tmp_z;
		}
	}
	return 0;
}

// ---------------------------------------
//	ルンゲクッタ法
// ---------------------------------------
int RK4() {
	int i, j, k;
	for (i = 0; i < nz; i++) {
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
	for (i = 0; i < nz; i++) {
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
	for (i = 0; i < nz; i++) {
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
	for (i = 0; i < nz; i++) {
		for (j = 0; j < nx; j++) {
			moment[i][j].H_ef[2] += 2 * ku * moment[i][j].m[2] / M;
		}
	}
	return 0;
}

int HA() {
	int i, j, k;
	double mp, mm;
	for (i = 0; i < nz; i++) {
		for (j = 0; j <nx; j++) {
			for (k = 0; k < 3; k++) {
				HA_sub_x(i, j, k, &mp, &mm);
				moment[i][j].H_ef[k] += 2*A*(mp - 2*moment[i][j].m[k] + mm) / (M*dx*dx);
				HA_sub_z(i, j, k, &mp, &mm);
				moment[i][j].H_ef[k] += 2*A*(mp - 2*moment[i][j].m[k] + mm) / (M*dz*dz);
			}
		}
	}
	return 0;
}

int HA_sub_x(int i, int j, int k, double* mp, double* mm) {
	if (j == 0) {
		*mp = moment[i][j+1].m[k];
		*mm = 0;
		if (j == 1) {
			*mm = -1;
		}
	} else if (j == nx-1) {
		*mp = 0;
		*mm = moment[i][j-1].m[k];
		if (j == 1) {
			*mp = 1;
		}
	} else {
		*mp = moment[i][j+1].m[k];
		*mm = moment[i][j-1].m[k];
	}
	return 0;
}

int HA_sub_z(int i, int j, int k, double* mp, double* mm) {
	if (i == 0) {
		*mp = moment[i+1][j].m[k];
		*mm = moment[i][j].m[k];

	} else if (i == nz-1) {
		*mp = moment[i][j].m[k];
		*mm = moment[i-1][j].m[k];
	} else {
		*mp = moment[i+1][j].m[k];
		*mm = moment[i-1][j].m[k];
	}
	return 0;
}

int HD() {
	int io_x, io_z, is_x, is_z;
	int ib_x = nx-1;
	int ib_z = nz-1;
	for (io_z = 0; io_z < nz; io_z++) {
		for (io_x = 0; io_x < nx; io_x++) {
			for (is_z = 0; is_z < nz; is_z++) {
				for (is_x = 0; is_x < nx; is_x++) {
					moment[io_z][io_x].H_ef[0] += qxx[is_z-io_z+ib_z][is_x-io_x+ib_x] * moment[is_z][is_x].m[0]  +
					  qxz[is_z-io_z+ib_z][is_x-io_x+ib_x] * moment[is_z][is_x].m[2];
					moment[io_z][io_x].H_ef[2] += qzz[is_z-io_z+ib_z][is_x-io_x+ib_x] * moment[is_z][is_x].m[2]  +
					  qxz[is_z-io_z+ib_z][is_x-io_x+ib_x] * moment[is_z][is_x].m[0];
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
	for (i = 0; i < nz; i++) {
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
	for (i = 0; i < nz; i++) {
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
		/*
			if (argv[i+1] == 'x') {
				printf("x(nm) mx:\n");
				for (j = 0; j < n; j++) {
					printf("%.6e %.6e\n", moment[j].x*10000000.0, moment[j].m[0]);
				}
				i += 1;
			}else if (argv[i+1] == 'y') {
				printf("x(nm) my:\n");
				for (j = 0; j < n; j++) {
					printf("%.6e %.6e\n", moment[j].x*10000000.0, moment[j].m[1]);
				}
				i += 1;
			}else if (argv[i+1] == 'z') {
				printf("x(nm) mz:\n");
				for (j = 0; j < n; j++) {
					printf("%.6e %.6e\n", moment[j].x*10000000.0, moment[j].m[2]);
				}
				i += 1;
			}else{*/
				printf("mx my mz:\n");
				for (j = 0; j < nz; j++) {
					for (k = 0; k < nx; k++) {
						printf("%.6e %.6e %.6e\n", moment[j][k].m[0], moment[j][k].m[1], moment[j][k].m[2]);
					}
				}
//			}
			printf("-----------------------------------------\n");
			break;
/*
		case 'H':
			printf("Hx Hy Hz:\n");
			for (j = 0; j < n; j++) {
				printf("%.6e %.6e %.6e\n", moment[j].H[0], moment[j].H[1], moment[j].H[2]);
			}
			printf("-----------------------------------------\n");
			break;

		case 'e':
			printf("Hx_ef Hy_ef Hz_ef:\n");
			for (j = 0; j < n; j++) {
				printf("%.6e %.6e %.6e\n", moment[j].H_ef[0], moment[j].H_ef[1], moment[j].H_ef[2]);
			}
			printf("-----------------------------------------\n");
			break;

		case '1':
			printf("k1x k1y k1z:\n");
			for (j = 0; j < n; j++) {
				printf("%.6e %.6e %.6e\n", moment[j].k1[0], moment[j].k1[1], moment[j].k1[2]);
			}
			printf("-----------------------------------------\n");
			break;

		case '2':
			printf("k2x k2y k2z:\n");
			for (j = 0; j < n; j++) {
				printf("%.6e %.6e %.6e\n", moment[j].k2[0], moment[j].k2[1], moment[j].k2[2]);
			}
			printf("-----------------------------------------\n");
			break;

		case '3':
			printf("k3x k3y k3z:\n");
			for (j = 0; j < n; j++) {
				printf("%.6e %.6e %.6e\n", moment[j].k3[0], moment[j].k3[1], moment[j].k3[2]);
			}
			printf("-----------------------------------------\n");
			break;

		case '4':
			printf("k4x k4y k4z:\n");
			for (j = 0; j < n; j++) {
				printf("%.6e %.6e %.6e\n", moment[j].k4[0], moment[j].k4[1], moment[j].k4[2]);
			}
			printf("-----------------------------------------\n");
			break;
*/

		// ここでqxx, qzz両方処理する。
		case 'q':
			if (argv[i+1] == 'x') {
				printf("k l qxx[k][l]:\n");
				for (j = 0; j < 2*nz-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e\n", j-nz+1, k-nx+1, qxx[j][k]);
					}
				}
				i += 1;
			}else if (argv[i+1] == 'z') {
				printf("k l qzz[k][l]:\n");
				for (j = 0; j < 2*nz-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e\n", j-nz+1, k-nx+1, qzz[j][k]);
					}
				}
				i += 1;
			}else if (argv[i+1] == '_') {
				printf("k l qxz[k][l]:\n");
				for (j = 0; j < 2*nz-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e\n", j-nz+1, k-nx+1, qxz[j][k]);
					}
				}
				i += 1;
			}else{
				printf("k l qxx qzz qxz:\n");
				for (j = 0; j < 2*nz-1; j++) {
					for (k = 0; k < 2*nx-1; k++) {
						printf("%d %d %.6e %.6e %.6e\n", j-nz+1, k-nx+1, qxx[j][k], qzz[j][k], qxz[j][k]);
					}
				}
			}
			printf("-----------------------------------------\n");
			break;

		case 'g':
			printf("#x(nm) z(nm) vx vz:\n");
			for (j = 0; j < nz; j++) {
				for (k = 0; k < nx; k++) {
					printf("%.6e %.6e %.6e %.6e\n", moment[j][k].x*10000000, moment[j][k].z*10000000, moment[j][k].m[0]*10, moment[j][k].m[2]*10);
				}
			}
			printf("#-----------------------------------------\n");
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
