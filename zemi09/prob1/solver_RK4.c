#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "solver_RK4.h"

M_moment moment[1024];

double dt = 0;
double dy = 0;
double dx = 0;
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
int n = 0;
double dn = 0;

double qxx[2047];
double qxz = 0;
double qzz[2047];

int main(int argc, char *argv[]) {
	if(argc < 3) {
		printf("the option is fault\n");
		printf("use: ./a.exe [save, load] [Bloch, Neel]\n");
		return 0;
	}
	init();

	int i, t;
	// 平衡状態計算
	if (strcmp(argv[1], "save") == 0) {
	//平衡状態セーブ
		// 初期条件設定
		for (i = 0; i < n/2; i++) {
			//myの設定
			moment[i].m[1] = -(dn/2-(1+(double)i))*(2/dn) - 1/dn;
			moment[n-i-1].m[1] = -moment[i].m[1];
			//mzの設定
			if (strcmp(argv[2], "Bloch") == 0) {
				moment[i].m[2] = sqrt(1-moment[i].m[1]*moment[i].m[1]);
				moment[n-i-1].m[2] = moment[i].m[2];
			}else if (strcmp(argv[2], "Neel") == 0) {
				moment[i].m[0] = sqrt(1-moment[i].m[1]*moment[i].m[1]);
				moment[n-i-1].m[0] = moment[i].m[0];
			}else{
				printf("use second option\n[Bloch, Neel]\n");
				return 0;
			}
		}

		// 平衡状態計算
		for (t = 0; t < loops; t++) {
			RK4();
			// 収束判定
			if (judge_break1()) {
				break;
			} else if (t == loops) {
				printf("timeout\n");
			}
			
		}
		tester(2,"ml");

		save_data("data.bin");
		printf("data is saved\ncalculate ended\n");
		return 0;

	}else if (strcmp(argv[1], "load") == 0){
		load_data("data.bin");
		printf("data is loaded\n");
		return 0;

	}else{
		printf("use option [save, load]\n");
		return 0;
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
	int i;
	for (i = 0; i < n; i++) {
		// 可読性のため一時的に別変数へ
		double* m = moment[i].m;
		double* H_ef = moment[i].H_ef;
		avgTmp += m[0]*H_ef[0] + m[1]*H_ef[1] + m[2]*H_ef[2];
	}
	avgTmp = avgTmp/n;
	count += 1;

	if(IsFirst) {
		avgFirst = avgTmp;
		avgLast = avgTmp;
		IsFirst = false;
		return 0;
	} else {
		if (fabs((avgLast-avgTmp)*1.0e12) <= fabs(avgFirst-avgTmp)) {
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
			sscanf(line, "%s %d", str, &n);
			dn = (double)n;
			break;
		default:
			break;
		}
	}

	// qxx, qzzの計算
	double k;
	for (k = 0; k < n; k++) {
		qxx[(int)k+(n-1)] = -2*M*(atan(0.5*dz/((k+0.5)*dx)) - atan(0.5*dz/((k-0.5)*dx)) - 
		  atan(-0.5*dz/((k+0.5)*dx)) + atan(-0.5*dz/((k-0.5)*dx)));
		qxx[-(int)k+(n-1)] = qxx[(int)k+(n-1)];

		qzz[(int)k+(n-1)] = -2*M*(atan((k+0.5)*dx/(0.5*dz)) - atan((k-0.5)*dx/(0.5*dz)) -
		  atan((k+0.5)*dx/(-0.5*dz)) + atan((k-0.5)*dx/(-0.5*dz)));
		qzz[-(int)k+(n-1)] = qzz[(int)k+(n-1)];
	}

	// 各位置でのxの値
	int i;
	for (i = 0; i < n/2; i++) {
		moment[i].x = -(dn/2 - (1+(double)i))*dx - dx/2;
		moment[n-i-1].x = -moment[i].x;
	}
	return 0;
}


// ---------------------------------------
//	ルンゲクッタ法
// ---------------------------------------
int RK4() {
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < 3; j++) {
			moment[i].m0[j] = moment[i].m[j];
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
	int i;
	for (i = 0; i < n; i++) {
		double* k = k_sub(i, target);
		double* m = moment[i].m;
		double* H_ef = moment[i].H_ef;
		double c = m[0]*H_ef[0] + m[1]*H_ef[1] + m[2]*H_ef[2];
		k[0] = dt * (-fabs(Gamma)) * (m[1]*H_ef[2] - m[2]*H_ef[1] + alpha*(c*m[0] - H_ef[0])) / (1 + alpha*alpha);
		k[1] = dt * (-fabs(Gamma)) * (m[2]*H_ef[0] - m[0]*H_ef[2] + alpha*(c*m[1] - H_ef[1])) / (1 + alpha*alpha);
		k[2] = dt * (-fabs(Gamma)) * (m[0]*H_ef[1] - m[1]*H_ef[0] + alpha*(c*m[2] - H_ef[2])) / (1 + alpha*alpha);
	}
	return 0;
}

double* k_sub(int i, int target) {
	switch (target) {
	case 1:
		return moment[i].k1;
		break;
	case 2:
		return moment[i].k2;
		break;
	case 3:
		return moment[i].k3;
		break;
	case 4:
		return moment[i].k4;
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
	int i;
	for (i = 0; i < n; i++) {
		moment[i].H_ef[0] = 0;
		moment[i].H_ef[1] = 0;
		moment[i].H_ef[2] = 0;
	}
	return 0;
}

int HK() {
	int i;
	for (i = 0; i < n; i++) {
		moment[i].H_ef[2] += 2 * ku * moment[i].m[2] / M;
	}
	return 0;
}

int HA() {
	int i, j;
	double mp, mm;
	for (i = 0; i < n; i++) {
		for (j = 0; j < 3; j++) {
			HA_sub(i, j, &mp, &mm);
			moment[i].H_ef[j] += 2*A*(mp - 2*moment[i].m[j] + mm) / (M*dx*dx);
		}
	}
	return 0;
}

int HA_sub(int i, int j, double* mp, double* mm) {
	if (i == 0) {
		*mp = moment[i+1].m[j];
		*mm = 0;
		if (j == 1) {
			*mm = -1;
		}
	} else if (i == n-1) {
		*mp = 0;
		*mm = moment[i-1].m[j];
		if (j == 1) {
			*mp = 1;
		}
	} else {
		*mp = moment[i+1].m[j];
		*mm = moment[i-1].m[j];
	}
	return 0;
}

int HD() {
	int io, is;
	int ib = n-1;
	for (io = 0; io < n; io++) {
		for (is = 0; is < n; is++) {
			moment[io].H_ef[0] += qxx[is-io+ib]*moment[is].m[0] + qxz*moment[is].m[2];
			moment[io].H_ef[2] += qxz*moment[is].m[0] 			+ qzz[is-io+ib]*moment[is].m[2];
		}
	}
	
	return 0;
}


// ---------------------------------------
//	vadd 
// ---------------------------------------
int vadd(int target, double r) {
	int i, j;
	for (i = 0; i < n; i++) {
		double* k = k_sub(i, target);
		for (j = 0; j < 3; j++) {
			moment[i].m[j] = moment[i].m0[j] + k[j]*r;
		}
		// 正規化
		double absm = sqrt(moment[i].m[0]*moment[i].m[0] + moment[i].m[1]*moment[i].m[1] + moment[i].m[2]*moment[i].m[2]);
		for (j = 0; j < 3; j++) {
			moment[i].m[j] = moment[i].m[j]/absm;
		}
	}
	return 0;
}

int vadd4() {
	int i, j;
	for (i = 0; i < n; i++) {
		double* k1 = moment[i].k1;
		double* k2 = moment[i].k2;
		double* k3 = moment[i].k3;
		double* k4 = moment[i].k4;

		for (j = 0; j < 3; j++) {
			moment[i].m[j] = moment[i].m0[j] + (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])/6;
		}
		// 正規化
		double absm = sqrt(moment[i].m[0]*moment[i].m[0] + moment[i].m[1]*moment[i].m[1] + moment[i].m[2]*moment[i].m[2]);
		for (j = 0; j < 3; j++) {
			moment[i].m[j] = moment[i].m[j]/absm;
		}
	}
	return 0;
}


// ---------------------------------------
//	各種値を返す関数
// ---------------------------------------

// 磁壁の中心の座標を返す関数
double calc_mid() {
	int i;
	for (i = 0; i < n; i++) {
		if (moment[i].m[1] < 0) {
			//printf("k: %.6e\nl: %.6e\n", moment[i-1].x, moment[i].x);
			return moment[i].x - dx * fabs(moment[i].m[1])/(fabs(moment[i].m[1])+moment[i-1].m[1]);
		}
	}
}

// 磁壁幅を返す関数
double calc_lw() {
	int i;
	for (i = 0; i < n; i++) {
		if (moment[i].m[1] > 0) {
			double grad = (moment[i].m[1] - moment[i-1].m[1])/dx;
			return 2/grad;
		}
	}
}


// ---------------------------------------
//	出力用関数
// ---------------------------------------

// テスト用関数
//  in: (オプション数, オプション名一文で)
// out: なし
// オプション名一覧
//	・c: 材料定数出力
//	・m: 磁気モーメント出力
//	・H: 外部磁界出力
//	・e: 実効磁界も含めた磁界の出力
//	・1: k1の出力
//	・2: k2の出力
//	・3: k3の出力
//	・4: k4の出力
//	・q: qxxとqzzの出力
//	・qx: qxxの出力
//	・qz: qzzの出力
int tester(int argc, char* argv) {
	int i, j;
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
			printf("       n = %d\n", n);
			printf("-----------------------------------------\n");
			break;

		case 'm':
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
			}else{
				printf("mx my mz:\n");
				for (j = 0; j < n; j++) {
					printf("%.6e %.6e %.6e\n", moment[j].m[0], moment[j].m[1], moment[j].m[2]);
				}
			}
			printf("-----------------------------------------\n");
			break;

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

		// ここでqxx, qzz両方処理する。
		case 'q':
			if (argv[i+1] == 'x') {
				printf("k qxx:\n");
				for (j = 0; j < 2*n-1; j++) {
					printf("%d %.6e\n", j-n+1, qxx[j]);
				}
				i += 1;
			}else if (argv[i+1] == 'z') {
				printf("k qzz:\n");
				for (j = 0; j < 2*n-1; j++) {
					printf("%d %.6e\n", j-n+1, qzz[j]);
				}
				i += 1;
			}else{
				printf("k qxx qzz:\n");
				for (j = 0; j < 2*n-1; j++) {
					printf("%d %.6e %.6e\n", j-n+1, qxx[j], qzz[j]);
				}
			}
			printf("-----------------------------------------\n");
			break;

		case 'l':;
			// 磁壁の中心を挟む二点から傾きを求める->磁壁を求める。
			printf("simlw(nm) = %.6e\n", calc_lw()*1e7);
			printf("-----------------------------------------\n");
			break;

/*	調整中
		case 't':
			printf("x, theta:\n");
			for (i = 0; i < interval*region; i++) {
				double x = ((double)i-interval)*dx + dx/2;
				double theta = acos(moment[i].m[2]);
				printf("%.8lf %lf\n", x, theta);
			}
			printf("-----------------------------------------\n");
			break;
*/

		default:
			break;
		}
	}
}

// データファイル書き込み関数
int save_data(char* filename) {
	FILE *to;
	double buf[3072];
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < 3; j++) {
			buf[3*i + j] = moment[i].m[j];
		}
	}

	to = fopen(filename, "wb");
	fwrite(buf, sizeof(double), 3072, to);
	fclose(to);

	return 0;
}

// ---------------------------------------
//	入力用関数
// ---------------------------------------
// データファイル読み込み関数
int load_data(char* filename) {
	FILE *from;
	double buf[3072];
	from = fopen(filename, "rb");
	fread(buf, sizeof(double), 3072, from);
	fclose(from);
	
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < 3; j++) {
			moment[i].m[j] = buf[3*i + j];
		}
	}
	return 0;
}
