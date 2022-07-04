#ifndef SUB_H
#define SUB_H

// 各セルの情報
struct m_moment {
	double m[3];
	double m0[3];
	double H[3];
	double H_ef[3];
	double Hd_std[3];
	double Hd_fftw[3];
	double k1[3];
	double k2[3];
	double k3[3];
	double k4[3];
	double x;
	double y;
	double z;
};
typedef struct m_moment M_moment;

// 関数群
int judge_break1();
int init();
int RK4();
int Euler(int target);
int llg(int target);
double* k_sub(int i, int j, int target);
int Heff();
int Hext();
int HK();
int HA();
int HA_sub_x(int i, int j, int k, double* mp, double* mm);
int HA_sub_y(int i, int j, int k, double* mp, double* mm);
int HA_sub_z(int i, int j, int k, double* mp, double* mm);
int HD();
int HD_std();
int HD_fftw();
int vadd(int target, double r);
int vadd4();
int tester(int argc, char argv[]);

#endif 