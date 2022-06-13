#ifndef SUB_H
#define SUB_H

struct m_moment {
	double m[3];
	double m0[3];
	double H[3];
	double H_ef[3];
	double k1[3];
	double k2[3];
	double k3[3];
	double k4[3];
	double x;
};
typedef struct m_moment M_moment;

int judge_break1();
int judge_break2();
int init();
int RK4();
int Euler(int target);
int llg(int target);
double* k_sub(int i, int target);
int Heff();
int Hext();
int HK();
int HA();
int HA_sub(int i, int j, double* mp, double* mm);
int HD();
int vadd(int target, double r);
int vadd4();
double calc_mid();
double calc_v();
int tester(int argc, char argv[]);
int print_v(double t);
int print_phi(double t);
int print_mid(double t);
int save_data(char* filename);
int load_data(char* filename);

#endif 