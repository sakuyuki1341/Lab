#ifndef SUB_H
#define SUB_H

int RK4(double* m, double* H, double dt, double alpha, double gamma, double ku);
int Euler(double* m, double* H, double* k, double dt, double alpha, double gamma, double ku);
int llg(double* k, double mx, double my, double mz, double Hx, double Hy, double Hz, double dt, double alpha, double gamma);
int vadd(double* m, double* m0, double* k, double r);
int vadd4(double* m, double* m0, double* k1, double* k2, double* k3, double* k4);
int Heff(double* m, double* H, double* Hx, double* Hy, double* Hz, double ku);
int Hext(double* H, double* Hx, double* Hy, double* Hz);
int Hanis(double* m, double* H, double* Hx, double* Hy, double* Hz, double ku);
int printout3(double* m, int t);
int printout2(double* m, int t, double dt);
int init(double* m, double* H, double* dt, double* alpha, double* gamma, double* ku, double* loops, int* plots);

#endif 