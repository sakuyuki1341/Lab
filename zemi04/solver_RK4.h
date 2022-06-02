#ifndef SUB_H
#define SUB_H

int judge_break();
int init();
int RK4();
int Euler(double* kx, double* ky, double* kz);
int llg(double* kx, double* ky, double* kz, double* Hx_ef, double* Hy_ef, double* Hz_ef);
int Heff(double* Hx_ef, double* Hy_ef, double* Hz_ef);
int Hext(double* Hx_ef, double* Hy_ef, double* Hz_ef);
int HK(double* Hx_ef, double* Hy_ef, double* Hz_ef);
int HA(double* Hx_ef, double* Hy_ef, double* Hz_ef);
int HD(double* Hx_ef, double* Hy_ef, double* Hz_ef);
int vadd(double* mx0, double* my0, double* mz0, double* kx, double* ky, double* kz, double r);
int vadd4(double* mx0, double* my0, double* mz0, double* k1x, double* k1y, double* k1z, double* k2x, double* k2y, double* k2z, double* k3x, double* k3y, double* k3z, double* k4x, double* k4y, double* k4z);

#endif 