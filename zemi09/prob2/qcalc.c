#include <math.h>
#include "qcalc.h"

double F1(double x,double y,double z,double r){
  double tmp=0.0;

  if( x==0 && y==0 && z==0 ){
    tmp=0.0;
  }else{
    tmp = (y*y+z*z-2*x*x)*r/6.0;
    if( fabs(y-r)!=0 )
      tmp += y*(z*z-x*x)*log(fabs(y-r))/2.0;
    if( fabs(z-r)!=0 )
      tmp += z*(y*y-x*x)*log(fabs(z-r))/2.0;
    if( x!=0 )
      tmp += x*y*z*atan((y*z)/(r*x));
  }
  return tmp;
}

double F2(double z, double x, double y, double r){
  double tmp=0.0;

  if(x==0&&y==0&&z==0){
    tmp = 0.0;
  }else{
    tmp = x*y*r/3.0;

    if(fabs(x+r)!=0)
      tmp += y*(y*y-3*z*z)*log(fabs(x+r))/6.0;
    if(fabs(y+r)!=0)
      tmp += x*(x*x-3*z*z)*log(fabs(y+r))/6.0;
    if(fabs(z+r)!=0)
      tmp +=-x*y*z*log(fabs(z+r));
    if(x!=0)
      tmp += x*x*z*atan((y*z)/(r*x))/2.0;
    if(y!=0)
      tmp += y*y*z*atan((z*x)/(r*y))/2.0;
    if(z!=0)
      tmp += z*z*z*atan((y*x)/(r*z))/6.0;
  }
  return tmp;
}

void Q_Calc(double  dxo, double  dyo, double  dzo, 
	    double  dxs, double  dys, double  dzs, 
	    double   xs, double   ys, double   zs,
	    double *qxx, double *qyy, double *qzz,
	    double *qxy, double *qxz, double *qyz, double M)
{
  int i0, j0, k0, i1, j1, k1;
  double v, s, r, in_x, in_y, in_z;

  *qxx = 0.0; *qyy = 0.0; *qzz = 0.0;
  *qxy = 0.0; *qxz = 0.0; *qyz = 0.0;
  v = dxo*dyo*dzo;

  for(i1=0; i1<2; i1++)
    for(j1=0; j1<2; j1++)
      for(k1=0; k1<2; k1++)
	for(i0=0; i0<2; i0++)
	  for(j0=0; j0<2; j0++)
	    for(k0=0; k0<2; k0++){
	      s = (double) M * pow(-1,i0+j0+k0+i1+j1+k1)/v;
	      in_x = xs + (i1-0.5)*dxo - (i0-0.5)*dxs;
	      in_y = ys + (j1-0.5)*dyo - (j0-0.5)*dys;
	      in_z = zs + (k1-0.5)*dzo - (k0-0.5)*dzs;
	      r = sqrt(pow(in_x,2.0)+pow(in_y,2.0)+pow(in_z,2.0));
	      *qxx += s * F1(in_x, in_y, in_z, r);
	      *qyy += s * F1(in_y, in_x, in_z, r);
	      *qzz += s * F1(in_z, in_x, in_y, r);
	      *qxy += s * F2(in_z, in_x, in_y, r);
	      *qxz += s * F2(in_y, in_x, in_z, r);
	      *qyz += s * F2(in_x, in_y, in_z, r);
	    }
}

