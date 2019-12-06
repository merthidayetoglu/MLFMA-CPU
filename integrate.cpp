//Integration Routines
//Mert Hidayetoglu, Oct 2015
#include "vars.h"

complex<double> hn(int order, double dist){
  return complex<double>(jn(order,dist),yn(order,dist));
}

static double integ(double x, double y){
  return -3*x*y+x*y*log(x*x+y*y)+x*x*atan(y/x)+y*y*atan(x/y);
}

complex<double> integrate(complex<double> post, complex<double> posb){

  extern double res;
  complex<double>numer(0,0);
  complex<double>anal(0,0);
  //NUMERICAL PART
  double dist = abs(post-posb);
  numer = complex<double>(j0(2*M_PI*dist),0);
  if(dist < res/8){
    numer = numer + complex<double>(0,0.5772156649015329*2/M_PI+2/M_PI*log(0.5));
  //else
    //numer = numer + complex<double>(0,y0(2*M_PI*dist)-2/M_PI*log(2*M_PI*dist));
    //ANALYTICAL PART
    double xcen=(posb-post).real();
    double ycen=(posb-post).imag();
    double xmin=xcen-res/2;
    double xmax=xcen+res/2;
    double ymin=ycen-res/2;
    double ymax=ycen+res/2;
    double analt=integ(xmax,ymax)-integ(xmin,ymax)-integ(xmax,ymin)+integ(xmin,ymin);
    anal=complex<double>(0,(analt/2+log(2*M_PI)*res*res)*2/M_PI);
  }
  else
    numer =  numer + complex<double>(0,y0(2*M_PI*dist));
  numer = numer*res*res;
  return (numer+anal)*complex<double>(0,0.25);
}

complex<double> integrate_multi(complex<double> center, complex<double> posb, int order){
  extern double res;
  extern int *numsamp;
  extern int level;
  int numang = numsamp[level-1];
  double angle = 2*M_PI*order/numang;
  complex<double> u = complex<double>(cos(angle),sin(angle));
  complex<double>numer(0,0);
  //NUMERICAL PART
  complex<double> c2s = posb-center;
  numer = exp(complex<double>(0,-2*M_PI*(u.real()*c2s.real()+u.imag()*c2s.imag())));
  return numer*res*res;
}

complex<double> integrate_local(complex<double> center, complex<double> post, int order){
  extern int *numsamp;
  extern int level;
  int numang = numsamp[level-1];
  double angle = 2*M_PI*order/numang;
  complex<double> u = complex<double>(cos(angle),sin(angle));
  complex<double>numer(0,0);
  //NUMERICAL PART
  complex<double> c2t = post-center;
  numer = exp(complex<double>(0,2*M_PI*(u.real()*c2t.real()+u.imag()*c2t.imag())));
  return numer;
}

/*
void setquad(){

  extern int quad;
  extern complex<double> *qlistp;
  extern double *qlistw;

  //1-POINT GAUSSIAN-LEGENDRE QUAD LIST
  double q1[1][2] = {+0.0000000000000000, 2.0000000000000000};
  //2-POINT GAUSSIAN-LEGENDRE QUAD LIST
  double q2[2][2] = {-0.5773502691896257, 1.0000000000000000,\
                     +0.5773502691896257, 1.0000000000000000};
  //3-POINT GAUSSIAN-LEGENDRE QUAD LIST
  double q3[3][2] = {+0.0000000000000000, 0.8888888888888888,\
                     -0.7745966692414834, 0.5555555555555556,\
                     +0.7745966692414834, 0.5555555555555556};
  //4-POINT GAUSSIAN-LEGENDRE QUAD LIST
  double q4[4][2] = {-0.3399810435848563, 0.6521451548625461,\
                     +0.3399810435848563, 0.6521451548625461,\
                     -0.8611363115940526, 0.3478548451374538,\
                     +0.8611363115940526, 0.3478548451374538};
  //5-POINT GAUSSIAN-LEGENDRE QUAD LIST
  double q5[5][2] = {+0.0000000000000000, 0.5688888888888889,\
                     -0.5384693101056831, 0.4786286704993665,\
                     +0.5384693101056831, 0.4786286704993665,\
                     -0.9061798459386640, 0.2369268850561891,\
                     +0.9061798459386640, 0.2369268850561891};
  //6-POINT GAUSSIAN-LEGENDRE QUAD LIST
  double q6[6][2] = {-0.2386191860831969, 0.4679139345726910,\
                     +0.2386191860831969, 0.4679139345726910,\
                     -0.6612093864662645, 0.3607615730481386,\
                     +0.6612093864662645, 0.3607615730481386,\
                     -0.9324695142031521, 0.1713244923791704,\
                     +0.9324695142031521, 0.1713244923791704};
  double *list = NULL;
  if(quad==1)
    list = *q1;
  else
  if(quad==4)
    list = *q2;
  else
  if(quad==9)
    list = *q3;
  else
  if(quad==16)
    list = *q4;
  else
  if(quad==25)
    list = *q5;
  else
  if(quad==36)
    list = *q6;
  else
    printf("ERROR!!! CHOOSE QUAD = 1, 4, 9, 16, 25, or 36 !!!ERROR\n");
  int ind = 0;
  for(int m = 0; m < sqrt(quad); m++)
    for(int n = 0; n < sqrt(quad); n++){
      qlistp[ind] = complex<double>(list[2*m],list[2*n]);
      qlistw[ind] = list[2*m+1]*list[2*n+1];
      ind=ind+1;
    }
}*/
