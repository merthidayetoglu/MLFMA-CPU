#include "vars.h"

extern int numunk;
extern int max_it;
extern double res_tol;
extern double mem;
extern int nummatvec;
extern double *hist;
complex<double> *buff;

extern complex<double> *xproc;
extern complex<double> *rproc;

double norm2(complex<double> *a){
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  double ret = 0;
  for(int m = myid*numunk/numproc; m < (myid+1)*numunk/numproc; m++)
    ret = ret + norm(a[m]);
  double rettot;
  MPI_Allreduce(&ret,&rettot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return rettot;
}
complex<double> inner(complex<double> *a, complex<double> *b){
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  complex<double> red = 0;
  for(int m = myid*numunk/numproc; m < (myid+1)*numunk/numproc; m++)
    red = red + conj(a[m])*b[m];
  complex<double> redtot;
  MPI_Allreduce(&red,&redtot,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
  return redtot;
}
void saxpy(complex<double> *a, complex<double> *b, complex<double> *c, complex<double> alpha){//c = a + alpha*b
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  #pragma omp parallel for
  for(int n = myid*numunk/numproc; n < (myid+1)*numunk/numproc; n++)
    c[n] = a[n] + alpha*b[n];
}
void matvec(complex<double> *x, complex<double> *o, complex<double> *r){
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  #pragma omp parallel for
  for(int n = myid*numunk/numproc; n < (myid+1)*numunk/numproc; n++)
    buff[n] = o[n]*x[n];
  mlfma(buff,r);
  #pragma omp parallel for
  for(int n = myid*numunk/numproc; n < (myid+1)*numunk/numproc; n++)
    r[n] = x[n]-r[n];
}
void bicgs(complex<double> *x, complex<double> *o, complex<double> *b){
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  //ALLOCATIONS
  complex<double> *p = new complex<double>[numunk];
  complex<double> *v = new complex<double>[numunk];
  complex<double> *s = new complex<double>[numunk];
  complex<double> *t = new complex<double>[numunk];
  complex<double> *r = new complex<double>[numunk];
  complex<double> *r_tld = new complex<double>[numunk];
  buff = new complex<double>[numunk];
  double mems = 0;
  mems = mems + (double)sizeof(complex<double>)*7*numunk/1024/1024;
  if(myid==0)printf("BICGS SOLVER MEM: %f MB\n",mems);
  if(myid==0)printf("TOTAL MEMORY: %f MB\n",mem+mems);

  complex<double> alpha;
  complex<double> beta;
  complex<double> omega;
  complex<double> rho;
  complex<double> rho_1;

  int iter = 0;
  double rnrm2;
  double snrm2;
  double bnrm2;
  double error;

  bnrm2 = sqrt(norm2(b));
  if(bnrm2 < 1e-20)
    bnrm2 = 1;
  //MATVEC 0 START
  matvec(x,o,r);
  //MATVEC 0 FINISH
  saxpy(b,r,r,complex<double>(-1,0));
  rnrm2 = sqrt(norm2(r));
  error = rnrm2/bnrm2;
  if(myid==0)printf("RES. ERROR: %e ITER: %d\n",error,iter);
  hist[iter]=error;
  if(error > res_tol){
    memcpy(&r_tld[myid*numunk/numproc],&r[myid*numunk/numproc],numunk/numproc*sizeof(complex<double>));
    //BEGIN ITERATIONS
    while(iter < max_it){
      iter++;
      rho = inner(r_tld,r);
      if(abs(rho) < 1e-20){
        printf("RHO BREAKDOWN\n");
        break;
      }
      if(iter == 1)
        memcpy(&p[myid*numunk/numproc],&r[myid*numunk/numproc],numunk/numproc*sizeof(complex<double>));
      else{
        beta = (rho/rho_1)*(alpha/omega);
        #pragma omp parallel for
        for(int n = myid*numunk/numproc; n < (myid+1)*numunk/numproc; n++)
          p[n] = r[n] + beta*(p[n]-omega*v[n]);
      }
      //PRECONDITIONER
      //MATVEC 1 START
      matvec(p,o,v);
      //MATVEC 1 FINISH
      alpha = rho/inner(r_tld,v);
      saxpy(r,v,s,alpha*complex<double>(-1,0));

      //snrm2 = sqrt(norm2(s));
      //if(snrm2/bnrm2 < res_tol){
      //  printf("EARLY CONVERGENCE!\n");
      //  saxpy(x,p,x,alpha);
      //  error = snrm2/bnrm2;
      //  break;
      //}

      //MATVEC 2 START
      matvec(s,o,t);
      //MATVEC 2 FINISH
      //STABILIZER
      omega = inner(t,s)/norm2(t);
      //UPDATE
      #pragma omp parallel for
      for(int n = myid*numunk/numproc; n < (myid+1)*numunk/numproc; n++){
        x[n] = x[n] + alpha*p[n] + omega*s[n];
        r[n] = s[n] - omega*t[n];
      }
      rnrm2 = sqrt(norm2(r));
      error = rnrm2/bnrm2;
      if(myid==0)printf("RES. ERROR: %e ITER: %d\n",error,iter);
      hist[iter]=error;
      if(error < res_tol)
        break;
      if(abs(omega) < 1e-20){
        printf("OMEGA BREAKDOWN\n");
        break;
      }
      rho_1 = rho;
    }
  }
  if(error < res_tol)
    ;//printf("CONVERGED!\n");
  else
    printf("NOT CONVERGED!\n");

  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0)printf("NUMBER OF ITERATIONS: %d\n",iter);
  if(myid==0)printf("NUMBER OF MATVECS: %d\n",nummatvec);
  if(myid==0)printf("RESIDUAL ERROR NORM: %e\n",error);

  delete[] p;
  delete[] v;
  delete[] s;
  delete[] t;
  delete[] r_tld;
  delete[] r;
  delete[] buff;
}
