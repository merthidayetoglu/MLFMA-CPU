//Linear Solution Routines
//Mert Hidayetoglu, Oct 2015
#include "vars.h"

void mlfma(complex<double> *x, complex<double> *r){

  extern int nummatvec;
  extern int level;
  extern int box;

  extern int ninter;
  extern int *numsamp;
  extern int *numclus;

  extern int **clusnear;
  extern int **clusfar;
  extern complex<double> *near;
  extern int *neid;
  extern complex<double> *coeff_multi;
  extern complex<double> *basis_local;
  extern complex<double> **trans;
  extern int **traid;
  extern complex<double> **aggmulti;
  extern complex<double> **agglocal;
  extern double **interp;
  extern double **anterp;
  extern int **intind;
  extern int **antind;
  extern complex<double> **shiftmul;
  extern complex<double> **shiftloc;
  extern complex<double> **temp;

  extern double matve_wall;
  extern double aggrh_wall;
  extern double *m2m_w;
  extern double trans_wall;
  extern double *m2l_w;
  extern double disah_wall;
  extern double *l2l_w;
  extern double nearf_wall;
  extern double commn_wall;
  extern double commf_wall;

  extern int **sendmap;
  extern int **recvmap;
  extern int **procmap;
  extern int **clusmap;
  extern int **procint;
  extern int *clcount;
  extern complex<double> **sendbuff;
  extern complex<double> **recvbuff;

  double matve_w;
  double aggrh_w;
  double trans_w;
  double disah_w;
  double nearf_w;
  double commn_w;
  double commf_w;

  int numproc;
  int numthread;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  nummatvec++;

  MPI_Barrier(MPI_COMM_WORLD);
  matve_w = get_wall_time();

  #pragma omp parallel
  numthread = omp_get_num_threads();
  MPI_Barrier(MPI_COMM_WORLD);
  aggrh_w = get_wall_time();
  //LOWEST-LEVEL AGGREGATION
  #pragma omp parallel for
  for(int clusm = myid*numclus[level-1]/numproc; clusm < (myid+1)*numclus[level-1]/numproc; clusm++){
    int indmulti = clusm*numsamp[level-1];
    int unk = clusm*box*box;
    for(int k = 0; k < numsamp[level-1]; k++){
      complex<double> reduce = 0;
      int indbasis = k*box*box;
      for(int n = 0; n < box*box; n++)
        reduce=reduce+coeff_multi[indbasis+n]*x[unk+n];
      aggmulti[level-1][indmulti+k]=reduce;
    }
  } 
  MPI_Barrier(MPI_COMM_WORLD);
  m2m_w[level-1] = m2m_w[level-1]+(get_wall_time()-aggrh_w);
  //HIGHER-LEVEL AGGREGATIONS
  for(int i = level-2; i > 1; i--){
    double walltemp = get_wall_time();
    if(numclus[i]>numproc*numthread)
      #pragma omp parallel for
      for(int clusm = myid*numclus[i]/numproc; clusm < (myid+1)*numclus[i]/numproc; clusm++){
        int indm = clusm*numsamp[i];
        for(int km = 0; km < numsamp[i]; km++){
          complex<double> temp1 = 0;
          for(int cn = 0; cn < 4; cn++){
            int indn = (clusm*4+cn)*numsamp[i+1];
            complex<double> reduce = 0;
            //INTERPOLATE
            for(int k = 0; k < ninter; k++)
              reduce = reduce+interp[i][km*ninter+k]*aggmulti[i+1][indn+intind[i][km*ninter+k]];
            //SHIFT
            temp1 = temp1 + reduce*shiftmul[i][cn*numsamp[i]+km];
          }
          aggmulti[i][indm+km]=temp1;
        }
      }
    else//DISTRIBUTE SAMPLES
      #pragma omp parallel
      for(int clusm = myid*numclus[i]/numproc; clusm < (myid+1)*numclus[i]/numproc; clusm++){
        #pragma omp for nowait
        for(int km = 0; km < numsamp[i]; km++){
          int indm = clusm*numsamp[i];
          complex<double> temp1 = 0;
          for(int cn = 0; cn < 4; cn++){
            int indn = (clusm*4+cn)*numsamp[i+1];
            complex<double> reduce = 0;
            //INTERPOLATE
            for(int k = 0; k < ninter; k++)
              reduce = reduce+interp[i][km*ninter+k]*aggmulti[i+1][indn+intind[i][km*ninter+k]];
            //SHIFT
            temp1 = temp1 + reduce*shiftmul[i][cn*numsamp[i]+km];
          }
          aggmulti[i][indm+km]=temp1;
        }
      }
    MPI_Barrier(MPI_COMM_WORLD);
    m2m_w[i] = m2m_w[i]+(get_wall_time()-walltemp);
  }
  aggrh_w = get_wall_time()-aggrh_w;
  commf_w = get_wall_time();
  for(int i = 2; i < level; i++)
    //#pragma omp parallel for
    for(int cl = 0; cl < clcount[i]; cl++)
      memcpy(&sendbuff[i][cl*numsamp[i]],&aggmulti[i][sendmap[i][cl]*numsamp[i]],numsamp[i]*sizeof(complex<double>));
  MPI_Barrier(MPI_COMM_WORLD);
  //FARFIELD COMMUNICATION
  for(int i = 2; i < level; i++){
    for(int p = 0; p < numproc; p++){
      int sib = procmap[i][p];
      if(sib != -1){
        complex<double> *send = &sendbuff[i][procint[i][p]*numsamp[i]];
        complex<double> *recv = &recvbuff[i][procint[i][p]*numsamp[i]];
        int amount = clusmap[i][p]*numsamp[i];
        MPI_Sendrecv(send,amount,MPI_DOUBLE_COMPLEX,sib,0,recv,amount,MPI_DOUBLE_COMPLEX,sib,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(int i = 2; i < level; i++)
    //#pragma omp parallel for
    for(int cl = 0; cl < clcount[i]; cl++)
      memcpy(&aggmulti[i][recvmap[i][cl]*numsamp[i]],&recvbuff[i][cl*numsamp[i]],numsamp[i]*sizeof(complex<double>));
  MPI_Barrier(MPI_COMM_WORLD);
  commf_w = get_wall_time()-commf_w;
  trans_w = get_wall_time();
  //TRANSLATION
  for(int i = 2; i < level; i++){
    double walltemp = get_wall_time();
    #pragma omp parallel for
    for(int clusm = myid*numclus[i]/numproc; clusm < (myid+1)*numclus[i]/numproc; clusm++){
      int indlocal = clusm*numsamp[i];
      memset(&agglocal[i][indlocal],0,numsamp[i]*sizeof(complex<double>));
      for(int cn = 0; cn < 27; cn++){
        int clusn = clusfar[i][clusm*27+cn];
        if(clusn!=-1){
          int indmulti = clusn*numsamp[i];
          int index = traid[i][clusm*27+cn]*numsamp[i];
          for(int k = 0; k < numsamp[i]; k++)
            agglocal[i][indlocal+k]=agglocal[i][indlocal+k]+trans[i][index+k]*aggmulti[i][indmulti+k];    
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    m2l_w[i] = m2l_w[i]+(get_wall_time()-walltemp);
  }
  trans_w = get_wall_time()-trans_w;
  disah_w = get_wall_time();
  //HIGHER-LEVEL DISAGGREGATIONS
  for(int i = 2; i < level-1; i++){
    double walltemp = get_wall_time();
    #pragma omp parallel for
    for(int clusn = myid*numclus[i]/numproc; clusn < (myid+1)*numclus[i]/numproc; clusn++){
      int indn = clusn*numsamp[i];
      int mythread = omp_get_thread_num();
      for(int cm = 0; cm < 4; cm++){
        int indm = (clusn*4+cm)*numsamp[i+1];
        //SHIFT
        for(int kn = 0; kn < numsamp[i]; kn++)
          temp[mythread][kn] = shiftloc[i][cm*numsamp[i]+kn]*agglocal[i][indn+kn];
        //ANTERPOLATE
        for(int km = 0; km < numsamp[i+1]; km++){
          complex<double> reduce = 0;
          for(int k = 0; k < 2*ninter; k++)
            reduce=reduce+anterp[i+1][km*2*ninter+k]*temp[mythread][antind[i+1][km*2*ninter+k]];
          agglocal[i+1][indm+km]=agglocal[i+1][indm+km]+reduce;
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    l2l_w[i] = l2l_w[i]+(get_wall_time()-walltemp);
  }
  double walltemp = get_wall_time();
  //LOWEST-LEVEL DISAGGREGATION
  #pragma omp parallel for
  for(int clusm = myid*numclus[level-1]/numproc; clusm < (myid+1)*numclus[level-1]/numproc; clusm++){
    int unk = clusm*box*box;
    int indlocal = clusm*numsamp[level-1];
    for(int n = 0; n < box*box; n++){
      complex<double> reduce = 0;
      int indbasis = n*numsamp[level-1];
      for(int k = 0; k < numsamp[level-1]; k++)
        reduce=reduce+basis_local[indbasis+k]*agglocal[level-1][indlocal+k];
      r[unk+n] = reduce;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  l2l_w[level-1] = l2l_w[level-1]+(get_wall_time()-walltemp);
  disah_w = get_wall_time()-disah_w;
  commn_w = get_wall_time();
  #pragma omp parallel for
  for(int cl = 0; cl < clcount[level]; cl++)
    memcpy(&sendbuff[level][cl*box*box],&x[sendmap[level][cl]*box*box],box*box*sizeof(complex<double>));
  //FARFIELD COMMUNICATION
  for(int p = 0; p < numproc; p++){
    int sib = procmap[level][p];
    if(sib != -1){
      complex<double> *send = &sendbuff[level][procint[level][p]*box*box];
      complex<double> *recv = &recvbuff[level][procint[level][p]*box*box];
      int amount = clusmap[level][p]*box*box;
      MPI_Sendrecv(send,amount,MPI_DOUBLE_COMPLEX,sib,0,recv,amount,MPI_DOUBLE_COMPLEX,sib,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
  }
  #pragma omp parallel for
  for(int cl = 0; cl < clcount[level]; cl++)
    memcpy(&x[recvmap[level][cl]*box*box],&recvbuff[level][cl*box*box],box*box*sizeof(complex<double>));
  MPI_Barrier(MPI_COMM_WORLD);
  commn_w = get_wall_time()-commn_w;
  nearf_w = get_wall_time();
  //NEARFIELD
  #pragma omp parallel for
  for(int clusm = myid*numclus[level-1]/numproc; clusm < (myid+1)*numclus[level-1]/numproc; clusm++){
    int testing = clusm*box*box;
    for(int m = 0; m < box*box; m++){
      complex<double> reduce = 0;
      for(int cn = 0; cn < 9; cn++){
        int clusn = clusnear[level-1][clusm*9+cn];
        if(clusn != -1){
          int ind = neid[clusm*9+cn];
          int basis = clusn*box*box;
          int indbox = ind*box*box*box*box+m*box*box;
          for(int n = 0; n < box*box; n++)
            reduce=reduce+x[basis+n]*near[indbox+n];
        }
      }
      r[testing+m]=r[testing+m]+reduce;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  nearf_w = get_wall_time()-nearf_w;
  matve_w = get_wall_time()-matve_w;
  matve_wall = matve_wall + matve_w;
  aggrh_wall = aggrh_wall + aggrh_w;
  trans_wall = trans_wall + trans_w;
  disah_wall = disah_wall + disah_w;
  nearf_wall = nearf_wall + nearf_w;
  commn_wall = commn_wall + commn_w;
  commf_wall = commf_wall + commf_w;
}
void mlfma_mpi(complex<double> *x, complex<double> *r){
  extern int numunk;
  extern complex<double> *rproc;
  extern complex<double> *xproc;
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0)printf("start\n");
  MPI_Scatter(x,numunk/numproc,MPI_DOUBLE_COMPLEX,&xproc[numunk/numproc*myid],numunk/numproc,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
  mlfma(xproc,rproc);
  MPI_Gather(&rproc[numunk/numproc*myid],numunk/numproc,MPI_DOUBLE_COMPLEX,r,numunk/numproc,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0)printf("end\n");
}

void nearfield(complex<double> *x, complex<double> *r){

  extern int box;
  extern int *numclus;
  extern complex<double> *near;
  extern int *neid;
  extern int **clusnear;
  extern int level;

  //NEARFIELD
  #pragma omp parallel for
  for(int clusm = 0; clusm < numclus[level-1]; clusm++){
    int testing = clusm*box*box;
    for(int m = 0; m < box*box; m++){
      complex<double> reduce = 0;
      for(int cn = 0; cn < 9; cn++){
        int clusn = clusnear[level-1][clusm*9+cn];
        if(clusn != -1){
          int ind = neid[clusm*9+cn];
          int basis = clusn*box*box;
          int indbox = ind*box*box*box*box+m*box*box;
          for(int n = 0; n < box*box; n++)
            reduce=reduce+x[basis+n]*near[indbox+n];
        }
      }
      r[testing+m]=r[testing+m]+reduce;
    }
  }
}
