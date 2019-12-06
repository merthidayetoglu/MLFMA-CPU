//Mert Hidayetoglu, May 2016
#include "vars.h"

extern double res;
extern int box;
extern int level;
extern int numunk;
extern int ninter;

extern int *numclus;
extern double *clusize;
extern int *numsamp;
extern int *numterm;

extern complex<double> *pos;
extern complex<double> **clusc;
extern int **clusnear;
extern int **clusfar;
extern int **cluspar;
extern int **cluschi;
extern int *unkpar;
extern int *unkmap;

extern complex<double> *coeff_multi;
extern complex<double> *basis_local;
extern double **interp;
extern double **anterp;
extern int **intind;
extern int **antind;
extern complex<double> **shiftmul;
extern complex<double> **shiftloc;
extern complex<double> **trans;
extern int **traid;
extern complex<double> *near;
extern int *neid;
extern complex<double> **aggmulti;
extern complex<double> **agglocal;
extern complex<double> **temp;

extern complex<double> *x;
extern complex<double> *b;
extern complex<double> *r;
extern complex<double> *o;
extern double *hist;
extern int max_it;

extern complex<double> *rproc;
extern complex<double> *xproc;
extern int **sendmap;
extern int **recvmap;
extern int **procmap;
extern int **clusmap;
extern int **procint;
extern int *clcount;
extern complex<double> **sendbuff;
extern complex<double> **recvbuff;

extern double *m2m_w;
extern double *m2l_w;
extern double *l2l_w;

extern double mem;

void preproc(){
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  if(myid==0)printf("ALLOCATION\n");
  //PREPROCESSING
  pos = new complex<double>[numunk];
  clusc = new complex<double>*[level];
  clusnear = new int*[level];
  clusfar = new int*[level];
  cluspar = new int*[level];
  cluschi = new int*[level];
  for(int i = 0; i < level; i++){
    clusc[i] = new complex<double>[numclus[i]];
    clusnear[i] = new int[numclus[i]*9];
    clusfar[i] = new int[numclus[i]*27];
    cluspar[i] = new int[numclus[i]];
    if(i < level-1){
      cluschi[i] = new int[numclus[i]*4];
      mem = mem + (double)sizeof(int)*4*numclus[i]/1024/1024;
    }
    else{
      cluschi[i] = new int[numclus[i]*box*box];
      mem = mem + (double)sizeof(int)*box*box*numclus[i]/1024/1024;
    }
    mem = mem + (double)sizeof(complex<double>)*numclus[i]/1024/1024;
    mem = mem + (double)sizeof(int)*37*numclus[i]/1024/1024;
  }
  unkpar = new int[numunk];
  unkmap = new int[numunk];
  mem = mem + (double)sizeof(complex<double>)*numunk/1024/1024;
  mem = mem + (double)sizeof(complex<double>*)*level/1024/1024;
  mem = mem + (double)sizeof(int*)*4*level/1024/1024;
  mem = mem + (double)sizeof(int)*2*numunk/1024/1024;
  if(myid==0){
  printf("NUMBER OF PROCSS.: %d\n",numproc);
  #pragma omp parallel
  if(omp_get_thread_num()==0)
    printf("NUMBER Of THREADS: %d\n",omp_get_num_threads());
  printf("TREE\n");
  }
  //TREE STRUCTURE
  cluspar[0][0]=-1;
  for(int i = 0; i < level; i++){
    int lenpar = pow(2,i);
    #pragma omp parallel for
    for(int m = 0; m < lenpar; m++)
      for(int n = 0; n < lenpar; n++){
        int par = m*lenpar+n;
        //HIGHER LEVELS
        if(i < level-1){
          for(int k = 0; k < 2; k++)
            for(int l = 0; l < 2; l++){
              //int chi = m*lenchi*2+n*2+lenchi*k+l;
              int chi = par*4+k*2+l;
              cluschi[i][par*4+k*2+l]=chi;
              cluspar[i+1][chi]=par;
            }
        }
        //LOWEST LEVEL
        else{
          for(int k = 0; k < box; k++)
            for(int l = 0; l < box; l++){
              //int chi = m*lenchi*box+n*box+lenchi*k+l;
              int chi = par*box*box+k*box+l;
              cluschi[i][par*box*box+k*box+l]=chi;
              unkpar[chi]=par;
            }
        }
      }
  }
  if(myid==0)printf("COORDINATES\n");
  //CLUSTER COORDINATES
  clusc[0][0]=0;
  for(int i = 0; i < level; i++){
    #pragma omp parallel for
    for(int m = 0; m < numclus[i]; m++){
      if(i < level-1){
        complex<double>cor=clusc[i][m]+clusize[i]/4*complex<double>(-1,+1);
        for(int k = 0; k < 2; k++)
          for(int l = 0; l < 2; l++)
            clusc[i+1][cluschi[i][m*4+k*2+l]]=cor+complex<double>(l*clusize[i+1],-k*clusize[i+1]);
      }
      //LOWEST LEVEL
      else{
        complex<double>cor=clusc[i][m]+(clusize[i]/2-res/2)*complex<double>(-1,+1);
        for(int k = 0; k < box; k++)
          for(int l = 0; l < box; l++){
            pos[cluschi[i][m*box*box+k*box+l]]=cor+complex<double>(l*res,-k*res);
          }
      }
    }
  }
  if(myid==0)printf("MAPPING\n");
  //UNKNOWN MAPPING
  complex<double> cor = clusize[0]/2*complex<double>(-1,1);
  int len = box*pow(2,level-1);
  #pragma omp parallel for
  for(int n = 0; n < numunk; n++){
    int l = floor((pos[n].real()-cor.real())/res);
    int k = floor((cor.imag()-pos[n].imag())/res);
    unkmap[k*len+l]=n;
  }
  if(myid==0)printf("NEAR-FAR\n");
  //NEAR-FIELD & FAR-FIELD CLUSTERS
  for(int i = 0; i < level; i++){
    #pragma omp parallel for
    for(int m = 0; m < numclus[i]; m++){
      for(int n = 0; n < 9; n++)
        clusnear[i][m*9+n]=-1;
      for(int n = 0; n < 27; n++)
        clusfar[i][m*27+n]=-1;
    }
  }
  clusnear[0][0]=0;
  for(int i = 1; i < level; i++){
    double clusdim = clusize[i];
    #pragma omp parallel for
    for(int m = 0; m < numclus[i]; m++){
      int nearid = 0;
      int farid = 0;
      complex<double>posm=clusc[i][m];
      int par = cluspar[i][m];
      for(int ind = 0; ind < 9; ind++){
        int uncle = clusnear[i-1][par*9+ind];
        if(uncle != -1){
          for(int ind2 = 0; ind2 < 4; ind2++){
            int n=cluschi[i-1][uncle*4+ind2];
            if(abs(clusc[i][n]-posm) < clusdim*sqrt(3)){
              clusnear[i][m*9+nearid]=n;
              nearid++;
            }
            else{
              clusfar[i][m*27+farid]=n;
              farid++;
            }
          }
        }
      }
    }
  }
}
void setup_mpi(){
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  rproc = new complex<double>[numunk];
  xproc = new complex<double>[numunk];

  int *proctemp = new int[numproc];
  int **friends = new int*[level+1];
  int **clustemp = new int*[level+1];

  for(int i = 0; i < level; i++){
    friends[i] = new int[numproc];
    clustemp[i] = new int[numclus[i]];
    fill_n(friends[i],numproc,0);
    fill_n(clustemp[i],numclus[i],-1);
  }
  friends[level] = new int[numproc];
  clustemp[level] = new int[numclus[level-1]];
  fill_n(friends[level],numproc,0);
  fill_n(clustemp[level],numclus[level-1],-1);
  //CLUSTEMP & FRIENDS
  for(int i = 0; i < level; i++){
    for(int clusm = 0; clusm < numclus[i]/numproc; clusm++)
      for(int cn = 0; cn < 27; cn++){
        int clusn = clusfar[i][(numclus[i]/numproc*myid+clusm)*27+cn];
        if(clusn!=-1)clustemp[i][clusn] = clusn/(numclus[i]/numproc);
      }
    for(int m = 0; m < numclus[i]; m++){
      int p = clustemp[i][m];
      if(p!=-1)friends[i][p]++;
    }
  }
  for(int clusm = 0; clusm < numclus[level-1]/numproc; clusm++)
    for(int cn = 0; cn < 9; cn++){
      int clusn = clusnear[level-1][(numclus[level-1]/numproc*myid+clusm)*9+cn];
      if(clusn!=-1)clustemp[level][clusn] = clusn/(numclus[level-1]/numproc);
    }
  for(int m = 0; m < numclus[level-1]; m++){
    int p = clustemp[level][m];
    if(p!=-1)friends[level][p]++;
  }
  sendmap = new int*[level+1];
  recvmap = new int*[level+1];
  procint = new int*[level+1];
  procmap = new int*[level+1];
  sendbuff = new complex<double>*[level+1];
  recvbuff = new complex<double>*[level+1];
  clusmap = new int*[level+1];
  clcount = new int[level+1];
  for(int i = 0; i < level+1; i++){
    procint[i] = new int[numproc];
    procmap[i] = new int[numproc];
    clusmap[i] = new int[numproc];
    fill_n(procint[i],numproc,-1);
    fill_n(procmap[i],numproc,-1);
    fill_n(clusmap[i],numproc,0);
  }
  //PROCTEMP
  {
    int count;
    int me,sib;
    char strtemp[80];
    if(numproc>9)
      sprintf(strtemp,"internodes0%d.dat",numproc);
    else
      sprintf(strtemp,"internodes00%d.dat",numproc);
    if(myid==0)printf("FILE %s\n",strtemp);
    FILE *mapf = fopen(strtemp,"r");
    fscanf(mapf,"%d",&count);
    count = 0;
    proctemp[count]=myid;
    for(int i = 1; i < numproc; i++)
      for(int n = 0; n < numproc; n++){
        fscanf(mapf,"%d %d",&me,&sib);
        if(me == myid){
          count++;
          proctemp[count] = sib;
        }
      }
    fclose(mapf);
  }
  //PROCMAP
  if(myid==0)printf("PROCMAP\n");
  for(int i = 0; i < level+1; i++){
    int numsib = 0;
    int numcls = 0;
    for(int p = 0; p < numproc; p++){
      if(friends[i][proctemp[p]] > 0 && proctemp[p]!=myid){
        procmap[i][numsib] = proctemp[p];
        clusmap[i][numsib] = friends[i][proctemp[p]];
        numcls = numcls + friends[i][proctemp[p]];
        numsib++;
      }
    }
    clcount[i] = numcls;
  }
  if(myid==0){
    printf("PROCMAP %d\n",myid);
    for(int i = 0; i < level+1; i++){
      for(int p = 0; p < numproc; p++)
        printf("%d ",procmap[i][p]);
      printf("\n");
    }
    printf("CLUSMAP %d\n",myid);
    for(int i = 0; i < level; i++){
      for(int p = 0; p < numproc; p++)
        printf("%d ",clusmap[i][p]);
      printf("| %d | %d\n",clcount[i],clcount[i]*numsamp[i]);
    }
    for(int p = 0; p < numproc; p++)
      printf("%d ",clusmap[level][p]);
    printf("| %d | %d\n",clcount[level],clcount[level]*box*box);
  }
  for(int i = 2; i < level+1; i++){
    recvmap[i] = new int[clcount[i]];
    sendmap[i] = new int[clcount[i]];
    if(i < level){
      sendbuff[i] = new complex<double>[clcount[i]*numsamp[i]];
      recvbuff[i] = new complex<double>[clcount[i]*numsamp[i]];
    }else{
      sendbuff[i] = new complex<double>[clcount[i]*box*box];
      recvbuff[i] = new complex<double>[clcount[i]*box*box];
    }
  }
  //RECVMAP
  if(myid==0)printf("RECVMAP\n");
  for(int i = 0; i < level+1; i++){
    int count = 0;
    for(int p = 0; p < numproc; p++){
      int procm = procmap[i][p];
      if(procm != -1){
        procint[i][p] = count;
        if(i < level)
          for(int m = 0; m < numclus[i]/numproc; m++){
            int clusn = clustemp[i][numclus[i]/numproc*procm+m];
            if(clusn != -1){
              recvmap[i][count] = procm*numclus[i]/numproc+m;
              count++;
            }
          }
        else
          for(int m = 0; m < numclus[level-1]/numproc; m++){
            int clusn = clustemp[i][numclus[level-1]/numproc*procm+m];
            if(clusn != -1){
              recvmap[i][count] = procm*numclus[level-1]/numproc+m;
              count++;
            }
          }
      }
    }
  }
  if(myid==0){
    printf("PROCINT %d\n",myid);
    for(int i = 0; i < level+1; i++){
      for(int p = 0; p < numproc; p++)
        printf("%d ",procint[i][p]);
      printf("\n");
    }
    printf("RECVMAP %d\n",myid);
    for(int i = 0; i < 5; i++){
      for(int m = 0; m < clcount[i]; m++)
        printf("%d ",recvmap[i][m]);
      printf(" **********\n");
    }
  }
  //SENDMAP
  for(int i = 2; i < level+1; i++){
    if(myid==0)printf("SENDMAP  %d\n",i+1);
    for(int p = 0; p < numproc; p++){
      int sib = procmap[i][p];
      if(sib != -1){
        int *send = &recvmap[i][procint[i][p]];
        int *recv = &sendmap[i][procint[i][p]];
        int amount = clusmap[i][p];
        MPI_Sendrecv(send,amount,MPI_INTEGER,sib,0,recv,amount,MPI_INTEGER,sib,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    }
  }
  if(myid==0){
    printf("SENDMAP %d\n",myid);
    for(int i = 0; i < 5; i++){
      for(int m = 0; m < clcount[i]; m++)
        printf("%d ",sendmap[i][m]);
      printf(" **********\n");
    }
  }
  for(int i = 0; i < level; i++){
    delete[] friends[i];
    delete[] clustemp[i];
  }
  delete[] proctemp;
  delete[] friends;
  delete[] clustemp;
}
void setup(){
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  if(myid==0)printf("ALLOCATION\n");
  //SETUP
  complex<double> *hank = new complex<double>[numterm[0]+1];
  complex<double> *expn = new complex<double>[numterm[0]];
  int *list = new int[numsamp[0]];
  near = new complex<double>[9*box*box*box*box];
  neid = new int[numclus[level-1]*9];
  coeff_multi = new complex<double>[numsamp[level-1]*box*box];
  basis_local = new complex<double>[box*box*numsamp[level-1]];
  mem = mem + (double)sizeof(complex<double>)*(numterm[0]+1);
  mem = mem + (double)sizeof(complex<double>)*(numterm[0]);
  mem = mem + (double)sizeof(int)*numsamp[0];
  mem = mem + (double)sizeof(complex<double>)*9*box*box*box*box/1024/1024;
  mem = mem + (double)sizeof(int)*numclus[level-1]*9/1024/1024;
  mem = mem + (double)sizeof(complex<double>)*numsamp[level-1]*box*box/1024/1024;
  mem = mem + (double)sizeof(complex<double>)*box*box*numsamp[level-1]/1024/1024;

  mem = mem + (double)sizeof(double*)*2*level/1024/1024;
  mem = mem + (double)sizeof(int*)*3*level/1024/1024;
  mem = mem + (double)sizeof(complex<double>*)*3*level/1024/1024;
  interp = new double*[level];
  anterp = new double*[level];
  intind = new int*[level];
  antind = new int*[level];
  shiftmul = new complex<double>*[level];
  shiftloc = new complex<double>*[level];
  trans = new complex<double>*[level];
  traid = new int*[level];
  for(int i = 0; i < level; i++){
    interp[i] = new double[numsamp[i]*ninter];
    anterp[i] = new double[numsamp[i]*2*ninter];
    intind[i] = new int[numsamp[i]*ninter];
    antind[i] = new int[numsamp[i]*2*ninter];
    shiftmul[i] = new complex<double>[4*numsamp[i]];
    shiftloc[i] = new complex<double>[4*numsamp[i]];
    trans[i] = new complex<double>[49*numsamp[i]];
    traid[i] = new int[numclus[i]*27];
    mem = mem + (double)sizeof(double)*2*numsamp[i]*ninter/1024/1024;
    mem = mem + (double)sizeof(int)*2*numsamp[i]*ninter/1024/1024;
    mem = mem + (double)sizeof(complex<double>)*57*numsamp[i]/1024/1024;
    mem = mem + (double)sizeof(int)*27*numclus[i]/1024/1024;
  }
  if(myid==0)printf("SETUP MEMORY: %f MB\n",mem/1024.0);

  //FILL NEARFIELD MATRIX
  if(myid==0)printf("NEARFIELD MATRIX\n");
  #pragma omp parallel for
  for(int n = 0; n < box*box*box*box; n++)
    for(int m = 0; m < 9; m++)
      near[n*9+m]=0;
  #pragma omp parallel for
  for(int clusm = 0; clusm < numclus[level-1]; clusm++)
    for(int cn = 0; cn < 9; cn++){
      int clusn = clusnear[level-1][clusm*9+cn];
      if(clusn != -1){
        complex<double> m2l = clusc[level-1][clusm]-clusc[level-1][clusn];
        int m = round(m2l.imag()/clusize[level-1])+1;
        int n = round(-m2l.real()/clusize[level-1])+1;
        int ind = m*3+n;
        neid[clusm*9+cn]=ind;
      }
      else
        neid[clusm*9+cn]=-1;
    }
  #pragma omp parallel for
  for(int clusm = 0; clusm < numclus[level-1]; clusm++)
    for(int cn = 0; cn < 9; cn++){
      int clusn = clusnear[level-1][clusm*9+cn];
      int ind = neid[clusm*9+cn];
      if(ind != -1)
        if(near[ind*box*box*box*box]==complex<double>(0))
          for(int m = 0; m < box*box*box*box; m++){
            int k = m/(box*box);
            int l = m%(box*box);
            int testing = cluschi[level-1][clusm*box*box+k];
            int basis = cluschi[level-1][clusn*box*box+l];
            near[ind*box*box*box*box+k*box*box+l]=integrate(pos[testing],pos[basis]);
          }
    }
  if(myid==0)printf("TRANSLATION OPERATORS\n");
  #pragma omp parallel for
  for(int p = 0; p < numterm[0]; p++)
    expn[p] = exp(complex<double>(0,(p+1)*M_PI/2));
  //FILL TRANSLATION OPERATORS
  for(int i = 2; i < level; i++){
    if(myid==0)printf("level %d\n",i+1);
    double a = clusize[i];
    #pragma omp parallel for
    for(int clusm = 0; clusm < numclus[i]; clusm++)
      for(int cn = 0; cn < 27; cn++){
        int clusn = clusfar[i][clusm*27+cn];
        if(clusn != -1){
          complex<double> l2m = clusc[i][clusn]-clusc[i][clusm];
          int n = round(l2m.real()/a)+3;
          int m = round(-l2m.imag()/a)+3;
          int ind = m*7+n;
          traid[i][clusm*27+cn]=ind;
        }
      }
    for(int cm = 0; cm < 7; cm++)
      for(int cn = 0; cn < 7; cn++){
        complex<double> l2m = -complex<double>(-3*a+cn*a,3*a-cm*a);
        if(abs(l2m)>a*sqrt(3)){
          #pragma omp parallel for
          for(int p = 0; p < numterm[i]+1; p++)
            hank[p] = hn(p,2*M_PI*abs(l2m));
          int index = cm*7*numsamp[i]+cn*numsamp[i];
          #pragma omp parallel for
          for(int k = 0; k < numsamp[i]; k++){
            complex<double> reduce = hank[0];
            for(int p = 1; p < numterm[i]+1; p++)
              reduce=reduce+2*cos(p*(2*M_PI*k/numsamp[i]-arg(l2m)))*hank[p]*expn[p-1];
            trans[i][index+k] = reduce;
          }
        }
      }
  }
  //FILL COEFF FOR MULTIPOLE
  //FILL BASIS FOR LOCAL
  #pragma omp parallel for
  for(int n = 0; n < box*box; n++){
    int unk = cluschi[level-1][n];
    for(int k = 0; k < numsamp[level-1]; k++){
      coeff_multi[k*box*box+n] = integrate_multi(clusc[level-1][0],pos[unk],k);
      basis_local[n*numsamp[level-1]+k]=integrate_local(clusc[level-1][0],pos[unk],k)/complex<double>(numsamp[level-1])*complex<double>(0,0.25);
    }
  }
  //FILL SHIFTERS
  if(myid==0)printf("MULTIPOLE & LOCAL SHIFTERS\n");
  for(int i = 0; i < level-1; i++){
    int par = 0;
    for(int cn = 0; cn < 4; cn++){
      int chi = cluschi[i][par*numclus[i]+cn];
      complex<double> rho = clusc[i][par]-clusc[i+1][chi];
      #pragma omp parallel for
      for(int k = 0; k < numsamp[i]; k++){
        double ang = 2*M_PI*k/numsamp[i];
        complex<double> kw = 2*M_PI*complex<double>(cos(ang),sin(ang));
        double product = kw.real()*rho.real()+kw.imag()*rho.imag();
        shiftmul[i][cn*numsamp[i]+k]=exp(complex<double>(0,product));
        shiftloc[i][cn*numsamp[i]+k]=exp(complex<double>(0,-product));
      }
    }
  }
  if(myid==0)printf("INTERPOLATORS\n");
  //FILL INTERPOLATORS
  for(int i = 0; i < level-1; i++){
    double ratio = (double)numsamp[i+1]/numsamp[i];
    #pragma omp parallel for
    for(int m = 0; m < numsamp[i]; m++){
      int center = 0;
      if(ninter%2 == 0)
        center = ceil(m*ratio);
      else
        center = round(m*ratio);
      double xm = (double)m/numsamp[i];
      for(int n = 0; n < ninter; n++){
        int ind = center+n-ninter/2;
        double mul = 1;
        double xn = (double)ind/numsamp[i+1];
        for(int k = 0; k < ninter; k++)
          if(k != n){
            double xk = (double)(center+k-ninter/2)/numsamp[i+1];
            mul = mul*(xm-xk)/(xn-xk);
          }
        interp[i][m*ninter+n]=mul;
        ind = ind%numsamp[i+1];
        if(ind < 0)
          ind = ind + numsamp[i+1];
        intind[i][m*ninter+n] = ind;
      }
    }
  }
  if(myid==0)printf("ANTERPOLATORS\n");
  //FILL ANTERPOLATORS
  for(int i = 3; i < level; i++){
    for(int n = 0; n < numsamp[0]; n++)
      list[n] = 0;
    for(int m = 0; m < numsamp[i]; m++)
      for(int n = 0; n < 2*ninter; n++){
        antind[i][m*2*ninter+n] = 0;
        anterp[i][m*2*ninter+n] = 0;
      }
    for(int m = 0; m < numsamp[i-1]; m++)
      for(int n = 0; n < ninter; n++){
        int ind = intind[i-1][m*ninter+n];
        antind[i][ind*2*ninter+list[ind]] = m;
        anterp[i][ind*2*ninter+list[ind]] = interp[i-1][m*ninter+n]*numsamp[i]/numsamp[i-1];
        list[ind]++;
      }
  }

  //FILE *interf = fopen("interp.txt","w");
  //FILE *intinf = fopen("intind.txt","w");
  //for(int i = 0; i < numsamp[level-2]; i++){
  //  for(int j = 0; j < ninter; j++){
  //    fprintf(interf,"%e ",interp[level-2][i*ninter+j]);
  //    fprintf(intinf,"%d ",intind[level-2][i*ninter+j]);
  //  }
  //  fprintf(interf,"\n");
  //  fprintf(intinf,"\n");
  //}
  //fclose(interf);
  //fclose(intinf);

  delete[] hank;
  delete[] list;
  delete[] expn;
  mem = mem - (double)sizeof(complex<double>)*(numterm[0]+1);
  mem = mem - (double)sizeof(complex<double>)*(numterm[0]);
  mem = mem - (double)sizeof(int)*numsamp[0];

  o = new complex<double>[numunk];
  b = new complex<double>[numunk];
  x = new complex<double>[numunk];
  r = new complex<double>[numunk];
  hist = new double[max_it+1];
  mem = mem + (double)sizeof(complex<double>)*4*numunk/1024/1024;
  mem = mem + (double)sizeof(double)*(max_it+1)/1024/1024;
  for(int i = 0; i < max_it+1; i++)
    hist[i] = 0;

  aggmulti = new complex<double>*[level];
  agglocal = new complex<double>*[level];
  mem = mem + (double)sizeof(complex<double>*)*2*level/1024/1024;
  for(int i = 0; i < level; i++){
    aggmulti[i] = new complex<double>[numclus[i]*numsamp[i]];
    agglocal[i] = new complex<double>[numclus[i]*numsamp[i]];
    mem = mem + (double)sizeof(complex<double>)*numclus[i]*numsamp[i]*2/1024/1024;
  }
  #pragma omp parallel shared(temp,mem)
  {
    if(omp_get_thread_num() == 0){
      mem = mem + (double)sizeof(complex<double>*)*omp_get_num_threads()/1024/1024;
      temp = new complex<double>*[omp_get_num_threads()];
      for(int i = 0; i < omp_get_num_threads(); i++){
        temp[i] = new complex<double>[numsamp[0]];
         mem = mem + (double)sizeof(complex<double>)*numsamp[0]/1024/1024;
       }
    }
  }
  m2m_w = new double[level];
  m2l_w = new double[level];
  l2l_w = new double[level];
  for(int i = 0; i < level; i++){
    m2m_w[i] = 0;
    m2l_w[i] = 0;
    l2l_w[i] = 0;
  }
  mem = mem + (double)sizeof(double)*3*level/1024/1024;
}

void final(){
  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  //SAVE FILES
  if(myid==0){
  #pragma omp parallel for
  for(int n = 0; n < numunk; n++)
    r[n] = b[unkmap[n]];
  FILE *incb = fopen("inc.bin","wb");
  fwrite(r,sizeof(complex<double>),numunk,incb);
  fclose(incb);
  #pragma omp parallel for
  for(int n = 0; n < numunk; n++)
    r[n] = x[unkmap[n]]-b[unkmap[n]];
  FILE *scab = fopen("sca.bin","wb");
  fwrite(r,sizeof(complex<double>),numunk,scab);
  fclose(scab);
  #pragma omp parallel for
  for(int n = 0; n < numunk; n++)
    r[n] = x[unkmap[n]];
  FILE *totb = fopen("tot.bin","wb");
  fwrite(r,sizeof(complex<double>),numunk,totb);
  fclose(totb);/*
  for(int n = 0; n < numunk; n++)
    r[n] = o[unkmap[n]]/(4*M_PI*M_PI)+complex<double>(1);
  FILE *epsb = fopen("eps.bin","wb");
  fwrite(r,sizeof(complex<double>),numunk,epsb);
  fclose(epsb);
  for(int n = 0; n < numunk; n++)
    r[n] = pos[unkmap[n]];
  FILE *posb = fopen("pos.bin","wb");
  fwrite(r,sizeof(complex<double>),numunk,posb);
  fclose(posb);*/
  FILE *histf = fopen("hist.txt","w");
  for(int i = 0; i < max_it+1; i++)
    fprintf(histf,"%e\n",hist[i]);
  fclose(histf);
  printf("FILES DONE... EXIT!\n");
  }

  for(int i = 0; i < level; i++){
    delete[] clusc[i];
    delete[] cluspar[i];
    delete[] cluschi[i];
    if(i < level-1){
      mem = mem - (double)sizeof(int)*4*numclus[i]/1024/1024;
    }
    else{
      mem = mem - (double)sizeof(int)*box*box*numclus[i]/1024/1024;
    }
    mem = mem - (double)sizeof(complex<double>)*numclus[i]/1024/1024;
    mem = mem - (double)sizeof(int)*37*numclus[i]/1024/1024;
  }
  delete[] pos;
  delete[] clusc;
  delete[] clusnear;
  delete[] clusfar;
  delete[] cluspar;
  delete[] cluschi;
  delete[] unkpar;
  delete[] unkmap;
  mem = mem - (double)sizeof(complex<double>)*numunk/1024/1024;
  mem = mem - (double)sizeof(complex<double>*)*level/1024/1024;
  mem = mem - (double)sizeof(int*)*level*4/1024/1024;
  mem = mem - (double)sizeof(int)*2*numunk/1024/1024;

delete[] near;
  delete[] neid;
  delete[] coeff_multi;
  delete[] basis_local;
  mem = mem - (double)sizeof(complex<double>)*9*box*box*box*box/1024/1024;
  mem = mem - (double)sizeof(int)*numclus[level-1]*9/1024/1024;
  mem = mem - (double)sizeof(complex<double>)*numsamp[level-1]*box*box/1024/1024;
  mem = mem - (double)sizeof(complex<double>)*box*box*numsamp[level-1]/1024/1024;
  for(int i = 0; i < level; i++){
    delete[] interp[i];
    delete[] anterp[i];
    delete[] intind[i];
    delete[] antind[i];
    delete[] shiftmul[i];
    delete[] shiftloc[i];
    delete[] trans[i];
    delete[] traid[i];
    mem = mem - (double)sizeof(double)*2*numsamp[i]*ninter/1024/1024;
    mem = mem - (double)sizeof(int)*2*numsamp[i]*ninter/1024/1024;
    mem = mem - (double)sizeof(complex<double>)*57*numsamp[i]/1024/1024;
    mem = mem - (double)sizeof(int)*27*numclus[i]/1024/1024;
  }
  delete[] interp;
  delete[] anterp;
  delete[] intind;
  delete[] antind;
  delete[] shiftmul;
  delete[] shiftloc;
  delete[] trans;
  delete[] traid;
  mem = mem - (double)sizeof(double*)*2*level/1024/1024;
  mem = mem - (double)sizeof(int*)*3*level/1024/1024;
  mem = mem - (double)sizeof(complex<double>*)*3*level/1024/1024;
  for(int i = 0; i < level; i++){
    delete[] aggmulti[i];
    delete[] agglocal[i];
    mem = mem - (double)sizeof(complex<double>)*2*numclus[i]*numsamp[i]/1024/1024;
  }
  delete[] aggmulti;
  delete[] agglocal;
  mem = mem - (double)sizeof(complex<double>*)*2*level/1024/1024;
  #pragma omp parallel shared(mem)
  {
    if(omp_get_thread_num() == 0){
      for(int i = 0; i < omp_get_num_threads(); i++){
        delete[] temp[i];
        mem = mem - (double)sizeof(complex<double>)*numsamp[0]/1024/1024;
      }
      delete[] temp;
      mem = mem - (double)sizeof(complex<double>*)*omp_get_num_threads()/1024/1024;
    }
  }
  delete[] o;
  delete[] b;
  delete[] x;
  delete[] r;
  delete[] hist;
  mem = mem - (double)sizeof(complex<double>)*4*numunk/1024/1024;
  mem = mem - (double)sizeof(double)*(max_it+1)/1024/1024;
  delete[] numclus;
  delete[] clusize;
  delete[] numsamp;
  delete[] numterm;
  mem = mem - (double)sizeof(int)*3*level/1024/1024;
  mem = mem - (double)sizeof(double)*level/1024/1024;
}
