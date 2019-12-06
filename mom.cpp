//Mert Hidayetoglu, Sep 2015
#include "vars.h"

  //input data
  double res = 0.1;
  int box = 8;
  int level = 10;

  //solver data
  int beta = 5;
  int ninter = 30;
  int max_it = 300;
  double res_tol = pow(10,-4);

  double dim;
  int numunk;

  int *numclus;
  double *clusize;
  int *numsamp;
  int *numterm;

  complex<double> *pos;
  complex<double> **clusc;
  int **clusnear;
  int **clusfar;
  int **cluspar;
  int **cluschi;
  int *unkpar;
  int *unkmap;

  complex<double> *near;
  int *neid;
  complex<double> *coeff_multi;
  complex<double> *basis_local;
  complex<double> *o;
  complex<double> *b;
  complex<double> *x;
  complex<double> *r;
  double *hist;

  double **interp;
  double **anterp;
  int **intind;
  int **antind;
  complex<double> **shiftmul;
  complex<double> **shiftloc;
  complex<double> **trans;
  int **traid;
  complex<double> **aggmulti;
  complex<double> **agglocal;
  complex<double> **temp;

  complex<double> *rproc;
  complex<double> *xproc;
  int **sendmap;
  int **recvmap;
  int **procmap;
  int **clusmap;
  int **procint;
  int *clcount;
  complex<double> **sendbuff;
  complex<double> **recvbuff;

  int nummatvec = 0;
  double mem = 0;
  double gpumem = 0;
  double matve_wall = 0;
  double train_wall = 0;
  double aggrh_wall = 0;
  double *m2m_w;
  double trans_wall = 0;
  double *m2l_w;
  double disah_wall = 0;
  double *l2l_w;
  double nearf_wall = 0;
  double traut_wall = 0;
  double commn_wall = 0;
  double commf_wall = 0;

  double total_wall = get_wall_time();
  double prepr_wall = 0;
  double setup_wall = 0;
  double solut_wall = 0;
  double postp_wall = 0;

int main(int argc, char **argv){

  int numproc;
  int myid;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  numclus = new int[level];
  clusize = new double[level];
  numsamp = new int[level];
  numterm = new int[level];
  mem = mem + (double)sizeof(int)*3*level/1024/1024;
  mem = mem + (double)sizeof(double)*level/1024/1024;

  for(int i = 0; i < level; i++)
    numclus[i] = pow(4,i);
  for(int i = 0; i < level; i++)
    clusize[i] = res*box*pow(2,level-i-1);
  dim = clusize[0];
  numunk = pow(dim/res,2);
  for(int i = 0; i < level; i++){
    double kr = sqrt(2)*2*M_PI*clusize[i];
    double term = kr+1.8*pow(beta,2.0/3)*pow(kr,1.0/3);
    numterm[i] = floor(term);
    numsamp[i] = 4*numterm[i];
  }

  if(myid==0){
  printf("PIXEL SIZE      : %f\n",res);
  printf("PIXEL-SIZE      : %d\n",box);
  printf("LOW LEVEL SIZE  : %f\n",res*box);
  printf("NUM LEVEL       : %d\n",level);
  printf("DIMENSION       : %f\n",dim);
  printf("NUM UNKNOWN     : %d\n",numunk);
  printf("NUM INTER       : %d\n",ninter);
  printf("NUM CLUS        :\n");
  for(int i = 0; i < level; i++)
    printf("LEVEL %d: %d\n",i+1,numclus[i]);

  printf("SOLVER DATA\n");
  for(int i = 0; i < level; i++)
    printf("LEVEL: %d  CLUSSIZE: %f  NUM. SAMP: %d  TOTAL: %d\n",i+1,clusize[i],numsamp[i],numsamp[i]*numclus[i]);
  printf("NUM TERMS       : \n");
  for(int i = 0; i < level; i++)
    printf("LEVEL: %d  NUM.TERM: %d\n",i+1,numterm[i]);

  printf("INTEGER: %d, DOUBLE: %d, POINTER: %d\n",(int)sizeof(int),(int)sizeof(double),(int)sizeof(complex<double>*));
  }

  //PREPROCESSING
  prepr_wall = get_wall_time();
  preproc();
  setup_mpi();
  prepr_wall = get_wall_time()-prepr_wall;
  if(myid==0)printf("PRE-PROCESSING MEMORY: %f MB\n",mem);

  setup_wall = get_wall_time();
  setup();
  setup_wall = get_wall_time()-setup_wall;
  if(myid==0)printf("SETUP MEMORY: %f MB\n",mem);

  complex<double> epsr(0.02,0);
  fill_n(o,numunk,complex<double>(0,0));
  #pragma omp parallel for
  for(int n = 0; n < numunk; n++)
    if(abs(pos[n])<100)
      o[n]=4*M_PI*M_PI*epsr;
  //CALCULATE RHS
  fill_n(x,numunk,complex<double>(0,0));
  #pragma omp parallel for
  for(int m = 0; m < numunk; m++)
    //b[m]=integrate(pos[m],complex<double>(-3,0))/res/res;
    b[m]=exp(complex<double>(0,2*M_PI*pos[m].real()));
  //SOLUTION
  MPI_Barrier(MPI_COMM_WORLD);
  solut_wall = get_wall_time();
  if(myid==0)printf("ITERATIVE SOLUTION\n");
  bicgs(x,o,b);
  //born(x,o,b);
  //conjg(x,o,b);

  //x[numunk/numproc*4-1] = 1;
  //for(int n = 0; n < 10; n++)
  //mlfma_mpi(x,b);
  //MPI_Gather(&x[numunk/numproc*myid],numunk/numproc,MPI_DOUBLE_COMPLEX,r,numunk/numproc,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
  //memcpy(x,r,numunk*sizeof(complex<double>));

  //if(myid==0){
  //printf("FINISHED\n");
  //for(int n = 0; n < numunk; n++)
  //  r[n] = x[n]*o[n];
  //double fartime = get_wall_time();
  //farfield(r);
  //printf("farfield time: %e\n",get_wall_time()-fartime);
  //}

  /*
  if(myid==numproc-1){
  FILE *testf = fopen("agg_cpu.txt","w");
  int addres = myid*numclus[2]/numproc*numsamp[2];
  for(int n = 0; n < numclus[2]/numproc*numsamp[2]; n++)
    fprintf(testf,"%e %e\n",aggmulti[2][addres+n].real(),aggmulti[2][addres+n].imag());
  fclose(testf);
  testf = fopen("loc_cpu.txt","w");
  for(int n = 0; n < numclus[2]/numproc*numsamp[2]; n++)
    fprintf(testf,"%e %e\n",agglocal[2][addres+n].real(),agglocal[2][addres+n].imag());
  fclose(testf);
  addres = myid*numclus[level-1]/numproc*numsamp[level-1];
  testf = fopen("agg_cpu_loww.txt","w");
  for(int n = 0; n < numclus[level-1]/numproc*numsamp[level-1]; n++)
    fprintf(testf,"%e %e\n",aggmulti[level-1][addres+n].real(),aggmulti[level-1][addres+n].imag());
  fclose(testf);
  testf = fopen("loc_cpu_loww.txt","w");
  for(int n = 0; n < numclus[level-1]/numproc*numsamp[level-1]; n++)
   fprintf(testf,"%e %e\n",agglocal[level-1][addres+n].real(),agglocal[level-1][addres+n].imag());
  fclose(testf);
  addres = myid*numclus[3]/numproc*numsamp[3];
  testf = fopen("loc_cpu_low.txt","w");
  for(int n = 0; n < numclus[3]/numproc*numsamp[3]; n++)
    fprintf(testf,"%e %e\n",agglocal[3][addres+n].real(),agglocal[3][addres+n].imag());
  fclose(testf);
  }*/

  MPI_Barrier(MPI_COMM_WORLD);
  solut_wall = get_wall_time()-solut_wall;
  postp_wall = get_wall_time();

  //if(myid==0){
  //for(int n = 0; n < numunk; n++)
  //  r[n] = integrate(pos[place],pos[unkmap[n]])/(res*res);
  //FILE *refb = fopen("ref.bin","wb");
  //fwrite(r,sizeof(complex<double>),numunk,refb);
  //fclose(refb);
  //}

  final();
  if(myid==0)printf("FILES DONE... EXIT!\n");
  postp_wall = get_wall_time()-postp_wall;

  total_wall = get_wall_time()-total_wall;
  //PRINT TIME
  if(myid==0){
  printf("\n");
  printf("MAVEC TOTAL TIME\n");
  printf("TRANSFER       : %e\n",train_wall);
  printf("AGGREGATION    : %e\n",aggrh_wall);
  for(int i = 0; i < level; i++)
  printf("AGGREGATION %d : %e\n",i+1,m2m_w[i]);
  printf("FAR COMM       : %e\n",commf_wall);
  printf("TRANSLATION    : %e\n",trans_wall);
  for(int i = 0; i < level; i++)
  printf("TRANSLATION %d : %e\n",i+1,m2l_w[i]);
  printf("DISAGGREG      : %e\n",disah_wall);
  for(int i = 0; i < level; i++)
  printf("DISAGGREG   %d : %e\n",i+1,l2l_w[i]);
  printf("FARFIELD       : %e\n",aggrh_wall+commf_wall+trans_wall+disah_wall);
  printf("NEAR COMM      : %e\n",commn_wall);
  printf("NEARFIELD      : %e\n",nearf_wall);
  printf("TRANSFER       : %e\n",traut_wall);
  double other_wall = matve_wall-aggrh_wall-trans_wall-disah_wall-nearf_wall-commf_wall-commn_wall;
  printf("OTHERS         : %e\n",other_wall);
  printf("TOTAL          : %e\n",matve_wall);
  printf("\n");
  printf("MAVEC SINGLE TIME\n");
  printf("TRANSFER       : %e\n",train_wall/nummatvec);
  printf("AGGREGATION    : %e\n",aggrh_wall/nummatvec);
  for(int i = 0; i < level; i++)
  printf("AGGREGATION %d : %e\n",i+1,m2m_w[i]/nummatvec);
  printf("FAR COMM       : %e\n",commf_wall/nummatvec);
  printf("TRANSLATION    : %e\n",trans_wall/nummatvec);
  for(int i = 0; i < level; i++)
  printf("TRANSLATION %d : %e\n",i+1,m2l_w[i]/nummatvec);
  printf("DISAGGREG      : %e\n",disah_wall/nummatvec);
  for(int i = 0; i < level; i++)
  printf("DISAGGREG   %d : %e\n",i+1,l2l_w[i]/nummatvec);
  printf("FARFIELD       : %e\n",(aggrh_wall+commf_wall+trans_wall+disah_wall)/nummatvec);
  printf("NEAR COMM      : %e\n",commn_wall/nummatvec);
  printf("NEARFIELD      : %e\n",nearf_wall/nummatvec);
  printf("TRANSFER       : %e\n",traut_wall/nummatvec);
  printf("OTHERS         : %e\n",other_wall/nummatvec);
  printf("TOTAL          : %e\n",matve_wall/nummatvec);
  printf("\n");
  printf("SOLUTION TIME\n");
  printf("PREPROCESSING  : %e\n",prepr_wall);
  printf("SETUP          : %e\n",setup_wall);
  printf("SOLUTION       : %e\n",solut_wall);
  printf("MATVEC         : %e\n",matve_wall);
  printf("SOLVER         : %e\n",solut_wall-matve_wall);
  printf("POSTPROCESSING : %e\n",postp_wall);
  other_wall = total_wall-prepr_wall-setup_wall-solut_wall-postp_wall;
  printf("OTHER          : %e\n",other_wall);
  printf("TOTAL          : %e\n",total_wall);

  FILE *outf = fopen("times.txt","w");
  fprintf(outf,"%e\n",train_wall/nummatvec);
  fprintf(outf,"%e\n",aggrh_wall/nummatvec);
  for(int i = 0; i < level; i++)
  fprintf(outf,"%e\n",m2m_w[i]/nummatvec);
  fprintf(outf,"%e\n",commf_wall/nummatvec);
  fprintf(outf,"%e\n",trans_wall/nummatvec);
  for(int i = 0; i < level; i++)
  fprintf(outf,"%e\n",m2l_w[i]/nummatvec);
  fprintf(outf,"%e\n",disah_wall/nummatvec);
  for(int i = 0; i < level; i++)
  fprintf(outf,"%e\n",l2l_w[i]/nummatvec);
  fprintf(outf,"%e\n",(aggrh_wall+commf_wall+trans_wall+disah_wall)/nummatvec);
  fprintf(outf,"%e\n",commn_wall/nummatvec);
  fprintf(outf,"%e\n",nearf_wall/nummatvec);
  fprintf(outf,"%e\n",traut_wall/nummatvec);
  other_wall = matve_wall-train_wall-aggrh_wall-trans_wall-disah_wall-nearf_wall-traut_wall;
  fprintf(outf,"%e\n",other_wall/nummatvec);
  fprintf(outf,"%e\n",matve_wall/nummatvec);

  fprintf(outf,"%e\n",prepr_wall);
  fprintf(outf,"%e\n",setup_wall);
  fprintf(outf,"%e\n",solut_wall);
  fprintf(outf,"%e\n",matve_wall);
  fprintf(outf,"%e\n",solut_wall-matve_wall);
  fprintf(outf,"%e\n",postp_wall);
  other_wall = total_wall-prepr_wall-setup_wall-solut_wall-postp_wall;
  fprintf(outf,"%e\n",other_wall);
  fprintf(outf,"%e\n",total_wall);
  fclose(outf);
  }

  delete[] m2m_w;
  delete[] m2l_w;
  delete[] l2l_w;
  mem = mem - (double)sizeof(double)*3*level/1024/1024;
  if(myid==0)printf("FINAL MEMORY: %e GB\n",mem/1024);

  MPI_Finalize();
}
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return ((double)time.tv_sec + (double)time.tv_usec * .000001)*1000;
}
