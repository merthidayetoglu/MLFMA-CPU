#include <cstdio>
#include <cmath>
#include <complex>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#include <sys/time.h>

using namespace std;

void preproc();
void setup();
void final();
void setup_mpi();

complex<double> integrate(complex<double>, complex<double>);
complex<double> integrate_multi(complex<double>, complex<double>, int);
complex<double> integrate_local(complex<double>, complex<double>, int);
complex<double> hn(int, double);

void bicgs(complex<double>*,complex<double>*,complex<double>*);
void bicgs_gpu(complex<double>*,complex<double>*,complex<double>*);
void born(complex<double>*,complex<double>*,complex<double>*);
void conjg(complex<double>*,complex<double>*,complex<double>*);

void mlfma(complex<double>*,complex<double>*);
void mlfma_mpi(complex<double>*,complex<double>*);

void farfield(complex<double>*);

double get_wall_time();
double get_cpu_time();
