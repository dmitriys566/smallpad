#pragma once
//Constants
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <vector>
#include <ios>
#include <string>
#include <array>
#include <assert.h>
#include <cassert>
#include <stdint.h>
#include <complex.h>
#include <functional>
#include "libamos.h"
#include "libeispack.h"
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <omp.h>
const double pi=3.14159265358979323846264338327;
const double e_const=     2.71828182845904523536028747135;
const double EulerGamma= 0.577215664901532860606512090082;
const double e0= 4.803242e-10;
const double m0= 9.10938188e-28;
const double hbar= 1.05457159642e-27;
const double kB= 1.3806503e-16;
const double evolt= 1.602176462e-12;
const double angstrem= 1e-8;
const double c_light= 2.99792458e10;
const double NaN= 0.0/0.0;
const double Inf= 1.0/0.0;
typedef double __complex__ double_complex_t;
#define NCPU 8

#ifdef __cplusplus
extern "C" {
#endif
typedef struct
{
  double re;
  double im;
}dcomplex;

extern void zgesv_( int* n, int* nrhs, dcomplex* a, int* lda, int* ipiv,
                dcomplex* b, int* ldb, int* info );

extern void dgtsv_(int*, int*, double*, double*, double*, double*, int*, int*);

#ifdef __cplusplus
}
#endif

double Zeta(double s);
double PolyLog(double s,double z);
double ClebschGordan(double j1,double m1,double j2,double m2,double j3,double m3);
int KronekerDelta(int l1,int l2);
double min(double a,double b);
double max(double a,double b);
double ChebyshevT(int n,double x);
double ChebyshevU(int n,double x);
double chebyshevt(double n,double x);
double chebyshevu(double n,double x);
double HermiteH(int n,double x);
double hermiteh(double n,double x);
double LaguerreL(int n,int m,double x);
double laguerrel(double n,double m,double x);
double EllipticE(double m);
double EllipticK(double m);
double_complex_t SphericalHarmonicY(int l,int m,double theta,double phi);
double Gamma(double x);
double ExpIntegralEi(double x);
double BinomCoeff(unsigned int n,unsigned int k);
double binom(double n,double k);
double LegendrePlm(int l,int m,double x);
double gcd(double a,double b);
double isprime(double x);
double pl(double n,double x);
double ql(double n,double x);
inline double plm(double l,double m,double x){
    return LegendrePlm((int)l,(int)m,x);
}
double GaussIntegrateParallel(double (*f)(double [],void *),void *serviceData,int ndim,double a[],double b[],int m,int nproc);
double_complex_t ZGaussIntegrateParallel(double_complex_t(*f)(double [],void *),void *serviceData,int ndim,double a[],double b[],int m,int nproc);
double GaussIntegrateElem(std::function<double(double *x,void *user_data)> f,void * serviceData,int ndim,double a[],double b[]);
double_complex_t ZGaussIntegrateElem(double_complex_t(*f)(double[],void *),void * serviceData,int ndim,double a[],double b[]);
double GaussIntegrate(std::function<double(double *x,void *user_data)> f,void * serviceData,int ndim,double a[],double b[],int m);
double_complex_t ZGaussIntegrate(double_complex_t(*f)(double[],void *),void * serviceData,int ndim,double a[],double b[],int m);
#ifndef CROSS_COMPILE
//double GaussIntegrateCuba(double (*f)(double[],void *),void * serviceData,int ndim,double a[],double b[],double epsrel,double epsabs);
//double_complex_t ZGaussIntegrateCuba(double_complex_t (*f)(double[],void *),void * serviceData,int ndim,double a[],double b[],double epsrel,double epsabs);
#endif
int IsNaN(double x);
int IsInf(double x);
int sign(double x);
double FindZero(std::function<double(double x,void *user_data)>,double a,double b,void *serviceData,int *errcode);
double FindNZero(std::function<double(double x,void *user_data)>,double x0,double x1,void *serviceData,double sep,int i,int * nf);
double EulerBeta(double x,double y);
double LegendreP(int n,double x);
double DLegendreP(int n,double x);
void MatrixMatrixMultiply(double *r,double *a,double *b,int m,int n,int k);
double SQR(double x);
double pythag_m(double a,double b );
void QR_decompose_rotation(double *a,double *q,double *r,int n);
void QR_decompose_reflection(double *a,double *q,double *r,int n);
void QR_solve(int N,double *q,double *r,double *b,double *x);
void printMatrix(double *a,char * title,int m,int n);
void printMatrixZ(double_complex_t *a,char * title,int m,int n);
#ifdef CONFIG_CRYPTO
double random_double();
#endif
int rk4(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t0,double dt,int nsteps,void *serviceData);
int rk4_step(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t,double dt,void *serviceData);
int rk5(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t0,double dt,int nsteps,void *serviceData);
int rk5_step(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t,double dt,void *serviceData);
double ipow(double x,int k);
double_complex_t ZGammaIncomplete(double s,double_complex_t z);
double_complex_t ZExpIntegralE(double s,double_complex_t z);
void print_complex(char * msg,double_complex_t z);
double_complex_t besselj(double nu, double_complex_t z);
double_complex_t bessely(double nu, double_complex_t z);
double_complex_t besseli(double nu, double_complex_t z);
double_complex_t besselk(double nu, double_complex_t z);
double dbesselj(double nu, double z);
double dbessely(double nu, double z);
double dbesseli(double nu, double z);
double dbesselk(double nu, double z);
double sj(int l,double x);
double sy(int l,double x);
double si(int l,double x);
double sk(int l,double x);
double AiryAi(double x,int derivative);
double AiryBi(double x,int derivative);
double ai(double x);
double bi(double x);
double aiprime(double x);
double biprime(double x);

double Pochgammer(int n,double alpha);
double Hypergeometric1F1(double a,double c,double z);
double WhittakerM(double lambda,double mu,double z);
double ifact(int k);
int EigenSystem(int n,double *a,std::vector<double_complex_t> &w,std::vector<std::vector<double_complex_t>> &zc,int eigenvectors,int nproc=1);
int EigenSystemOrig(int n,double *a,double_complex_t *w,double_complex_t *zc,int eigen_vectors);
void print_complex_vector(char * msg,int n,double_complex_t *z);
void print_complex_matrix(char * msg,int n,double_complex_t *z);

bool find_root(std::function<std::vector<double>(std::vector<double> &x,void *user_data)> f, const std::vector<double> &x_init, std::vector<double> &x_out, void *user_data, double epsilon=1e-10, int max_iter=10000);
bool find_minimum_descent(std::function<double(std::vector<double> &x,void *user_data)> f, const std::vector<double> &x_init, std::vector<double> &x_out, void *user_data, double lambda_init=0.01, double epsilon=1e-10, int max_iter=10000);
bool find_minimum_lbfgs(std::function<double(std::vector<double> &x,void *user_data)> f, const std::vector<double> &x_init, std::vector<double> &x_out, void *user_data, double lambda_init=0.01, double epsilon=1e-10, int max_iter=10000);
bool find_minimum_1d(std::function<double (double &, void *)> f, const double &x_init, double &x_out, void *user_data, double epsilon=1e-10, int max_iter=10000);

