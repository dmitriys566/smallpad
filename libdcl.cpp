#include "libdcl.h"
#include <stdint.h>
#include <cmath>
#include "libamos.h"
#define N_max 200
#define N_trunc 15
#define nu_max 15

//Cuba defines
#define VERBOSE 0
#define LAST 4
#define MINEVAL 200000
#define MAXEVAL 1000000000
#define KEY 0
#define NCOMP 1
#define STATEFILE NULL
#define SPIN NULL

typedef struct{
    double s;
    double z;
}polylog_param;

typedef struct{
    double_complex_t p;
    double k;
}serviceData;

typedef struct{
    double mu;
    double lambda;
    double z;
}w_params;

typedef struct{
    double alpha;
    double beta1;
    double beta2;
    double gamma1;
    double gamma2;
    double x;
    double y;
}f2_params;

typedef struct{
    double (*f)(double *,void *);
    double **x_glob;
    double *c_glob;
    void *serviceData;
    int ndim;
    long long int size;
    int nproc;
    int tid;
}pt_params;

typedef struct{
    double_complex_t(*f)(double *,void *);
    double **x_glob;
    double *c_glob;
    void *serviceData;
    int ndim;
    long long int size;
    int nproc;
    int tid;
}pt_params_z;

//This function is used to increment indexes in Gauss integrate//
//TODO fix recursion in orthogonal polynoms.
//TODO add parallelism in integration
int increment(int *array,int N,int nmax){
    int i,carry;
    carry=1;
    for(i=0;(i<N)&&carry;i++){
        if(array[i]==nmax){
            array[i]=0;
            carry=1;
        }else{
            array[i]++;
            carry=0;
        }
    }
    return carry;
}

//USE
//Use the Lanczcos approximation
//TODO: fix the coefficients
double Gamma(double x){
    static const double d[15]=
    {0.011748482397059865854194497518649980956864516461236210282662113866005188,
     0.67149902860259779045990250413501961046434632677594141601778811910717214,
     -0.70018558813697759056405597501837691149234060570108385068405304408801,
     0.1660776982193675198571572139135855322686605807017872543338482203672,
     -0.00577924080929345179860093218662419128059386859373306997961415072,
     3.993855469407750719798353396430172688531170874259718037933e-7,
     5.46582035496776954711515712918817773446478503007554770768e-7,
     -1.155750791439381406933196280628910251054174145609947511436e-6,
     1.85730234701191102358179289192258099062275696389504662339e-6,
     -2.4702880923232733952229889270160370646909100217280405525e-6,
     2.55458552585000269416928843452846310351308677477529019125e-6,
     -1.93048838216538385557473513898173458259530760512112460056e-6,
     9.9178601846535382019593795774495103059935300410647600095e-7,
     -3.07702603925219628485419424700998229231072885973688462704e-7,
     4.3350939794140517541113396050506568511274887170176440299e-8};
    static const double r=607.0/128.0;
    double rv;
    double sum;
    int i;
    if(x<0.5){
        return (pi/sin(pi*x))/Gamma(1-x);
    }
    //Uncomment if you want better precision. This is slower.
    /*if(x>3.0){
    return (x-1)*Gamma(x-1);
  }*/
    rv=2.0*sqrt(e_const/pi)*pow((x+r-0.5)/e_const,x-0.5);
    sum=d[0];
    for(i=1;i<15;i++){
        sum+=d[i]/(x-1+i);
    }
    rv=rv*sum;
    return rv;
}
//USE
double LegendrePlm(int l,int m,double x){
    double fact1,fact2,rv,sq;
    int i;
    if((abs(m)>abs(l)&&l>=0)||(abs(m)>abs(-l-1)&&l<0)){
        return 0;
    }
    if(l<0){
        return LegendrePlm(-l-1,m,x);
    }
    if(m<0){
        fact1=1.0;
        fact2=1.0;
        rv=1.0;
        for(i=1;i<=l+m;i++){
            fact1*=i;
        }
        for(i=1;i<=l-m;i++){
            fact2*=i;
        }
        if(m%2){
            rv=-rv;
        }
        return rv*fact1*LegendrePlm(l,-m,x)/fact2;
    }
    if(m==l){
        rv=1.0;
        if(m%2){
            rv=-rv;
        }
        for(i=1;i<=(2*m-1);i+=2){
            rv*=i;
        }
        sq=sqrt(1-x*x);
        for(i=0;i<m;i++){
            rv*=sq;
        }
        return rv;
    }
    return (x*(2.0*l-1)*LegendrePlm(l-1,m,x)-(l+m-1)*LegendrePlm(l-2,m,x))/(l-m);
}
//USE
double EllipticE(double m){
    double a[20];
    double b[20];
    double c[20];
    double k,e,epsilon,sum,pow;
    int N=19;
    int i;
    if(m<0||m>1){
        return 0;
    }
    if(m==1.0){
        k=1.0/0.0;
        e=1.0;
    }else{
        a[0]=1.0;
        b[0]=sqrt(1.0-m);
        c[0]=sqrt(m);
        //err=1.0;
        epsilon=1e-16;
        for(i=1;i<=N;i++){
            a[i]=(a[i-1]+b[i-1])*0.5;
            b[i]=sqrt(a[i-1]*b[i-1]);
            c[i]=(a[i-1]-b[i-1])*0.5;
            //err=c[i];
        }
        i--;
        if(fabs(c[i])>epsilon){
            printf("Elliptic does not converge after %d steps. STOP\n",N);
            return 0;
        }
        k=pi/(2*a[i]);
        sum=0;
        pow=1;
        for(i=0;i<=N;i++){
            sum+=c[i]*c[i]*pow;
            pow*=2;
        }
        sum*=0.5;
        e=(1-sum)*k;
    }
    return e;
}

//USE
double EllipticK(double m){
    double a[20];
    double b[20];
    double c[20];
    double k,e,epsilon,sum,pow;
    int N=19;
    int i;
    if(m<0||m>1){
        return 0;
    }
    if(m==1.0){
        k=1.0/0.0;
        e=1.0;
    }else{
        a[0]=1.0;
        b[0]=sqrt(1.0-m);
        c[0]=sqrt(m);
        //err=1.0;
        epsilon=1e-16;
        for(i=1;i<=N;i++){
            a[i]=(a[i-1]+b[i-1])*0.5;
            b[i]=sqrt(a[i-1]*b[i-1]);
            c[i]=(a[i-1]-b[i-1])*0.5;
            //err=c[i];
        }
        i--;
        if(fabs(c[i])>epsilon){
            printf("Elliptic does not converge after %d steps. STOP\n",N);
            return 0;
        }
        k=pi/(2*a[i]);
        sum=0;
        pow=1;
        for(i=0;i<=N;i++){
            sum+=c[i]*c[i]*pow;
            pow*=2;
        }
        sum*=0.5;
        e=(1-sum)*k;
    }
    return k;
}

//NOT USE
double_complex_t SphericalHarmonicY(int l,int m,double theta,double phi){
    double fact1,fact2;
    int i;
    if(l<0){
        return 0;
    }
    fact1=1;
    fact2=1;
    for(i=1;i<=(l-m);i++){
        fact1*=i;
    }
    for(i=1;i<=(l+m);i++){
        fact2*=i;
    }
    return sqrt((2.0*l+1)*fact1/(4*pi*fact2))*LegendrePlm(l,m,cos(theta))*cexp(I*m*phi);
}

//USE
double HermiteH(int n,double x){
    if(n<0){
        return 0;
    }
    if(n==0){
        return 1;
    }
    if(n==1){
        return 2*x;
    }
    return 2*x*HermiteH(n-1,x)-2*(n-1)*HermiteH(n-2,x);
}

//USE
double hermiteh(double n,double x){
    return HermiteH((int)n,x);
}
//USE
double ChebyshevT(int n,double x){
    if(n<0){
        return 0;
    }
    if(n==0){
        return 1;
    }
    if(n==1){
        return x;
    }
    return 2*x*ChebyshevT(n-1,x)-ChebyshevT(n-2,x);
}
//USE
double ChebyshevU(int n,double x){
    if(n<0){
        return 0;
    }
    if(n==0){
        return 1;
    }
    if(n==1){
        return 2*x;
    }
    return 2*x*ChebyshevU(n-1,x)-ChebyshevU(n-2,x);
}
//USE
double chebyshevt(double n,double x){
    return ChebyshevT((int)n,x);
}
//USE
double chebyshevu(double n,double x){
    return ChebyshevU((int)n,x);
}
//USE
double BinomCoeff(unsigned int n,unsigned int k){
    int i;
    double rv;
    rv=1.0;
    if(k<=n){
        for(i=1;i<=n;i++){
            rv*=i;
        }
        for(i=1;i<=k;i++){
            rv/=i;
        }
        for(i=1;i<=(n-k);i++){
            rv/=i;
        }
    }else{
        rv=0.0;
    }
    return rv;
}

double binom(double n,double k){
    return BinomCoeff((unsigned int)n,(unsigned int)k);
}
//USE
double ExpIntegralEi(double x){
    double epsilon,xm,rv,sum,fact,term,prev;
    int i,maxiter;
    epsilon=1e-17;
    xm=fabs(log(epsilon));
    maxiter=120;
    if(fabs(x)<xm){
        //Do power series
        rv=EulerGamma+log(fabs(x));
        sum=0;
        fact=1.0;
        for(i=1;i<maxiter;i++){
            fact*=x/i;
            term=fact/i;
            sum+=term;
            if(fabs(term)<epsilon*sum) break;
        }
        rv+=sum;
    }else{
        //Do asymptotic expansion
        sum=0;
        term=1.0;
        for(i=1;i<maxiter;i++){
            prev=term;
            term*=i/x;
            if(fabs(term)<epsilon) break;
            if(fabs(term)<fabs(prev)){
                sum+=term;
            }else{
                sum-=prev;
                break;
            }
        }
        rv=exp(x)*(1.0+sum)/x;
    }
    return rv;
}

//m difines number of sections;
double GaussIntegrate(std::function<double(double *x,void *user_data)> f,void *serviceData,int ndim,double a[],double b[],int m){
    int *idx;
    int carry;
    double rv;
    double *a1,*b1;
    double delta;
    int i;
    idx=(int *)malloc(ndim*sizeof(int));
    a1=(double *)malloc(ndim*sizeof(double));
    b1=(double *)malloc(ndim*sizeof(double));
    carry=0;
    rv=0;
    for(i=0;i<ndim;i++){
        idx[i]=0;
    }
    delta=1.0/(double)m;
    while(!carry){
        for(i=0;i<ndim;i++){
            a1[i]=delta*idx[i]*(b[i]-a[i])+a[i];
            b1[i]=(idx[i]+1)*delta*(b[i]-a[i])+a[i];
        }
        rv+=GaussIntegrateElem(f,serviceData,ndim,a1,b1);
        carry=increment(idx,ndim,m-1);
    }
    free(a1);
    free(b1);
    free(idx);
    return rv;
}

//m difines number of sections;
double_complex_t ZGaussIntegrate(double_complex_t(*f)(double[],void *),void *serviceData,int ndim,double a[],double b[],int m){
    int *idx;
    int carry;
    double_complex_t rv;
    double *a1,*b1;
    double delta;
    int i;
    idx=(int *)malloc(ndim*sizeof(int));
    a1=(double *)malloc(ndim*sizeof(double));
    b1=(double *)malloc(ndim*sizeof(double));
    carry=0;
    rv=0;
    for(i=0;i<ndim;i++){
        idx[i]=0;
    }
    delta=1.0/(double)m;
    while(!carry){
        for(i=0;i<ndim;i++){
            a1[i]=delta*idx[i]*(b[i]-a[i])+a[i];
            b1[i]=(idx[i]+1)*delta*(b[i]-a[i])+a[i];
        }
        rv+=ZGaussIntegrateElem(f,serviceData,ndim,a1,b1);
        carry=increment(idx,ndim,m-1);
    }
    free(a1);
    free(b1);
    free(idx);
    return rv;
}


#ifndef CROSS_COMPILE
/*double GaussIntegrateCuba(double (*f)(double[],void *),void *serviceData,int ndim,double a[],double b[],double epsrel,double epsabs){
  int nvec_glob=1;
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  int integrand(const int *ndim_intern, const cubareal xx[],const int *ncomp, cubareal ff[], void *userdata) {
    cubareal *t;
    cubareal jacobian=1;
    int i;
    t=(cubareal *)malloc(ndim*sizeof(cubareal));//Translated to [a[i],b[i]]
    for(i=0;i<ndim;i++){
      t[i]=(xx[i])*(b[i]-a[i])+a[i];
      jacobian*=b[i]-a[i];
    }
    ff[0]=jacobian*f(t,userdata);
  }
  Cuhre(ndim, NCOMP, integrand, serviceData, nvec_glob,epsrel, epsabs, VERBOSE | LAST,MINEVAL,MAXEVAL, KEY,STATEFILE, SPIN,&nregions, &neval, &fail, integral, error, prob);
  return integral[0];
}

double_complex_t ZGaussIntegrateCuba(double_complex_t (*f)(double[],void *),void *serviceData,int ndim,double a[],double b[],double epsrel,double epsabs){
  int nvec_glob=1;
  int comp, nregions, neval, fail;
  cubareal integral_r[NCOMP],integral_i[NCOMP], error[NCOMP], prob[NCOMP];

  int integrand_real(const int *ndim_intern, const cubareal xx[],const int *ncomp, cubareal ff[], void *userdata) {
    cubareal *t;
    cubareal jacobian=1;
    int i;
    t=(cubareal *)malloc(ndim*sizeof(cubareal));//Translated to [a[i],b[i]]
    for(i=0;i<ndim;i++){
      t[i]=(xx[i])*(b[i]-a[i])+a[i];
      jacobian*=b[i]-a[i];
    }
    ff[0]=jacobian*creal(f(t,userdata));
  }
  int integrand_imag(const int *ndim_intern, const cubareal xx[],const int *ncomp, cubareal ff[], void *userdata) {
    cubareal *t;
    cubareal jacobian=1;
    int i;
    t=(cubareal *)malloc(ndim*sizeof(cubareal));//Translated to [a[i],b[i]]
    for(i=0;i<ndim;i++){
      t[i]=(xx[i])*(b[i]-a[i])+a[i];
      jacobian*=b[i]-a[i];
    }
    ff[0]=jacobian*cimag(f(t,userdata));
  }
  Cuhre(ndim, NCOMP, integrand_real, serviceData, nvec_glob,epsrel, epsabs, VERBOSE | LAST,MINEVAL,MAXEVAL, KEY,STATEFILE, SPIN,&nregions, &neval, &fail, integral_r, error, prob);
  Cuhre(ndim, NCOMP, integrand_imag, serviceData, nvec_glob,epsrel, epsabs, VERBOSE | LAST,MINEVAL,MAXEVAL, KEY,STATEFILE, SPIN,&nregions, &neval, &fail, integral_i, error, prob);
  return integral_r[0]+I*integral_i[0];
}*/
#endif

double GaussIntegrateElem(std::function<double(double *x,void *user_data)> f,void *serviceData,int ndim,double a[],double b[]){
    //define Gauss rule;
    static const double ksi[15]={-0.9879925180204854,-0.937273392400706,-0.8482065834104272,-0.7244177313601701,-0.5709721726085388,-0.3941513470775634,-0.2011940939974345,0,0.2011940939974345,0.3941513470775634,0.5709721726085388,0.7244177313601701,0.8482065834104272,0.937273392400706,0.9879925180204854};
    static const double an[15]={0.03075324199611749,0.07036604748810815,0.107159220467172,0.1395706779261543,0.1662692058169939,0.1861610000155622,0.1984314853271116,0.2025782419255613,0.1984314853271116,0.1861610000155622,0.1662692058169939,0.1395706779261543,0.107159220467172,0.07036604748810815,0.03075324199611749};
    int nmax=14;
    int *idx;
    double *x;
    double rv,eta;
    int i,j,carry;
    idx=(int *)malloc(ndim*sizeof(int));
    x=(double *)malloc(ndim*sizeof(double));
    for(i=0;i<ndim;i++){
        idx[i]=0;
    }
    rv=0.0;
    carry=0;
    while(!carry){
        for(i=0;i<ndim;i++){
            x[i]=ksi[idx[i]]*(b[i]-a[i])/2.0+(a[i]+b[i])/2.0;
        }
        eta=f(x,serviceData);
        for(i=0;i<ndim;i++){
            eta=eta*an[idx[i]]*(b[i]-a[i])/2.0;
        }
        rv+=eta;
        carry=increment(idx,ndim,nmax);
    }
    free(idx);
    free(x);
    return rv;
}

double_complex_t ZGaussIntegrateElem(double_complex_t(*f)(double[],void *),void *serviceData,int ndim,double a[],double b[]){
    //define Gauss rule;
    static const double ksi[15]={-0.9879925180204854,-0.937273392400706,-0.8482065834104272,-0.7244177313601701,-0.5709721726085388,-0.3941513470775634,-0.2011940939974345,0,0.2011940939974345,0.3941513470775634,0.5709721726085388,0.7244177313601701,0.8482065834104272,0.937273392400706,0.9879925180204854};
    static const double an[15]={0.03075324199611749,0.07036604748810815,0.107159220467172,0.1395706779261543,0.1662692058169939,0.1861610000155622,0.1984314853271116,0.2025782419255613,0.1984314853271116,0.1861610000155622,0.1662692058169939,0.1395706779261543,0.107159220467172,0.07036604748810815,0.03075324199611749};
    int nmax=14;
    int *idx;
    double *x;
    double_complex_t rv,eta;
    int i,j,carry;
    idx=(int *)malloc(ndim*sizeof(int));
    x=(double *)malloc(ndim*sizeof(double));
    for(i=0;i<ndim;i++){
        idx[i]=0;
    }
    rv=0.0;
    carry=0;
    while(!carry){
        for(i=0;i<ndim;i++){
            x[i]=ksi[idx[i]]*(b[i]-a[i])/2.0+(a[i]+b[i])/2.0;
        }
        eta=f(x,serviceData);
        for(i=0;i<ndim;i++){
            eta=eta*an[idx[i]]*(b[i]-a[i])/2.0;
        }
        rv+=eta;
        carry=increment(idx,ndim,nmax);
    }
    free(idx);
    free(x);
    return rv;
}

void * pt_proc(void *params){
    pt_params *p;
    double **x_glob;
    double *c_glob;
    double (*f)(double *,void *);
    void *serviceData;
    int ndim,tid,nproc;
    long long int size,i;
    double rv=0;
    double *rvptr;
    p=(pt_params *)params;
    x_glob=p->x_glob;
    c_glob=p->c_glob;
    f=p->f;
    ndim=p->ndim;
    size=p->size;
    tid=p->tid;
    nproc=p->nproc;
    serviceData=p->serviceData;
    for(i=tid;i<size;i+=nproc){
        rv+=f(x_glob[i],serviceData)*c_glob[i];
    }
    rvptr=(double *)malloc(sizeof(double));
    *rvptr=rv;
    return (void *)rvptr;
}

void * pt_proc_z(void *params){
    pt_params_z *p;
    double **x_glob;
    double *c_glob;
    double_complex_t (*f)(double *,void *);
    void *serviceData;
    int ndim,tid,nproc;
    long long int size,i;
    double_complex_t rv=0;
    double_complex_t *rvptr;
    p=(pt_params_z *)params;
    x_glob=p->x_glob;
    c_glob=p->c_glob;
    f=p->f;
    ndim=p->ndim;
    size=p->size;
    tid=p->tid;
    nproc=p->nproc;
    serviceData=p->serviceData;
    for(i=tid;i<size;i+=nproc){
        rv+=f(x_glob[i],serviceData)*c_glob[i];
    }
    rvptr=(double_complex_t *)malloc(sizeof(double_complex_t));
    *rvptr=rv;
    return (void *)rvptr;
}

//Realize cubature of Gauss Rule using pthreads.
double GaussIntegrateParallel(double (*f)(double [],void *),void *serviceData,int ndim,double a[],double b[],int m,int nproc){
    pthread_t *tids;
    pt_params *tparams;
    int *idx,*idx_current;
    int carry,carry_current;
    double rv,eta;
    double *a1,*b1;
    double *x;
    double **x_glob,*c_glob;
    double **retval;
    long long int size;
    double delta;
    int i,j;
    //define Gauss rule;
    static const double ksi[15]={-0.9879925180204854,-0.937273392400706,-0.8482065834104272,-0.7244177313601701,-0.5709721726085388,-0.3941513470775634,-0.2011940939974345,0,0.2011940939974345,0.3941513470775634,0.5709721726085388,0.7244177313601701,0.8482065834104272,0.937273392400706,0.9879925180204854};
    static const double an[15]={0.03075324199611749,0.07036604748810815,0.107159220467172,0.1395706779261543,0.1662692058169939,0.1861610000155622,0.1984314853271116,0.2025782419255613,0.1984314853271116,0.1861610000155622,0.1662692058169939,0.1395706779261543,0.107159220467172,0.07036604748810815,0.03075324199611749};
    int nmax=14;
    long long int counter=0;
    idx_current=(int *)malloc(ndim*sizeof(int));
    x=(double *)malloc(ndim*sizeof(double));
    idx=(int *)malloc(ndim*sizeof(int));
    a1=(double *)malloc(ndim*sizeof(double));
    b1=(double *)malloc(ndim*sizeof(double));
    retval=(double **)malloc(nproc*sizeof(double *));
    tids=(pthread_t *)malloc(nproc*sizeof(pthread_t));
    tparams=(pt_params *)malloc(nproc*sizeof(pt_params));
    carry=0;
    carry_current=0;
    rv=0;
    for(i=0;i<ndim;i++){
        idx[i]=0;
        idx_current[i]=0;
    }
    size=1;
    for(i=0;i<ndim;i++){
        size*=m;
        size*=nmax+1;
    }
    x_glob=(double **)malloc(size*sizeof(double *));
    c_glob=(double *)malloc(size*sizeof(double));
    for(counter=0;counter<size;counter++){
        x_glob[counter]=(double *)malloc(ndim*sizeof(double));
        c_glob[counter]=1.0;
    }
    delta=1.0/(double)m;
    counter=0;
    while(!carry){
        //Initial array setup;
        for(i=0;i<ndim;i++){
            a1[i]=delta*idx[i]*(b[i]-a[i])+a[i];
            b1[i]=(idx[i]+1)*delta*(b[i]-a[i])+a[i];
        }
        carry_current=0;
        for(i=0;i<ndim;i++){
            idx_current[i]=0;//Initial setup;
        }
        while(!carry_current){
            for(i=0;i<ndim;i++){
                x[i]=ksi[idx_current[i]]*(b1[i]-a1[i])/2.0+(a1[i]+b1[i])/2.0;
                x_glob[counter][i]=x[i];
            }
            //      eta=f(x,serviceData);
            for(i=0;i<ndim;i++){
                c_glob[counter]*=an[idx_current[i]]*(b1[i]-a1[i])/2.0;
            }
            //      rv+=eta;
            counter++;
            carry_current=increment(idx_current,ndim,nmax);
        }
        carry=increment(idx,ndim,m-1);
    }
    //Fill the structures;
    for(i=0;i<nproc;i++){
        tparams[i].f=f;
        tparams[i].x_glob=x_glob;
        tparams[i].c_glob=c_glob;
        tparams[i].serviceData=serviceData;
        tparams[i].ndim=ndim;
        tparams[i].size=size;
        tparams[i].nproc=nproc;
        tparams[i].tid=i;
        pthread_create(&tids[i],0,pt_proc,(void *)&tparams[i]);
    }
    for(i=0;i<nproc;i++){
        pthread_join(tids[i],(void **)&retval[i]);
    }
    rv=0;
    for(i=0;i<nproc;i++){
        rv+=*retval[i];
    }
    free(a1);
    free(b1);
    free(idx);
    free(idx_current);
    free(x_glob);
    free(c_glob);
    free(tparams);
    free(tids);
    free(x);
    return rv;
}

double gcd(double a, double b)
{
    int a_,b_,r,gc;
    a_=(int)a;
    b_=(int)b;
    if(abs(a_)<abs(b_)){
        r=a_;
        a_=b_;
        b_=r;
    }
    gc=abs(b);
    while(r!=0){
        r=a_%b_;
        if(r!=0){
            gc=r;
        }
        a_=b_;
        b_=r;
    }
    return gc;
}

double isprime(double x)
{
    int x_;
    double rv=1;
    x_=(int)x;
    for(int i=2;i*i<=x_;i++){
        if(x_%i==0){
            rv=0;
        }
    }
    return rv;
}
//Realize cubature of Gauss Rule using pthreads.
double_complex_t ZGaussIntegrateParallel(double_complex_t(*f)(double [],void *),void *serviceData,int ndim,double a[],double b[],int m,int nproc){
    pthread_t *tids;
    pt_params_z *tparams;
    int *idx,*idx_current;
    int carry,carry_current;
    double_complex_t rv,eta;
    double *a1,*b1;
    double *x;
    double **x_glob,*c_glob;
    double_complex_t **retval;
    long long int size;
    double delta;
    int i,j;
    //define Gauss rule;
    static const double ksi[15]={-0.9879925180204854,-0.937273392400706,-0.8482065834104272,-0.7244177313601701,-0.5709721726085388,-0.3941513470775634,-0.2011940939974345,0,0.2011940939974345,0.3941513470775634,0.5709721726085388,0.7244177313601701,0.8482065834104272,0.937273392400706,0.9879925180204854};
    static const double an[15]={0.03075324199611749,0.07036604748810815,0.107159220467172,0.1395706779261543,0.1662692058169939,0.1861610000155622,0.1984314853271116,0.2025782419255613,0.1984314853271116,0.1861610000155622,0.1662692058169939,0.1395706779261543,0.107159220467172,0.07036604748810815,0.03075324199611749};
    int nmax=14;
    long long int counter=0;
    idx_current=(int *)malloc(ndim*sizeof(int));
    x=(double *)malloc(ndim*sizeof(double));
    idx=(int *)malloc(ndim*sizeof(int));
    a1=(double *)malloc(ndim*sizeof(double));
    b1=(double *)malloc(ndim*sizeof(double));
    retval=(double_complex_t **)malloc(nproc*sizeof(double_complex_t*));
    tids=(pthread_t *)malloc(nproc*sizeof(pthread_t));
    tparams=(pt_params_z *)malloc(nproc*sizeof(pt_params_z));
    carry=0;
    carry_current=0;
    rv=0;
    for(i=0;i<ndim;i++){
        idx[i]=0;
        idx_current[i]=0;
    }
    size=1;
    for(i=0;i<ndim;i++){
        size*=m;
        size*=nmax+1;
    }
    x_glob=(double **)malloc(size*sizeof(double *));
    c_glob=(double *)malloc(size*sizeof(double));
    for(counter=0;counter<size;counter++){
        x_glob[counter]=(double *)malloc(ndim*sizeof(double));
        c_glob[counter]=1.0;
    }
    delta=1.0/(double)m;
    counter=0;
    while(!carry){
        //Initial array setup;
        for(i=0;i<ndim;i++){
            a1[i]=delta*idx[i]*(b[i]-a[i])+a[i];
            b1[i]=(idx[i]+1)*delta*(b[i]-a[i])+a[i];
        }
        carry_current=0;
        for(i=0;i<ndim;i++){
            idx_current[i]=0;//Initial setup;
        }
        while(!carry_current){
            for(i=0;i<ndim;i++){
                x[i]=ksi[idx_current[i]]*(b1[i]-a1[i])/2.0+(a1[i]+b1[i])/2.0;
                x_glob[counter][i]=x[i];
            }
            //      eta=f(x,serviceData);
            for(i=0;i<ndim;i++){
                c_glob[counter]*=an[idx_current[i]]*(b1[i]-a1[i])/2.0;
            }
            //      rv+=eta;
            counter++;
            carry_current=increment(idx_current,ndim,nmax);
        }
        carry=increment(idx,ndim,m-1);
    }
    //Fill the structures;
    for(i=0;i<nproc;i++){
        tparams[i].f=f;
        tparams[i].x_glob=x_glob;
        tparams[i].c_glob=c_glob;
        tparams[i].serviceData=serviceData;
        tparams[i].ndim=ndim;
        tparams[i].size=size;
        tparams[i].nproc=nproc;
        tparams[i].tid=i;
        pthread_create(&tids[i],0,pt_proc_z,(void *)&tparams[i]);
    }
    for(i=0;i<nproc;i++){
        pthread_join(tids[i],(void **)&retval[i]);
    }
    rv=0;
    for(i=0;i<nproc;i++){
        rv+=*retval[i];
    }
    free(a1);
    free(b1);
    free(idx);
    free(idx_current);
    free(x_glob);
    free(c_glob);
    free(tparams);
    free(tids);
    free(x);
    return rv;
}

int IsNaN(double x){
    int rv;
    rv=0;
    if(x>0.0){
        rv=0;
    }else{
        if(x<=0.0){
            rv=0;
        }else{
            rv=1;
        }
    }
    return rv;
}

int IsInf(double x){
    return 1.0/x==0;
}

int sign(double x){
    if(x>0.0){
        return 1;
    }
    if(x<0.0){
        return -1;
    }
    if(x==0.0){
        return 0;
    }
}

int KronekerDelta(int l1,int l2){
    if(l1==l2){
        return 1;
    }
    return 0;
}

//TODO:fix n<0;
double LegendreP(int n,double x){
    int i;
    double arr[3];
    arr[0]=1.0;
    arr[1]=x;
    for(i=1;i<n;i++){
        arr[2]=((2.0*i+1.0)*x*arr[1]-i*arr[0])/(i+1.0);
        arr[0]=arr[1];
        arr[1]=arr[2];
    }
    if(n==0){
        return arr[0];
    }else{
        return arr[1];
    }
}

double LegendreQ(int n,double x){
    int i;
    double arr[3];
    arr[0]=1.0;
    arr[1]=x;
    for(i=1;i<n;i++){
        arr[2]=((2.0*i+1.0)*x*arr[1]-i*arr[0])/(i+1.0);
        arr[0]=arr[1];
        arr[1]=arr[2];
    }
    if(n==0){
        return arr[0];
    }else{
        return arr[1];
    }
}

double pl(double n,double x){
    return LegendreP((int)n,x);
}

double ql(double n,double x){
    return LegendreQ((int)n,x);
}

double DLegendreP(int n,double x){
    return (LegendreP(n+1,x)-x*LegendreP(n,x))*(n+1)/(x*x-1.0);
}

double FindZero(std::function<double(double x,void *user_data)> f,double a,double b,void *serviceData,int * errcode){
    double z1,z2,z3,znew,zold;
    double f1,f2,f3;
    double epsilon,er_rel;
    double qj,Aj,Bj,Cj;
    int maxiter,cnt,cont_iter;
    maxiter=100;
    cnt=0;
    *errcode=0;
    epsilon=1e-18;
    cont_iter=1;
    if((!IsNaN(f(a,serviceData)))&&(!IsNaN(f(b,serviceData)))){
        if(f(a,serviceData)*f(b,serviceData)>=0.0){
            if(f(a,serviceData)==0.0){
                return a;
            }
            if(f(b,serviceData)==0.0){
                return b;
            }
            *errcode=1;
            return 0.0/0.0;
        }else{
            z1=a;
            z2=b;
            z3=(a+b)/2.0;
            zold=z1;
            for(cnt=0;cnt<maxiter&&cont_iter;cnt++){
                f1=f(z1,serviceData);
                f2=f(z2,serviceData);
                f3=f(z3,serviceData);
                qj=(z3-z2)/(z2-z1);
                Aj=qj*f3-qj*(1.0+qj)*f2+f1*qj*qj;
                Bj=(2.0*qj+1.0)*f3-(1.0+qj)*(1.0+qj)*f2+f1*qj*qj;
                Cj=f3*(1.0+qj);
                znew=z3-sign(Bj)*(z3-z2)*(2.0*Cj)/(fabs(Bj)+sqrt(Bj*Bj-4.0*Aj*Cj));
                if(IsNaN(znew)){
                    cont_iter=0;
                }else{
                    if(f1*f3>0.0){
                        z1=z2;
                        z2=z3;
                        z3=znew;
                    }else{
                        z2=z3;
                        z3=znew;
                    }
                }
            }
            er_rel=fabs((znew-zold)/znew);
            if(er_rel<epsilon){
                cont_iter=0;
            }
            zold=znew;
        }
        return z2;
    }else{
        *errcode=1;
        return 0.0/0.0;
    }
}

//Find i-th zero of function f;
double FindNZero(std::function<double(double x,void *user_data)> f,double x0,double x1,void *serviceData,double sep,int i,int * nf){
    double x;
    double a,b,eps;
    int cnt;
    *nf=0;
    a=x0;
    b=x0+sep;
    cnt=0;
    eps=1e-10;
    if(x0>x1){
        *nf=1;
        printf("Error. Check bounds. Exit\n");
        return 0.0/0.0;
    }
    if(i==0){
        if(f(x0,serviceData)==0.0){
            return x0;
        }else{
            *nf=1;
            return 0.0/0.0;
        }
    }
    while((cnt<i)&&(a<x1)){
        if(b>=x1){
            b=x1;
            if(f(a,serviceData)*f(b,serviceData)<0 || f(b,serviceData)==0){
                cnt++;
            }
        }else{
            if(f(a,serviceData)*f(b,serviceData)<0 || f(b,serviceData)==0){
                cnt++;
            }
            b=b+sep;
        }
        a=a+sep;
    }
    if(cnt<i){
        *nf=1;
        return 0.0/0.0;
    }else{
        if(a<x1){
            b=b-sep;
        }
        a=a-sep;
        return FindZero(f,a,b,serviceData,nf);
    }
}

double LaguerreL(int n,int m,double x){
    if(n<0){
        return 0;
    }
    if(n==0){
        return 1;
    }
    if(n==1){
        return m+1-x;
    }
    return ((2*(n-1)+m+1-x)*(LaguerreL(n-1,m,x))-(n-1+m)*LaguerreL(n-2,m,x))/(n);
}

double laguerrel(double n,double m,double x){
    return LaguerreL((int)n,(int)m,x);
}

//Input: a[m,n];b[n,k]
//Output: r[m,k]
void MatrixMatrixMultiply(double *r,double *a,double *b,int m,int n,int k){
    int i,j,cnt;
    double sum;
    for(i=0;i<m;i++){
        for(j=0;j<k;j++){
            sum=0;
            for(cnt=0;cnt<n;cnt++){
                sum+=a[i*n+cnt]*b[cnt*k+j];
            }
            r[i*k+j]=sum;
        }
    }
    return;
}

double SQR(double x){
    return x*x;
}

double pythag_m ( double a, double b )
{
    double a1,b1;
    a1=fabs(a);
    b1=fabs(b);
    if(a1>b1){
        return a1*sqrt(1.0+SQR(b1/a1));
    }else{
        return (b1==0.0? 0.0 : b1*sqrt(1.0+SQR(a1/b1)));
    }
}

void QR_decompose_rotation(double *a,double *q,double *r,int n){
    int i,j,k;
    int idxi,idxj;
    double c,s,sr,ri,rj,qi,qj,tmp;
    for(i=0;i<n*n;i++){
        r[i]=a[i];
        q[i]=0;
    }
    for(i=0;i<n;i++){
        q[i*n+i]=1.0;
    }
    for(j=0;j<n;j++){
        for(i=j+1;i<n;i++){
            sr=pythag_m(r[j*n+j],r[i*n+j]);
            if(sr==0.0){
                c=1.0;
                s=0.0;
            }else{
                c=r[j*n+j]/sr;
                s=r[i*n+j]/sr;
            }
            //Do rotation Tij
            for(k=0;k<n;k++){
                idxi=i*n+k;
                idxj=j*n+k;
                ri=r[idxi]*c-r[idxj]*s;
                rj=r[idxi]*s+r[idxj]*c;
                qi=q[idxi]*c-q[idxj]*s;
                qj=q[idxi]*s+q[idxj]*c;
                r[idxi]=ri;
                r[idxj]=rj;
                q[idxi]=qi;
                q[idxj]=qj;
            }
        }
    }
    //Do final transpose
    for(i=0;i<n;i++){
        for(j=i;j<n;j++){
            tmp=q[i*n+j];
            q[i*n+j]=q[j*n+i];
            q[j*n+i]=tmp;
        }
    }
    return;
}

//Решатель СЛАУ
void QR_solve(int N,double *q,double *r,double *b,double *x){
    int i,j;
    double *bn;
    double sum;
    bn=(double *)malloc(N*sizeof(double));
    for(i=0;i<N;i++){
        bn[i]=0;
        for(j=0;j<N;j++){
            bn[i]+=q[j*N+i]*b[j];//bn=Qt*b;
        }
    }
    for(i=N-1;i>=0;i--){
        sum=0;
        for(j=i+1;j<N;j++){
            sum+=r[i*N+j]*x[j];
        }
        x[i]=(bn[i]-sum)/r[i*N+i];
    }
    free(bn);
    return;
}

void QR_decompose_reflection(double *a,double *q,double *r,int n){
    int i,j,k;
    double *v;
    double norm,dot1,dot2,tmp;
    v=(double *)malloc(n*sizeof(double));
    for(i=0;i<n*n;i++){
        r[i]=a[i];
        q[i]=0;
    }
    for(i=0;i<n;i++){
        q[i*n+i]=1.0;
    }
    for(j=0;j<n-1;j++){
        norm=0;
        for(i=0;i<n;i++){
            v[i]=0;
            if(i>=j){
                v[i]=r[i*n+j];
            }
            norm+=v[i]*v[i];
        }
        norm=sqrt(norm);
        if(v[j]>0){
            v[j]=v[j]+norm;
        }else{
            v[j]=v[j]-norm;
        }
        norm=0;
        for(i=0;i<n;i++){
            norm+=v[i]*v[i];
        }
        norm=sqrt(norm);
        if(norm!=0.0){
            for(i=0;i<n;i++){
                v[i]=v[i]/norm;
            }
        }
        for(k=0;k<n;k++){
            dot1=0;
            dot2=0;
            for(i=0;i<n;i++){
                dot1+=r[i*n+k]*v[i];
                dot2+=q[i*n+k]*v[i];
            }
            dot1=dot1*2.0;
            dot2=dot2*2.0;
            for(i=0;i<n;i++){
                r[i*n+k]=r[i*n+k]-v[i]*dot1;
                q[i*n+k]=q[i*n+k]-v[i]*dot2;
            }
        }
    }
    //Do final transpose
    for(i=0;i<n;i++){
        for(j=i;j<n;j++){
            tmp=q[i*n+j];
            q[i*n+j]=q[j*n+i];
            q[j*n+i]=tmp;
        }
    }
    free(v);
}

#ifdef CONFIG_CRYPTO
double random_double() {
    union {
        uint64_t i;
        unsigned char c[sizeof(uint64_t)];
    } u;

    if (!RAND_bytes(u.c, sizeof(u.c))) {
        fprintf(stderr, "Can't get random bytes!\n");
        exit(1);
    }
    return (u.i >> 11) * (1.0/9007199254740992.0);
}
#endif

void printMatrix(double *a,char * title,int m,int n){
    int i,j;
    if(title!=NULL){
        printf("%s\n",title);
    }
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            printf("%20.13lg ",a[i*n+j]);
        }
        printf("\n");
    }
}

void printMatrixZ(double_complex_t *a,char * title,int m,int n){
    int i,j;
    if(title!=NULL){
        printf("%s\n",title);
    }
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            printf(" (%20.13lg + I %20.13lg) ",creal(a[i*n+j]),cimag(a[i*n+j]));
        }
        printf("\n");
    }
}

//Realisation of Runge-Kutta method of 4-th order
//Integrand calling convention:
//F(neq,t,y[],rv[]);
int rk4_step(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t,double dt,void *serviceData){
    int i;
    double *k1,*k2,*k3,*k4;
    double *y1,*rv;
    //Allocating memory
    k1=(double *)malloc(neq*sizeof(double));
    k2=(double *)malloc(neq*sizeof(double));
    k3=(double *)malloc(neq*sizeof(double));
    k4=(double *)malloc(neq*sizeof(double));
    y1=(double *)malloc(neq*sizeof(double));
    rv=(double *)malloc(neq*sizeof(double));
    F(neq,t,y,rv,serviceData);
    for(i=0;i<neq;i++){
        k1[i]=dt*rv[i];
        y1[i]=y[i]+0.5*k1[i];
    }
    F(neq,t+0.5*dt,y1,rv,serviceData);
    for(i=0;i<neq;i++){
        k2[i]=dt*rv[i];
        y1[i]=y[i]+0.5*k2[i];
    }
    F(neq,t+0.5*dt,y1,rv,serviceData);
    for(i=0;i<neq;i++){
        k3[i]=dt*rv[i];
        y1[i]=y[i]+k3[i];
    }
    F(neq,t+dt,y1,rv,serviceData);
    for(i=0;i<neq;i++){
        k4[i]=dt*rv[i];
        y[i]=y[i]+(1.0/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
    }
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(y1);
    free(rv);
}

int rk4(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t0,double dt,int nsteps,void *serviceData){
    double t;
    int i;
    t=t0;
    for(i=0;i<nsteps;i++){
        rk4_step(F,neq,y,t,dt,serviceData);
        t+=dt;
    }
}

int rk5_step(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t,double dt,void *serviceData){
    int i;
    double *k1,*k2,*k3,*k4,*k5,*k6;
    double *y1,*rv;
    //Allocating memory
    k1=(double *)malloc(neq*sizeof(double));
    k2=(double *)malloc(neq*sizeof(double));
    k3=(double *)malloc(neq*sizeof(double));
    k4=(double *)malloc(neq*sizeof(double));
    k5=(double *)malloc(neq*sizeof(double));
    k6=(double *)malloc(neq*sizeof(double));
    y1=(double *)malloc(neq*sizeof(double));
    rv=(double *)malloc(neq*sizeof(double));
    F(neq,t,y,rv,serviceData);
    for(i=0;i<neq;i++){
        k1[i]=dt*rv[i];
        y1[i]=y[i]+0.25*k1[i];
    }
    F(neq,t+0.25*dt,y1,rv,serviceData);
    for(i=0;i<neq;i++){
        k2[i]=dt*rv[i];
        y1[i]=y[i]+(3.0/32.0)*k1[i]+(9.0/32.0)*k2[i];
    }
    F(neq,t+(3.0/8.0)*dt,y1,rv,serviceData);
    for(i=0;i<neq;i++){
        k3[i]=dt*rv[i];
        y1[i]=y[i]+(1932.0/2197.0)*k1[i]-(7200.0/2197.0)*k2[i]+(7296.0/2197.0)*k3[i];
    }
    F(neq,t+(12.0/13.0)*dt,y1,rv,serviceData);///i am here
    for(i=0;i<neq;i++){
        k4[i]=dt*rv[i];
        y1[i]=y[i]+(439.0/216.0)*k1[i]-8.0*k2[i]+(3680.0/513.0)*k3[i]-(845.0/4104.0)*k4[i];
    }
    F(neq,t+dt,y1,rv,serviceData);
    for(i=0;i<neq;i++){
        k5[i]=dt*rv[i];
        y1[i]=y[i]-(8.0/27.0)*k1[i]+2.0*k2[i]-(3544.0/2565.0)*k3[i]+(1859.0/4104.0)*k4[i]-(11.0/40.0)*k5[i];
    }
    F(neq,t+0.5*dt,y1,rv,serviceData);
    for(i=0;i<neq;i++){
        k6[i]=dt*rv[i];
        y[i]=y[i]+(16.0/135.0)*k1[i]+(6656.0/12825.0)*k3[i]+(28561.0/56430.0)*k4[i]-(9.0/50.0)*k5[i]+(2.0/55.0)*k6[i];
    }
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(k6);
    free(y1);
    free(rv);
}

int rk5(void (*F)(int,double,double [],double[],void *),int neq,double y[],double t0,double dt,int nsteps,void *serviceData){
    double t;
    int i;
    t=t0;
    for(i=0;i<nsteps;i++){
        rk5_step(F,neq,y,t,dt,serviceData);
        t+=dt;
    }
}

double ipow(double x,int k){
    int i;
    double res;
    res=1.0;
    if(k>0){
        for(i=0;i<k;i++){
            res*=x;
        }
    }
    if(k<0){
        for(i=0;i<-k;i++){
            res/=x;
        }
    }
    return res;
}

double ifact(int k){
    double rv;
    int i;
    rv=1.0;
    for(i=1;i<=k;i++){
        rv*=i;
    }
    return rv;
}

void print_complex(char * msg,double_complex_t z){
    printf("%s (%.16lg %.16lg)\n",msg,creal(z),cimag(z));
}

void print_complex_vector(char * msg,int n,double_complex_t *z){
    int i;
    printf("%s\n",msg);
    for(i=0;i<n;i++){
        printf("(%.16lg %.16lg)\n",creal(z[i]),cimag(z[i]));
    }
}

void print_complex_matrix(char * msg,int n,double_complex_t *z){
    int i,j;
    printf("%s\n",msg);
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("(%.16lg %.16lg) ",creal(z[i*n+j]),cimag(z[i*n+j]));
        }
        printf("\n");
    }
}


//Use the continued fraction representation
double_complex_t _zgamcontinued(double_complex_t z,double_complex_t a,double_complex_t b,double_complex_t c,double_complex_t d,int deep){
    double_complex_t term;
    if(deep>0){
        return a-b/(c+d/_zgamcontinued(z,a+2,b+z,c+2,d+z,deep-1));
    }else{
        return a;
    }
}


double_complex_t ZGammaIncomplete(double s,double_complex_t z){
    double_complex_t a,b,c,d;
    double_complex_t term;
    int deep;
    a=s;
    b=s*z;
    c=s+1;
    d=z;
    deep=150;
    return cpow(z,s)*cexp(-z)/_zgamcontinued(z,a,b,c,d,deep);
}

double_complex_t _zexpintegrale(double_complex_t z,double_complex_t a,double_complex_t b,double_complex_t c,double_complex_t d,int deep){
    double_complex_t term;
    if(deep>0){
        return a+b/(c+d/_zexpintegrale(z,a,b+1,c,d+1,deep-1));
    }else{
        return a;
    }
}


double_complex_t ZExpIntegralE(double nu,double_complex_t z){
    double_complex_t a,b,c,d;
    double_complex_t term;
    int deep;
    a=z;
    b=nu;
    c=1;
    d=1;
    deep=150;
    return cexp(-z)/_zexpintegrale(z,a,b,c,d,deep);
}

double min(double a,double b){
    if(a<b){
        return a;
    }
    return b;
}

double max(double a,double b){
    if(a<b){
        return b;
    }
    return a;
}

double ClebschGordan(double j1,double m1,double j2,double m2,double j3,double m3){
    double rval,a,b,n,xmin,xmax;
    int nmin,nmax,i,v;
    a=0;
    b=0;
    if((m1+m2!=m3)||(fabs(m1)>j1)||(fabs(m2)>j2)||(fabs(m3)>j3)){
        rval=0;
    }else{
        if((fabs(j1-j2)<=j3)&&(fabs(j1+j2)>=j3)){
            xmin=0;
            xmin=max(xmin,j2-j3-m1);
            xmin=max(xmin,j1-j3+m2);
            xmax=min(j2+m2,j1+j2-j3);
            xmax=min(xmax,j1-m1);
            nmin=round(xmin);
            nmax=round(xmax);
            a=sqrt((2.0*j3+1.0)*Gamma(j1+j2-j3+1)*Gamma(j3+j1-j2+1)*Gamma(j3+j2-j1+1)/Gamma(j1+j2+j3+2));
            b=0;
            for(i=nmin;i<=nmax;i++){
                n=i;
                b=b+(ipow(-1,i))*sqrt(Gamma(j1+m1+1)*Gamma(j1-m1+1)*Gamma(j2+m2+1)*Gamma(j2-m2+1)*Gamma(j3+m3+1)*Gamma(j3-m3+1))/(Gamma(n+1)*Gamma(j1+j2-j3-i+1)*Gamma(j1-m1-i+1)*Gamma(j2+m2-i+1)*Gamma(j3-j2+m1+i+1)*Gamma(j3-j1-m2+i+1));
            }
            rval=a*b;
        }else{
            rval=0;
        }
    }
    return rval;
}


//Implementation of Borwein Arc Bessel technique.
//TODO.

//Use external implementation


double dbesselj(double nu, double z){
    return creal(besselj(nu,z));
}

double dbessely(double nu, double z){
    return creal(bessely(nu,z));
}

double dbesseli(double nu, double z){
    return creal(besseli(nu,z));
}

double dbesselk(double nu, double z){
    return creal(besselk(nu,z));
}



double_complex_t besselj(double nu, double_complex_t z){
    int kode=1;
    int n=1;
    double zr=creal(z);
    double zi=cimag(z);
    int nz,ierr;
    double cyr[1],cyi[1];
    double_complex_t res;
    if(nu<0){
        return (2*nu+2)*besselj(nu+1,z)/z-besselj(nu+2,z);
    }
    zbesj_(&zr,&zi,&nu,&kode,&n,cyr,cyi,&nz,&ierr);
    if(ierr!=0){
        printf("error!\n");
    }
    res=cyr[0]+I*cyi[0];
    return res;
}

double_complex_t bessely(double nu, double_complex_t z){
    int kode=1;
    int n=1;
    double zr=creal(z);
    double zi=cimag(z);
    int nz,ierr;
    double cyr[1],cyi[1];
    double wrkr[1],wrki[1];
    double_complex_t res;
    if(nu<0){
        return (2*nu+2)*bessely(nu+1,z)/z-bessely(nu+2,z);
    }
    zbesy_(&zr,&zi,&nu,&kode,&n,cyr,cyi,&nz,wrkr,wrki,&ierr);
    if(ierr!=0){
        printf("error!\n");
    }
    res=cyr[0]+I*cyi[0];
    return res;
}

double_complex_t besseli(double nu, double_complex_t z){
    int kode=1;
    int n=1;
    double zr=creal(z);
    double zi=cimag(z);
    int nz,ierr;
    double cyr[1],cyi[1];
    double_complex_t res;
    if(nu<0){
        return (2*nu+2)*besseli(nu+1,z)/z+besseli(nu+2,z);
    }
    zbesi_(&zr,&zi,&nu,&kode,&n,cyr,cyi,&nz,&ierr);
    if(ierr!=0){
        printf("error!\n");
    }
    res=cyr[0]+I*cyi[0];
    return res;
}

double_complex_t besselk(double nu, double_complex_t z){
    int kode=1;
    int n=1;
    double zr=creal(z);
    double zi=cimag(z);
    int nz,ierr;
    double cyr[1],cyi[1];
    double_complex_t res;
    if(nu<0){
        return -(2*nu+2)*besselk(nu+1,z)/z+besselk(nu+2,z);
    }
    zbesk_(&zr,&zi,&nu,&kode,&n,cyr,cyi,&nz,&ierr);
    if(ierr!=0){
        printf("error!\n");
    }
    res=cyr[0]+I*cyi[0];
    return res;
}

//Spherical Bessel Functions
double sj(int l,double x){
    return sqrt(pi/(2*x))*dbesselj(l+0.5,x);
}
double sy(int l,double x){
    return sqrt(pi/(2*x))*dbessely(l+0.5,x);
}
double si(int l,double x){
    return sqrt(pi/(2*x))*dbesseli(l+0.5,x);
}
double sk(int l,double x){
    return sqrt(pi/(2*x))*dbesselk(l+0.5,x);
}


//Add new code
double Pochgammer(int n,double alpha){
    int i;
    double prod=1;
    for(i=0;i<n;i++){
        prod*=(alpha+i);
    }
    return prod;
}

double Hypergeometric1F1(double a,double c,double z){
    double sum=0;
    double eps=1e-16;
    double term=1;
    int i,maxiter;
    maxiter=100;
    for(i=0;i<maxiter && fabs(term) >fabs(sum)*eps;i++){
        term=ipow(z,i)*Pochgammer(i,a)/Pochgammer(i,c)/ifact(i);
        sum+=term;
    }
    return sum;
}

double EulerBeta(double x,double y){
    return Gamma(x)*Gamma(y)/Gamma(x+y);
}

double to_m(double *tp,void *sd){
    w_params *ptr;
    double t,lambda,mu,z;
    ptr=(w_params *)sd;
    t=*tp;
    lambda=ptr->lambda;
    mu=ptr->mu;
    z=ptr->z;
    return pow((1+t),mu-lambda-1.0/2.0)*pow((1-t),mu+lambda-1.0/2.0)*exp(z*t/2.0);
}

double WhittakerM(double lambda,double mu,double z){
    double result,a,b;
    w_params p;
    p.mu=mu;
    p.lambda=lambda;
    p.z=z;
    a=-1;
    b=1;
    result=pow(z,mu+1.0/2.0)/(pow(2.0,2*mu)*EulerBeta(mu+lambda+1.0/2.0,mu-lambda+1.0/2.0));
    result=result*GaussIntegrate(to_m,(void *)&p,1,&a,&b,1000);
    return result;
}


double_complex_t to_f2(double *z,void *sd){
    f2_params *p;
    double u,v;
    double alpha,beta1,beta2,gamma1,gamma2,x,y;
    double_complex_t result;
    double_complex_t ca,cb;
    p=(f2_params *)sd;
    alpha=p->alpha;
    beta1=p->beta1;
    beta2=p->beta2;
    gamma1=p->gamma1;
    gamma2=p->gamma2;
    x=p->x;
    y=p->y;
    u=z[0];
    v=z[1];
    result=cpow(u,beta1-1.0)*cpow(v,beta2-1.0)*cpow(1-u,gamma1-beta1-1)*cpow(1-v,gamma2-beta2-1)*cpow(1-u*x-v*y,-alpha);
    if(IsInf(creal(result))){
        result=0;
    }
    if(IsNaN(creal(cpow(1-u*x-v*y,-alpha)))){
        result=cpow(u,beta1-1.0)*cpow(v,beta2-1.0)*cpow(1-u,gamma1-beta1-1)*cpow(1-v,gamma2-beta2-1);//Hack;; 0^0=1;
    }
    return result;
}

double_complex_t HypergeometricF2(double alpha,double beta1,double beta2,double gamma1,double gamma2,double x,double y){
    f2_params p;
    double_complex_t result;
    double a[2];
    double b[2];
    p.alpha=alpha;
    p.beta1=beta1;
    p.beta2=beta2;
    p.gamma1=gamma1;
    p.gamma2=gamma2;
    p.x=x;
    p.y=y;
    a[0]=0.0;
    b[0]=1.0;
    a[1]=0.0;
    b[1]=1.0;
    result=Gamma(gamma1)*Gamma(gamma2)/(Gamma(beta1)*Gamma(beta2)*Gamma(gamma1-beta1)*Gamma(gamma2-beta2));
    result*=ZGaussIntegrate(to_f2,(void *)&p,2,a,b,10);
    return result;
}

//простой ряд с ограниченной сходимостью
double HypergeometricFA(int n,double alpha[],double beta [],double gamma [],double z[]){
    int *mi;
    int i,carry,m_max,ms;
    double prod,sum;
    mi=(int *)malloc(n*sizeof(int));
    for(i=0;i<n;i++){
        mi[i]=0;
    }
    carry=0;
    m_max=100;
    sum=0;
    while(!carry){
        prod=1;
        ms=0;
        for(i=0;i<n;i++){
            prod*=Pochgammer(mi[i],beta[i])*ipow(z[i],mi[i])/(Pochgammer(mi[i],gamma[i])*ifact(mi[i]));
            ms+=mi[i];
        }
        prod*=Pochgammer(ms,alpha[0]);
        sum+=prod;
        carry=increment(mi,n,m_max);
    }
    return sum;
}

void swap_double(double &a,double &b){
    double tmp;
    tmp=b;
    b=a;
    a=tmp;
}

void swap_vector_double(std::vector<double> &a,std::vector<double> &b){
    size_t size,i;
    double tmp;
    assert(a.size()==b.size());
    size=a.size();
    for(i=0;i<size;i++){
        tmp=b[i];
        b[i]=a[i];
        a[i]=tmp;
    }
}

void swap_vector_complex(std::vector<double_complex_t> &a,std::vector<double_complex_t> &b){
    size_t size,i;
    double_complex_t tmp;
    assert(a.size()==b.size());
    size=a.size();
    for(i=0;i<size;i++){
        tmp=b[i];
        b[i]=a[i];
        a[i]=tmp;
    }
}


//Original version of EigenSystem;
int EigenSystemOrig(int n,double *a,double_complex_t *w,double_complex_t *zc,int eigen_vectors){
    double *wr,*wi,*z;
    int i,j;
    wr=(double *)malloc(n*sizeof(double));
    wi=(double *)malloc(n*sizeof(double));
    z=(double *)malloc(n*n*sizeof(double));
    rg_elm(n,a,wr,wi,eigen_vectors,z);
    for(i=0;i<n;i++){
        w[i]=wr[i]+I*wi[i];
        if(eigen_vectors){
            for(j=0;j<n;j++){
                if(wi[i]>0){
                    zc[j*n+i]=z[j*n+i]+I*z[j*n+i+1];
                    zc[j*n+i+1]=z[j*n+i]-I*z[j*n+i+1];
                }
                if(wi[i]==0){
                    zc[j*n+i]=z[j*n+i];
                }//wi<0 --skip
            }
        }
    }
    free(wr);
    free(wi);
    free(z);
}


int EigenSystem(int n,double *a,std::vector<double_complex_t> &w,std::vector<std::vector<double_complex_t>> &zc,int eigen_vectors,int nproc){
    Eigen::setNbThreads(nproc);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M;
    M.resize(n,n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            M(i,j)=a[i*n+j];
        }
    }
    Eigen::EigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> es(M);
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> ev=es.eigenvalues();
    w.clear();
    zc.clear();
    w.resize(n);
    std::vector<double_complex_t> cur;
    for(int i=0;i<n;i++){
        w[i]=ev(i,0).real()+I*ev(i,0).imag();
    }
    if(eigen_vectors){
        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> evec=es.eigenvectors();
        for(int j=0;j<n;j++){
            cur.clear();
            cur.resize(n);
            for(int i=0;i<n;i++){
                cur[i]=evec(i,j).real()+I*evec(i,j).imag();
            }
            zc.push_back(cur);
        }
    }
}

double AiryAi(double x,int derivative){
    double zr,zi;
    double fr,fi,zero=0.0;
    int id,kode,err,nz;
    id=(derivative!=0?1:0);
    kode=1;
    zr=x;
    zi=0;
    zairy_(&zr,&zi,&id,&kode,&fr,&fi,&nz,&err);
    if(err){
        printf("ERROR\n");
    }
    return fr;
}

double AiryBi(double x,int derivative){
    double zr,zi;
    double fr,fi,zero=0.0;
    int id,kode,err;
    id=(derivative!=0?1:0);
    kode=1;
    zr=x;
    zi=0;
    zbiry_(&zr,&zi,&id,&kode,&fr,&fi,&err);
    if(err){
        printf("ERROR\n");
    }
    return fr;
}

double ai(double x){
    return AiryAi(x,0);
}

double aiprime(double x){
    return AiryAi(x,1);
}
double bi(double x){
    return AiryBi(x,0);
}

double biprime(double x){
    return AiryBi(x,1);
}

double Zeta(double s){
    int n,k,sign;
    double result=0;
    double addition;
    for(n=0;n<100;n++){
        addition=0;
        for(k=0;k<=n;k++){
            sign=(k%2==0?1:-1);
            addition+=(1.0/(pow(2,n+1)))*BinomCoeff(n,k)*sign/(pow(k+1,s));
        }
        result+=addition;
        if(fabs(addition/result)<1e-16){
            break;
        }
    }
    result=result/(1-pow(2,1-s));
    return result;
}

double PolyLog(double s,double z){
    double result,add,fact=1.0;
    int n;
    result=Gamma(1-s)*pow(log(1.0/z),s-1);
    for(n=0;n<100;n++){
        if(n>=1){
            fact*=n;
        }
        add=Zeta(s-n)*pow(log(z),n)/fact;
        result+=add;
        if(fabs(add/result)<1e-17){//Convergence
            break;
        }
    }
    return result;
}


/**
 * @brief find_root
 * Процедура поиска корней системы нелинейных уравнений.
 * @param f
 * Тестируемая функция
 * @param x_init
 * Вектор начального приближения (размер -- число неизвестных.)
 * @param x_out
 * Вывод: решение системы
 * @param user_data
 * Пользовательские данные, которые передаются в функцию.
 * @param epsilon
 * Задаваемая погрешность
 * @param max_iter
 * Максимальное число итераций.
 * @return
 */

bool find_root(std::function<std::vector<double>(std::vector<double> &x,void *user_data)> f, const std::vector<double> &x_init, std::vector<double> &x_out, void *user_data, double epsilon, int max_iter)
{
    double dx = 0, f_norm = 0, x_norm = 0;
    double x_min = 0, x_max = 0;
    int iter = 0;
    bool flag = true;
    size_t num_dimensions;
    std::vector<double> x, x_shift, f_out;
    num_dimensions = x_init.size();
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> V;
    M.resize(num_dimensions, num_dimensions);
    V.resize(num_dimensions, 1);
    x.resize(num_dimensions);
    x_shift.resize(num_dimensions);
    f_out.resize(num_dimensions);
    for(int i = 0; i < num_dimensions; i ++)
    {
        x_shift[i] = x_init[i];
        x[i] = x_init[i];
    }
    //j -- function
    //i -- xi;
    //M(j,i);
    iter  = 0;
    x_min = fabs(x[0]);
    x_max = 0;
    for(int i = 0;i < num_dimensions;i ++)
    {
        x_min = (fabs(x[i]) < x_min ? fabs(x[i]) : x_min);
        x_max = (fabs(x[i]) > x_max ? fabs(x[i]) : x_max);
    }
    double x_mean = (x_min + x_max)/2.0;
    dx = (x_mean > 0 ? x_mean * epsilon : 0);
    do
    {
        for(int i = 0;i < num_dimensions;i ++)
        {
            x_shift[i] = x[i];
        }
        for(int i = 0;i < num_dimensions;i ++)
        {
            x_shift[i] = x[i] + dx;
            for(int j = 0;j < num_dimensions;j ++)
            {
                M(j,i) = (f(x_shift, user_data)[j] - f(x,user_data)[j])/dx;
            }
            x_shift[i] = x[i];
        }
        M = M.inverse();
        for(int i = 0; i < num_dimensions; i ++)
        {
            V(i,0) = f(x,user_data)[i];
        }
        V = M * V;
        for(int i = 0; i < num_dimensions; i ++)
        {
            x[i]=x[i]-V(i,0);
        }
        f_out=f(x,user_data);
        f_norm=0;
        x_norm=0;
        for(int i = 0; i < num_dimensions; i ++)
        {
            f_norm += f_out[i] * f_out[i];
        }
        f_norm = sqrt(f_norm);
        x_norm = V.norm();
        if(f_norm < epsilon || x_norm < epsilon)
        {
            flag = false;//STOP;
        }
        iter++;
        if(iter >= max_iter)
        {
            //Convergence doesn't achieved
            x_out = x;
            return false;
        }
    }while(flag);
    x_out = x;
    return true;
}

bool find_minimum_1d(std::function<double (double &, void *)> f, const double &x_init, double &x_out, void *user_data, double epsilon, int max_iter)
{
    double step;
    double a,b,c,x,y;//Концы интервала.
    double p,q,alpha;
    double s,fy,fx,fa,fb,fc,delta;
    double l,r;
    bool flag,direct;
    int cnt;
    step=epsilon;
//    step=(fabs(x_init)>0? x_init*epsilon : epsilon);
    flag=true;
    //Итерации пока не определим интервал.
    //Далее -- автоматический выбор шага.
    l=x_init-step;
    r=x_init+step;
    while(f(l,user_data)==f(r,user_data))
    {
        step*=2;
        l=x_init-step;
        r=x_init+step;
    }
    step*=2.0;
    step*=5;
    c=x_init;
    b=c+step;
    fb=f(b,user_data);
    fc=f(c,user_data);
    cnt=0;
    direct=(fb<fc);
    if(direct){//DIRECT
        while(fb<fc){
            a=c;
            fa=fc;
            c=b;
            fc=fb;
            step*=2;
            b=c+step;
            fb=f(b,user_data);
        }
    }else{//INVERSE
        a=c-step;
        fa=f(a,user_data);
        while(fa<fc){
            b=c;
            fb=fc;
            c=a;
            fc=fa;
            step*=2;
            a=c-step;
            fa=f(a,user_data);
        }
    }
    //Есть интервал.
    if(b<a){
        c=a;
        a=b;
        b=c;
    }
    //Используем метод аппроксимации параболами.
    flag=true;
    x=(a+b)/2.0;
    fa=f(a,user_data);
    fb=f(b,user_data);
    fx=f(x,user_data);
    cnt=0;
    do{
        cnt++;
        p=(x-a)*(fb-fx);
        q=(b-x)*(fa-fx);
        s=(p+q);
        y=(p*(x+a)+q*(b+x))/(2*s);
        alpha=-s/((x*x-a*a)*(b-x)-(b*b-x*x)*(x-a));
        delta=fabs(x-y);
        if((y<=a) || (y>=b) || std::isnan(y) || std::isinf(y)){
            break;
        }
        fy=f(y,user_data);
        if(fy>fx){
            if(y<x){
                a=y;
                fa=fy;
            }else{
                b=y;
                fb=fy;
            }
        }else{
            if(y<x){
                b=x;
                fb=fx;
                x=y;
                fx=fy;
            }else{
                a=x;
                fa=fx;
                x=y;
                fx=fy;
            }
        }
        if(delta<epsilon){
            flag=false;
        }
        if(cnt>max_iter){
            return false;
        }
    }while(flag);
    if(fa<fx){
        x=a;
        fx=fa;
    }
    if(fb<fx){
        x=b;
        fx=fa;
    }
    x_out=x;
    return true;
}

/**
 * @brief find_minimum
 * Процедура поиска минимума методом градиентного спуска.
 * @param f
 * @param x_init
 * @param x_out
 * @param user_data
 * @param epsilon
 * @param max_iter
 * @return
 */

bool find_minimum_descent(std::function<double (std::vector<double> &, void *)> f, const std::vector<double> &x_init, std::vector<double> &x_out, void *user_data,double lambda_init, double epsilon, int max_iter)
{
    double dx = 0, f_norm = 0, x_norm = 0;
    double x_min = 0, x_max = 0;
    int iter = 0;
    bool flag = true;
    size_t num_dimensions;
    std::vector<double> x, x_shift, x_tmp, f_grad;
    num_dimensions = x_init.size();
    x.resize(num_dimensions);
    x_shift.resize(num_dimensions);
    x_tmp.resize(num_dimensions);
    f_grad.resize(num_dimensions);
    for(int i = 0; i < num_dimensions; i ++)
    {
        x_shift[i] = x_init[i];
        x[i] = x_init[i];
    }
    x_min = fabs(x[0]);
    x_max = 0;
    for(int i = 0;i < num_dimensions;i ++)
    {
        x_min = (fabs(x[i]) < x_min ? fabs(x[i]) : x_min);
        x_max = (fabs(x[i]) > x_max ? fabs(x[i]) : x_max);
    }
    double x_mean = (x_min + x_max)/2.0;
    dx = (x_mean > 0? x_mean * epsilon : epsilon);//Выбираем шаг.
    double lambda=lambda_init;//Начальное значение параметра обучения.
    do{
        iter++;
        for(int i = 0;i < num_dimensions;i ++)
        {
            x_shift[i] = x[i];
        }
        for(int i = 0;i < num_dimensions;i ++)
        {
            x_shift[i] = x[i] + dx;
            f_grad [i] = (f(x_shift, user_data) - f(x,user_data))/dx;
            x_shift[i] = x[i];
        }
        auto caller_function =
                [x,dx,f_grad,num_dimensions,f]  (double &lambda,void *user_data)
        {
            double result;
            std::vector<double> x_tmp;
            x_tmp.resize(num_dimensions);
            for(int i=0;i<num_dimensions;i++){
                x_tmp[i]=x[i]-lambda*f_grad[i];
            }
            result=f(x_tmp,user_data);
            return result;
        };
        double lambda_out;
        lambda=lambda_init;
        if(!find_minimum_1d(caller_function,lambda,lambda_out,user_data)){
            //Convergence doesn't achievied;
            return false;
        }
        lambda=lambda_out;
        printf("iter= %d lambda= %lg\n",iter,lambda_out);
        //        lambda_out=0.01;
        for(int i=0;i<num_dimensions;i++){
            x_tmp[i]=x[i]-lambda*f_grad[i];
        }
        x_norm=0;
        for(int i=0;i<num_dimensions;i++){
            x_norm+=(x_tmp[i]-x[i])*(x_tmp[i]-x[i]);
        }
        for(int i=0;i<num_dimensions;i++){
            x[i]=x_tmp[i];
        }
        x_norm=sqrt(x_norm);
        f_norm=fabs(f(x,user_data)-f(x_tmp,user_data));
        double grad_norm=0;
        for(int i=0;i<num_dimensions;i++){
            grad_norm+=f_grad[i]*f_grad[i];
        }
        grad_norm=sqrt(grad_norm);
        if((x_norm<epsilon || f_norm <epsilon)&&(grad_norm<epsilon)){
            flag=false;
        }
        if(iter>max_iter){
            return false;//Convergence doesn't achieved.
        }
    }while(flag);
    x_out=x;
    return true;
}
