/******************************************************
This library allows us to manipulate periodic functions
by approximating them using discrete Fourier transform (DFT).

The class "grid" stores either a discrete (equispaced) sampling
of points on the standard torus or the corresponding DFT.

We notice that the real object that approximate the function
is the corresponding truncated Fourier series rather than the
interpolating trigonometric polynomial. This is important in
order to approximate the derivatives of the funcion.

Some elementary operators (=,+,-,*,/) have been overloaded
in order to simplify the coding of the parameterization
method.

IMPORTANT REMARK: Notice that the same class is used to
store both the sampled points of the funcion and the Fourier coefficients.
Notice also that most of the routines below (deriva, cohomological, etc)
act on the Fourier cofficients. This has de advantage of simplifying
the implementation of the class. Be careful!
******************************************************/
#ifndef GRID_H
#define GRID_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include "complex.h"

using namespace std;

class grid
{
    public:

    int ndim;      // Number of angular variables
    int nelem;     // Total number of points in the grid
    int *nn;       // Total number of points for each angle
    complex *elem; // Sampled points or Fourier coefficients

    grid();
    grid(const grid& y);
    grid(int n, int *index);
    ~grid();
    
    grid& operator = (const grid& y);

    void alloc_grid(){
        nn=new int[ndim];
        elem=new complex[nelem];
    }
};

grid operator - (const grid& g);
grid operator + (const grid& g1, const grid& g2);
grid operator - (const grid& g1, const grid& g2);
grid operator * (const grid& g1, const grid& g2);
grid operator + (const grid& g, const double& x);
grid operator * (const grid& g, const double& x);

void start_grid(const double& tol);
void four1(double *data, int nn, int isign);
void fourn(double *data, int *nn, int ndim, int isign);
grid fft_F(const grid& g);
grid fft_B(const grid& g);
void trigo_to_series(int *nn, int *index, int *indexserie, int ndim);
void series_to_trigo(int *nn, int *indexserie, int *index, int ndim);
int position(int *nn, int *index, int ndim);
void indices(int pos, int *nn, int *index, int ndim);
grid deriva (const grid& g, int k);
grid shift (const grid& g, double *omega);
grid shiftc (const grid& g, const double &rho);
grid cohomological (const grid& g, double *omega);
void zerogrid(grid& g);
void tail(const grid& g, double *tg);
void clean(grid& g);
void aver(const grid& g, double &avg);
void print_gridR(const grid& g);
void print_gridF(const grid& g);
double norm(const grid& g);
double normsobo(const grid& g, const double& r);
double normsup(const grid& g, const double& rho);
double normexp(const grid& g, const double& rho);
int compare (const void * a, const void * b);

double tolgrid;

grid::grid()
{
    ndim=0;
    nelem=0;
    nn=NULL;
    elem=NULL;
}

grid::grid(int n, int *index)
{
    ndim=n;
    nelem=1;
    for(int i=0;i<ndim;i++) nelem = nelem*index[i];

    alloc_grid();
    for(int i=0;i<ndim;i++) nn[i]=index[i];
    for(int i=0;i<nelem;i++) elem[i]=0.0;
}

grid::grid(const grid& y)
{
    ndim=y.ndim;
    nelem=y.nelem;

    alloc_grid();
    for(int i=0;i<ndim;i++) nn[i]=y.nn[i];
    for(int i=0;i<nelem;i++) elem[i]=y.elem[i];
}

grid::~grid()
{
    if(nn!=NULL){
        delete [] nn;
    }
    if(elem!=NULL){
        delete [] elem;
    }
}

inline grid& grid::operator = (const grid& y)
{
    if (this == &y){
        return *this;
    }
    else {
        if (ndim!=y.ndim || nelem!=y.nelem){
            this->~grid();
            ndim=y.ndim;
            nelem=y.nelem;
            alloc_grid();
        }
        for(int i=0;i<ndim;i++) nn[i]=y.nn[i];
        for(int i=0;i<nelem;i++) elem[i]=y.elem[i];
        return *this;
    }
}

grid operator + (const grid& g1, const grid& g2)
{
    if (g1.ndim!=g2.ndim || g1.nelem!=g2.nelem){
        cout << "grid operator +: error in dimensions" << endl;
    }

    grid suma (g1.ndim,g1.nn);
    for(int i=0;i<g2.nelem;i++) suma.elem[i]=g1.elem[i]+g2.elem[i];
    return suma;
}

grid operator - (const grid& g)
{
    grid menos (g.ndim,g.nn);
    for(int i=0;i<g.nelem;i++) menos.elem[i]=-g.elem[i];
    return menos;
}

grid operator - (const grid& g1, const grid& g2)
{
    if (g1.ndim!=g2.ndim || g1.nelem!=g2.nelem){
        cout << "grid operator -: error in dimensions" << endl;
    }

    grid resta (g1.ndim,g1.nn);
    for(int i=0;i<g2.nelem;i++) resta.elem[i]=g1.elem[i]-g2.elem[i];
    return resta;
}

grid operator * (const grid& g1, const grid& g2)
{
    if (g1.ndim!=g2.ndim || g1.nelem!=g2.nelem){
        cout << "grid operator *: error in dimensions" << endl;
    }

    grid prod (g1.ndim,g1.nn);
    for(int i=0;i<g2.nelem;i++) prod.elem[i]=g1.elem[i]*g2.elem[i];
    return prod;
}

grid operator + (const grid& g, const double& x)
{
    grid suma (g.ndim,g.nn);
    for(int i=0;i<g.nelem;i++) suma.elem[i]=g.elem[i]+x;
    return suma;
}

grid operator * (const grid& g, const double& x)
{
    grid prod (g.ndim,g.nn);
    for(int i=0;i<g.nelem;i++) prod.elem[i]=g.elem[i]*x;
    return prod;
}

void start_grid(const double& tol)
/*
    This defines the global variable tolgrid that we use in
    the solution of the 1-bite cohomnological equation
    and in the clean function
*/
{
    tolgrid=tol;
}

grid fft_F(const grid& g)
/*
    Forward discrete Fourier transform using the FFT
    routines in:
    "WH Press, SA Teukolsky, WT Vetterling, BP Flannery.
    Numerical Recipes in C: the art of scientific computing"

    We can use also FFTW3 rourtines, which are faster.
    However, the routines four1 and fourn have the advantage
    that they can be extended to higher precision arithmetics
    in a straightforward way.
*/
{
    grid g2(g.ndim,g.nn);
    double *data, tot;

    data = new double [2*g.nelem];
    for (int i=0;i<g.nelem;i++){
        data[2*i]   = g.elem[i].real;
        data[2*i+1] = g.elem[i].imag;
    }
    if (g.ndim==1) four1(data,g.nelem,-1);
    else fourn(data,g.nn,g.ndim,-1);
    tot=g.nelem;
    for (int i=0;i<g.nelem;i++){
        g2.elem[i].real = data[2*i]/tot;
        g2.elem[i].imag = data[2*i+1]/tot;
    }

    delete [] data;
    return g2;
}

grid fft_B(const grid& g)
/*
    Backward discrete Fourier transform using FFT.
    See comments in
    grid fft_F(const grid& g)
*/
{
    grid g2(g.ndim,g.nn);
    double *data;

    data = new double [2*g.nelem];
    for (int i=0;i<g.nelem;i++){
        data[2*i]   = g.elem[i].real;
        data[2*i+1] = g.elem[i].imag;
    }
    if (g.ndim==1) four1(data,g.nelem,1);
    else fourn(data,g.nn,g.ndim,1);
    for (int i=0;i<g.nelem;i++){
        g2.elem[i].real = data[2*i];
        g2.elem[i].imag = data[2*i+1];
    }

    delete [] data;
    return g2;
}

void tail(const grid& g, double *tg)
/*  Given a DFT stored in the grid g, this routines
    computes the "averaged" tails defined as the sum
    of the Fourier coefficients on the sets
    nn[i]/4 \leq index[i] \leq 3*nn[i]/4

    The reason to average is to be able to consider
    tolerances independent of the size of the grid.
*/
{
    int *index;
    index=new int[g.ndim];

    complex **copyg;
    copyg=new complex*[g.ndim];
    for(int i=0;i<g.ndim;i++) copyg[i] =new complex[g.nelem];
    for(int l=0;l<g.nelem;l++){
    indices(l,g.nn,index,g.ndim);
    for(int i=0;i<g.ndim;i++){
        if (g.nn[i]<=4*index[i] && 4*index[i]<=3*g.nn[i])
               copyg[i][l]=g.elem[l];
    }
    }
    for(int i=0;i<g.ndim;i++) {
        tg[i]=0.0;
        qsort (copyg[i], g.nelem, sizeof(complex), compare);
        for(int l=0;l<g.nelem;l++){
            tg[i]=tg[i]+abs(copyg[i][l]);
        }    
    }
    for(int i=0;i<g.ndim;i++) tg[i] = 2.0*tg[i]/g.nelem;
    for(int i=0;i<g.ndim;i++) delete [] copyg[i];
    delete [] copyg;

    delete [] index;
}

void clean(grid& g)
/*
   Given a grid of Fourier coefficients of a ndim-periodic function
   we set to 0.0 those coefficients that are below to tolerance tolgrid
*/
{
    int *index;

    index=new int[g.ndim];

    for(int l=0;l<g.nelem;l++){ 
        if (abs(g.elem[l])<tolgrid) g.elem[l]=0.0;
    }
    delete [] index;

}

int compare (const void * a, const void * b)
{
    if (abs(*(const complex*)a) < abs(*(const complex*)b)) return -1;
    return abs(*(const complex*)a) > abs(*(const complex*)b);
}

double norm(const grid& g)
/*
   Given a grid of Fourier coefficients of a ndim-periodic function
   we compute the l1 norm
*/
{
    double normg;
    complex *copyg;

    copyg=new complex[g.nelem];
    for(int l=0;l<g.nelem;l++) copyg[l]=g.elem[l];
    qsort(copyg, g.nelem, sizeof(complex), compare);

    normg=0.0;
    for(int l=0;l<g.nelem;l++){
        normg=normg+abs(copyg[l]);
    }

    delete [] copyg;
    return normg;
}

double normsup(const grid& g, const double& rho)
/*
   Given a grid of Fourier coefficients of a ndim-periodic function
   we compute the supremum norm by approximating the maximum
   of the funcion at the boundary
*/
{
    double normg;
    grid gR;

    gR=shiftc(g,rho); // Shift to the boundary
    gR=fft_B(gR);     // Evaluation at the boundary

    normg=0.0;
    for(int l=0;l<gR.nelem;l++){
        if (normg<abs(gR.elem[l])) normg=abs(gR.elem[l]);
    }
    return normg;
}

double normexp(const grid& g, const double& rho)
/*
   Given a grid of Fourier coefficients of a ndim-periodic function
   we compute the "Fourier norm", i.e., the l1-norm of the
   weighted Fouier coefficients.
*/
{
    double norm, aux1, aux2, pi2;
    int *index;
    int *index_serie;

    index=new int[g.ndim];
    index_serie=new int[g.ndim];

    norm=0.0;
    aux1=0.0;
    aux2=0.0;
    pi2=6.2831853071795864769252867665590057684;

    for(int l=0;l<g.nelem;l++){
        indices(l,g.nn,index,g.ndim);
        trigo_to_series(g.nn,index,index_serie,g.ndim);
        aux1=abs(g.elem[l]);
        aux2=0.0;
        for(int i=0;i<g.ndim;i++){
            if(index_serie[i]>0) aux2 = aux2 + index_serie[i];
            if(index_serie[i]<0) aux2 = aux2 - index_serie[i];
        }
        norm = norm + aux1*exp(pi2*aux2*rho);
    }

    delete [] index;
    delete [] index_serie;

    return norm;

}

double normsobo(const grid& g, const double& r)
/*
   Given a grid of Fourier coefficients of a ndim-periodic function
   we compute the r-Sobolev semi-norm.
*/
{
    double norm, aux1, aux2;
    int *index;
    int *index_serie;

    index=new int[g.ndim];
    index_serie=new int[g.ndim];

    norm=0.0;
    aux1=0.0;
    aux2=0.0;

    for(int l=0;l<g.nelem;l++){
    indices(l,g.nn,index,g.ndim);
    trigo_to_series(g.nn,index,index_serie,g.ndim);
    aux1=abs(g.elem[l])*abs(g.elem[l]);
    aux2=0.0;
    for(int i=0;i<g.ndim;i++){
        if(index_serie[i]>0) aux2 = aux2 + index_serie[i];
        if(index_serie[i]<0) aux2 = aux2 - index_serie[i];
    }
    aux1 = pow(aux2,2.0*r)*aux1;
        norm=norm+aux1;
    }
    norm=sqrt(norm);

    delete [] index;
    delete [] index_serie;

    return norm;
}

void aver(const grid& g, double &avg)
/*
   Given a grid of sampled points of a ndim-periodic function
   we compute the average.
*/
{
    double aux;
    avg=0.0;
    for(int l=0;l<g.nelem;l++){
        avg = avg + g.elem[l].real;
    }
    aux=g.nelem;
    avg = avg/aux;
}

void four1(double *data, int nn, int isign)
/*
   THIS ROUTINE IS PROVIDED IN "Numerical Recipes in C". Cambridge

   Replaces data[0,...2*nn-1] by its discrete Fourier transform if isign=1
                              by its inverse "    "       "     if isign=-1
   nn MUST be an integer power of 2
   *data has the format data[0]=Re(point0), data[1]=Im(point0)
                        data[2]=Re(point1), data[3]=Im(point1)
                        data[4]=Re(point2), data[5]=Im(point2)
                        ...
                        data[2*j]=Re(pointj), data[2*j+1]=Im(pointj)
 */
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi, aux;
    double val0, val1, val2, pi2;

    val0=0.0;
    val1=1.0;
    val2=2.0;
    pi2=6.2831853071795864769252867665590057684;
    
    n = nn << 1;
    j = 1;
    
    for(i = 1; i < n; i += 2){
    if(j > i){
        aux = data[j-1];
        data[j-1] = data[i-1];
        data[i-1] = aux;
        aux = data[j];
        data[j] = data[i];
        data[i] = aux;
    }
    m = n >>1;
    while(m >= 2 && j > m){
        j -= m;
        m >>= 1;
    }
        j += m;
    }
    mmax = 2;
    while(n > mmax){
        istep = mmax << 1;
        theta = isign*(pi2/(int)mmax);
        wtemp = sin(theta/val2);
        wpr = val0-val2*wtemp*wtemp;
        wpi = sin(theta);
        wr = val1;
        wi = val0;
        for(m = 1; m < mmax; m += 2){
            for(i = m; i <= n; i += istep){
                j = i+mmax;
                tempr = wr*data[j-1]-wi*data[j];
                tempi = wr*data[j]+wi*data[j-1];
                data[j-1] = data[i-1]-tempr;
                data[j] = data[i]-tempi;
                data[i-1] = data[i-1]+tempr;
                data[i] = data[i]+tempi;
            }
            wtemp = wr;
            wr = wtemp*wpr-wi*wpi+wr;
            wi = wi*wpr+wtemp*wpi+wi;
        }
        mmax = istep;
    }
    return ;
}

void fourn(double *data, int *nn, int ndim, int isign)
/*
    THIS ROUTINE IS PROVIDED IN "Numerical Recipes in C". Cambridge

    Replaces data[0,...2*L-1] by its discrete Fourier transform if isign=1
                              by its inverse "    "       "     if isign=-1
    The inter L that defines the lengh of the array *data is L=nn[0]*nn[1]*...*nn[ndim-1]
    nn[0..ndim-1] is an integer array containing the lenghs of each dimension,
                  which MUST be all integer powers of 2
    *data has the format of row-major order using that two consecutive numbers
          define a complex number.
*/
{
    int idim;
    int i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
    int ibit,k1,k2,n,nprev,nrem,ntot;
    double tempi,tempr,aux;
    double theta,wi,wpi,wpr,wr,wtemp;
    double pi2;
    
    pi2=6.2831853071795864769252867665590057684;
    
    for (ntot=1,idim=1;idim<=ndim;idim++)
    ntot *= nn[idim-1];
    nprev=1;
    for (idim=ndim;idim>=1;idim--) {
        n=nn[idim-1];
        nrem=ntot/(n*nprev);
        ip1=nprev << 1;
        ip2=ip1*n;
        ip3=ip2*nrem;
        i2rev=1;
        for (i2=1;i2<=ip2;i2+=ip1) {
            if (i2 < i2rev) {
                for (i1=i2;i1<=i2+ip1-2;i1+=2) {
                    for (i3=i1;i3<=ip3;i3+=ip2) {
                        i3rev=i2rev+i3-i2;
                        aux=data[i3-1];
                        data[i3-1]=data[i3rev-1];
                        data[i3rev-1]=aux;

                        aux=data[i3];
                        data[i3]=data[i3rev];
                        data[i3rev]=aux;
                    }
                }
            }
            ibit=ip2 >> 1;
            while (ibit >= ip1 && i2rev > ibit) {
                i2rev -= ibit;
                ibit >>= 1;
            }
            i2rev += ibit;

        }
        ifp1=ip1;
        while (ifp1 < ip2) {
            ifp2=ifp1 << 1;
            theta=isign*pi2/(ifp2/ip1);
            wtemp=sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi=sin(theta);
            wr=1.0;
            wi=0.0;
            for (i3=1;i3<=ifp1;i3+=ip1) {
                for (i1=i3;i1<=i3+ip1-2;i1+=2) {
                    for (i2=i1;i2<=ip3;i2+=ifp2) {
                    k1=i2;
                    k2=k1+ifp1;
                    tempr=wr*data[k2-1]-wi*data[k2];
                    tempi=wr*data[k2]+wi*data[k2-1];
                    data[k2-1]=data[k1-1]-tempr;
                    data[k2]=data[k1]-tempi;
                    data[k1-1] = data[k1-1]+tempr;
                    data[k1] = data[k1]+tempi;
                    }
                }
                wtemp=wr;
                wr=wtemp*wpr-wi*wpi+wr;
                wi=wi*wpr+wtemp*wpi+wi;
            }
            ifp1=ifp2;
        }
        nprev *= n;
    }
}

void trigo_to_series(int *nn, int *index, int *indexserie, int ndim)
/*
    If *nn is the number of points in each angle and
    *index is a multiindex (in DFT) of the form "0<=index[i]<nn[i]" for
    "0<=i<ndim", then this routine compute the corresponding
    multiindex of the Fourier series.
*/
{
    for(int i=0;i<ndim;i++){
        if (index[i]==0) indexserie[i]=0;
        else if(index[i]<nn[i]/2) indexserie[i]=index[i];
        else indexserie[i]=index[i]-nn[i];
    }
}

void series_to_trigo(int *nn, int *indexserie, int *index, int ndim)
/*
    Inverse of the routine
    void trigo_to_series(int *nn, int *index, int *indexserie, int ndim)
*/
{
    for(int i=0;i<ndim;i++){
        if (indexserie[i]==0) index[i]=0;
        else if(indexserie[i]>0) index[i]=indexserie[i];
        else index[i]=indexserie[i]+nn[i];
    }
}

int position(int *nn, int *index, int ndim)
/*
    Storage of a coefficient. Given a multiindex *index
    related to the discretization *nn, this routine
    computes the position of the corresponding point or
    Fourier coefficient in the array *elem of the class
*/
{
    int pos, prod;
    pos=1;
    for(int k=1;k<=ndim;k++){
        prod=1;
        for(int l=k+1;l<=ndim;l++){
            prod = prod * nn[l-1];
        }
        pos = pos + prod*index[k-1];
    }
    return pos-1;
}

void indices(int pos, int *nn, int *index, int ndim)
/*
    Given a position "pos" of a point or Fourier coefficient 
    in the array *elem of the class, this routine computes
    the corresponding multiindex *index related to the
    discretization *nn
*/
{
    int pos0=pos;
    
    for(int k=ndim;k>1;k--){
        index[k-1] = pos0 % nn[k-1];
        pos0 = pos0/nn[k-1];
    }
    index[0]=pos0;
}

grid deriva (const grid& g, int k)
/*
   Given a grid of Fourier coefficients of a ndim-periodic function
   we compute the derivative with respect to theta_k
*/
{
    if (k>=g.ndim) cout << "Problem with index in routine deriva" << endl;
    int *index;
    int *index_serie;
    grid dg(g.ndim,g.nn);

    double pi2,aux;
    pi2=6.2831853071795864769252867665590057684;

    index=new int[g.ndim];
    index_serie=new int[g.ndim];

    for(int i=0;i<g.nelem;i++){
        indices(i,g.nn,index,g.ndim);
        trigo_to_series(g.nn,index,index_serie,g.ndim);
        aux=index_serie[k];
        dg.elem[i]=complex(0.0,pi2*aux)*g.elem[i];
    }

    delete [] index;    
    delete [] index_serie;

    return dg;
}

grid shift (const grid& g, double *omega)
/*
   Given a grid of Fourier coefficients of a ndim-periodic function
   we compute the composition f(theta+omega)
*/
{
    int *index;
    int *index_serie;
    grid shiftg(g.ndim,g.nn);

    double pi2,aux;
    pi2=6.2831853071795864769252867665590057684;

    index=new int[g.ndim];
    index_serie=new int[g.ndim];

    for(int i=0;i<g.nelem;i++){
        indices(i,g.nn,index,g.ndim);
        trigo_to_series(g.nn,index,index_serie,g.ndim);
        aux=0.0;
        for(int j=0;j<g.ndim;j++) aux = aux + index_serie[j]*omega[j];
        shiftg.elem[i]=complex(cos(pi2*aux),sin(pi2*aux))*g.elem[i];
    }

    delete [] index;    
    delete [] index_serie;

    return shiftg;
}

grid shiftc (const grid& g, const double &rho)
/*
    Given a grid of Fourier coefficients of a ndim-periodic function
    we compute the composition f(theta+i*rho), that is, we obtain
    the Fourier coefficient of the function shifted to the boundary
*/
{
    int *index;
    int *index_serie;
    grid shiftg(g.ndim,g.nn);

    double pi2,aux;
    pi2=6.2831853071795864769252867665590057684;

    index=new int[g.ndim];
    index_serie=new int[g.ndim];

    for(int i=0;i<g.nelem;i++){
        indices(i,g.nn,index,g.ndim);
        trigo_to_series(g.nn,index,index_serie,g.ndim);
        aux=0.0;
        for(int j=0;j<g.ndim;j++) aux = aux + index_serie[j];
        aux = aux*rho;
        shiftg.elem[i]=exp(-pi2*aux)*g.elem[i];
    }

    delete [] index;    
    delete [] index_serie;

    return shiftg;
}

grid cohomological (const grid& g, double *omega)
/*
   Given a grid of Fourier coefficients of a ndim-periodic function
   we compute the solution of the 1-bite cohomological equation
*/
{
    int *index;
    int *index_serie;
    grid coho(g.ndim,g.nn);

    double pi2, aux, aux1;
    pi2=6.2831853071795864769252867665590057684;

    index=new int[g.ndim];
    index_serie=new int[g.ndim];

    coho.elem[0]=0.0;
    for(int i=1;i<g.nelem;i++){
        indices(i,g.nn,index,g.ndim);
        trigo_to_series(g.nn,index,index_serie,g.ndim);
        aux=0.0;
        for(int j=0;j<g.ndim;j++) aux = aux + index_serie[j]*omega[j];
        
        aux1=abs(g.elem[i]);
        if (aux1<tolgrid) coho.elem[i]=complex(0.0,0.0);
        else coho.elem[i]=g.elem[i]/(1.0-complex(cos(pi2*aux),sin(pi2*aux)));
        //coho.elem[i]=g.elem[i]/(1.0-complex(cos(pi2*aux),sin(pi2*aux)));
    }

    delete [] index;    
    delete [] index_serie;

    return coho;
}

void zerogrid(grid& g)
/*
    No explanation need
*/
{
    for(int i=0;i<g.nelem;i++)
        g.elem[i]=complex(0.0,0.0);
}

void print_gridR(const grid& g)
/*
    Auxiliar function to be used to check intermediate computation.
    It prints the sampled function at the grid
*/
{
    int *index;
    double aux;
    
    index = new int[g.ndim];
    for(int l=0;l<g.nelem;l++){
        indices(l,g.nn,index,g.ndim);
        for (int i=0;i<g.ndim;i++) {
            aux=index[i];
            aux=aux/g.nn[i];
                cout << aux << " " ;
        }
        cout << g.elem[l].real << " " << g.elem[l].imag << endl;
    }
    delete [] index;
}

void print_gridF(const grid& g)
/*
    Auxiliar function to be used to check intermediate computation.
    It prints the Fourier coefficients of the function
*/
{
    int *index;
    int *index_serie;
    
    index = new int[g.ndim];
    index_serie = new int[g.ndim];
    for(int l=0;l<g.nelem;l++){
        indices(l,g.nn,index,g.ndim);
        trigo_to_series(g.nn,index,index_serie,g.ndim);
        
        for (int i=0;i<g.ndim;i++) cout << index_serie[i] << " ";
            cout << g.elem[l].real << " " << g.elem[l].imag << endl;
    }
    delete [] index;
    delete [] index_serie;
}

#endif

