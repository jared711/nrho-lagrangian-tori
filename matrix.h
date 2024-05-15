/******************************************************
This library allows us to manipulate matrices whose entries are
periodic functions. This functions are coded using the library grid.h

See the description at the head of "grid.h" for details on the
class grid, which stores either a discrete (equispaced) sampling
of points on the standard torus or the corresponding DFT.

Some elementary operators (=,+,-,*) have been overloaded
in order to simplify the coding of the parameterization
method.

REMARK: Most of the routines in this library are just an extension
to the corresponding routine in "grid.h" by applying it to the
elements of the matrix. In these cases, we refer to the explanations
in the library "grid.h"

IMPORTANT REMARK: Notice that the same class is used to
store both the sampled points of the funcion and the Fourier coefficients.
Notice also that most of the routines below (deriva, cohomological, etc)
act on the Fourier cofficients. This has de advantage of simplifying
the implementation of the class. Be careful!
******************************************************/
#ifndef MATRIX_H
#define MATRIX_H

#include "complex.h"
#include "grid.h"

class matrix
{
    public:

    int nrows;   // Number of rows
    int ncols;   // Number of columns
    grid **coef; // Sampled function or truncated Fourier series

    matrix();
    matrix(int nr, int nc);
    matrix(int nr, int nc, int ndim, int *nn);
    matrix(const matrix&m);
    ~matrix();

    matrix& operator = (const matrix& m);

    void alloc_matrix(){
        coef=new grid*[nrows];
        for(int i=0;i<nrows;i++){
            coef[i] = new grid[ncols];
        }
    }
};

matrix operator - (const matrix& m);
matrix operator + (const matrix& m1, const matrix& m2);
matrix operator - (const matrix& m1, const matrix& m2);
matrix operator * (const matrix& m1, const matrix& m2);
matrix operator * (const matrix& m, const double& x);
void start_matrix(const double& tol);
matrix fft_F(const matrix& m);
matrix fft_B(const matrix& m);
matrix shift(const matrix& m, double *omega);
matrix cohomological (const matrix& m, double *omega);
matrix diff(const matrix& m);
matrix trans(const matrix& m);
matrix projX(const matrix& m);
matrix projY(const matrix& m);
matrix inv(const matrix& m);
void tail(const matrix& m, double *tg);
void clean(matrix& m);
double norm(const matrix& m);
double normsup(const matrix& m, const double &rho);
double normexp(const matrix& m, const double &rho);
int qrbksb(double **a ,int m, int n, double *b, double *x, const double& tol);
void qrdcmp( double **a , int m, int n, const double& tol);
void aver(const matrix& m, double **avg);
void print_matrixR(const matrix& m);
void print_matrixF(const matrix& m);

double tolqr;

matrix::matrix()
{
    nrows=0;
    ncols=0;
    coef=NULL;
}

matrix::matrix(int nr, int nc)
{
    nrows=nr;
    ncols=nc;

    alloc_matrix();
}

matrix::matrix(int nr, int nc, int ndim, int *nn)
{
    nrows=nr;
    ncols=nc;

    alloc_matrix();
    for(int i=0;i<nrows;i++){
        for(int j=0;j<ncols;j++){
            coef[i][j] = grid(ndim,nn);
        }
    }
}

matrix::matrix(const matrix&m)
{
    nrows=m.nrows;
    ncols=m.ncols;

    alloc_matrix();
    for(int i=0;i<nrows;i++){
        for(int j=0;j<ncols;j++){
            coef[i][j]=m.coef[i][j];
        }
    }
}

matrix::~matrix()
{
    if (coef!=NULL){
        for(int i=0;i<nrows;i++){
            delete [] coef[i];
        }
        delete [] coef;
    }
}

inline matrix& matrix::operator = (const matrix& m)
{
    if (this == &m){
        return *this;
    }
    else{
        if (nrows!=m.nrows || ncols!=m.ncols){
            this->~matrix();
            nrows=m.nrows;
            ncols=m.ncols;
            alloc_matrix();
        }
        for(int i=0;i<nrows;i++){
            for(int j=0;j<ncols;j++){
                coef[i][j]=m.coef[i][j];
            }
        }
        return *this;
    }
}

matrix operator + (const matrix& m1, const matrix& m2)
{
    if ((m1.nrows!=m2.nrows)||(m1.ncols!=m2.ncols))
        cout << "Cuidado con las dimensiones en la + de matrices" << endl;

    matrix suma(m1.nrows,m1.ncols);

    for(int i=0;i<m1.nrows;i++){
        for(int j=0;j<m1.ncols;j++){
            suma.coef[i][j] = m1.coef[i][j] + m2.coef[i][j];
        }
    }
    return suma;
}

matrix operator - (const matrix& m)
{
    matrix menos(m.nrows,m.ncols);

    for(int i=0;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            menos.coef[i][j] = - m.coef[i][j];
        }
    }
    return menos;
}

matrix operator - (const matrix& m1, const matrix& m2)
{
    if ((m1.nrows!=m2.nrows)||(m1.ncols!=m2.ncols))
    cout << "Cuidado con las dimensiones en la - de matrices" << endl;

    matrix resta(m1.nrows,m1.ncols);

    for(int i=0;i<m1.nrows;i++){
        for(int j=0;j<m1.ncols;j++){
            resta.coef[i][j] = m1.coef[i][j] - m2.coef[i][j];
        }
    }
    return resta;
}

matrix operator * (const matrix& m1, const matrix& m2)
{
    if ((m1.ncols!=m2.nrows)) 
    cout << "Cuidado con las dimensiones en la * de matrices" << endl;

    matrix prod(m1.nrows,m2.ncols,m1.coef[0][0].ndim,m1.coef[0][0].nn);
    
    for(int i=0;i<m1.nrows;i++){
        for(int k=0;k<m2.ncols;k++){
            for(int j=0;j<m1.ncols;j++){
                prod.coef[i][k] = prod.coef[i][k] + m1.coef[i][j]*m2.coef[j][k];
            }
        }
    }

    return prod;
}

matrix operator * (const matrix& m, const double& x)
{
    matrix mx(m.nrows,m.ncols);

    for(int i=0;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            mx.coef[i][j] = m.coef[i][j]*x;
        }
    }
    return mx;
}

void start_matrix(const double& tol)
/*
    This defines the global variable tolqr that we use in
    the QR decomposition used to evaluate the inverse of a matrix
*/
{
    tolqr=tol;
}

matrix fft_F(const matrix& m)
{
    matrix fftm(m.nrows,m.ncols);

    for(int i=0;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            fftm.coef[i][j]=fft_F(m.coef[i][j]);
        }
    }
    return fftm;
}

matrix fft_B(const matrix& m)
{
    matrix fftm(m.nrows,m.ncols);

    for(int i=0;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            fftm.coef[i][j]=fft_B(m.coef[i][j]);
        }
    }
    return fftm;
}

matrix shift(const matrix& m, double *omega)
{
    matrix sg(m.nrows,m.ncols);

    for(int i=0;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            sg.coef[i][j]=shift(m.coef[i][j],omega);
        }
    }
    return sg;
}

matrix cohomological (const matrix& m, double *omega)
{
    matrix coho(m.nrows,m.ncols);

    for(int i=0;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            coho.coef[i][j]=cohomological(m.coef[i][j],omega);
        }
    }
    return coho;
}

void tail(const matrix& m, double *tg)
{
    double *tgi;
    int n=m.coef[0][0].ndim;

    tgi = new double[n];

    for (int k=0;k<n;k++) tg[k]=0.0;

    for(int i=0;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            tail(m.coef[i][j],tgi);
            for(int k=0;k<n;k++){
                if (tgi[k]>tg[k]) tg[k]=tgi[k];
            }
        }
    }

    delete [] tgi;
}

void clean(matrix& m)
{
    for(int i=0;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            clean(m.coef[i][j]);
        }
    }
}

double norm(const matrix& m)
{
    double normm, normm0;

    normm=0.0;
    normm0=0.0;
    for(int i=0;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            normm0=norm(m.coef[i][j]);
            if (normm0>normm) normm=normm0;
        }
    }
    return normm;
}

double normsup(const matrix& m, const double &rho)
{
    double normm, normm0;
    
    normm=0.0;
    for(int i=0;i<m.nrows;i++){
        normm0=0.0;
        for(int j=0;j<m.ncols;j++){
            normm0=normm0+normsup(m.coef[i][j],rho);
        }
        if (normm0>normm) normm=normm0;
    }
    return normm;
}

double normexp(const matrix& m, const double &rho)
{
    double normm, normm0;

    normm=0.0;
    for(int i=0;i<m.nrows;i++){
        normm0=0.0;
        for(int j=0;j<m.ncols;j++){
            normm0=normm0+normexp(m.coef[i][j],rho);
        }
        if (normm0>normm) normm=normm0;
    }
    return normm;
}

matrix diff(const matrix& m)
/*
    Given a vector of periodic functions, this routine
    computes the matrix of partial derivatives with respect
    to all the angles.
*/
{
    if (m.ncols!=1) cout << "diff - No estas calculando la derivada de un vector" << endl;

    int nvar=m.coef[0][0].ndim;

    matrix Dm(m.nrows,nvar);

    for (int i=0;i<m.nrows;i++){
        for(int k=0;k<nvar;k++){
            Dm.coef[i][k]=deriva(m.coef[i][0],k);
        }
    }
    return Dm;
}

matrix projX(const matrix& m)
/*
    Given a matrix of periodic funtions (both sampled points of Fourier modes)
    this routine set to 0.0 the elements corresponding to the y-variables
*/
{
    matrix Xm(m);

    for(int i=m.nrows/2;i<m.nrows;i++){
        for(int j=0;j<m.ncols;j++){
            zerogrid(Xm.coef[i][j]);
        }
    }

    return Xm;
}

matrix projY(const matrix& m)
/*
    Given a matrix of periodic funtions (both sampled points of Fourier modes)
    this routine set to 0.0 the elements corresponding to the x-variables
*/
{
    matrix Ym(m);

    for(int i=0;i<m.nrows/2;i++){
        for(int j=0;j<m.ncols;j++){
            zerogrid(Ym.coef[i][j]);
        }
    }

    return Ym;
}

matrix trans(const matrix& m)
/*
    Given a matrix of periodic funtions (both sampled points of Fourier modes)
    this computes the transposed matrix
*/
{
    matrix mt(m.ncols,m.nrows);

    for (int i=0;i<m.nrows;i++){
        for (int j=0;j<m.ncols;j++){
            mt.coef[j][i]=m.coef[i][j];
        }
    }
    return mt;
}

void aver(const matrix& m, double **avg)
{
    for (int i=0;i<m.nrows;i++){
        for (int j=0;j<m.ncols;j++){
            aver(m.coef[i][j],avg[i][j]);
        }
    }
}

void print_matrixR(const matrix& m)
{
    int *index;

    index = new int[m.coef[0][0].ndim];
    for(int l=0;l<m.coef[0][0].nelem;l++){
        indices(l,m.coef[0][0].nn,index,m.coef[0][0].ndim);
        for (int i=0;i<m.coef[0][0].ndim;i++) cout << index[i] << " ";
        cout << endl;
        for (int i=0;i<m.nrows;i++){
            for (int j=0;j<m.ncols;j++){
                cout << m.coef[i][j].elem[l].real << " + i * " << m.coef[i][j].elem[l].imag << ",    ";
            }
            cout << endl;
        }
        cout << "------------------" << endl;
    }
    delete [] index;
}

void print_matrixF(const matrix& m)
{
    int *index;
    int *index_serie;

    index = new int[m.coef[0][0].ndim];
    index_serie = new int[m.coef[0][0].ndim];
    for(int l=0;l<m.coef[0][0].nelem;l++){
        indices(l,m.coef[0][0].nn,index,m.coef[0][0].ndim);
        trigo_to_series(m.coef[0][0].nn,index,index_serie,m.coef[0][0].ndim);

        for (int i=0;i<m.coef[0][0].ndim;i++) cout << index_serie[i] << " ";
        cout << endl;
        for (int i=0;i<m.nrows;i++){
            for (int j=0;j<m.ncols;j++){
                    cout << m.coef[i][j].elem[l].real << " + i * " << m.coef[i][j].elem[l].imag << ",    ";
            }
            cout << endl;
        }
        cout << "------------------" << endl;
    }
    delete [] index;
    delete [] index_serie;
}

matrix inv(const matrix& m)
/*
    Given a matrix of periodic funtions (for sampled points only!)
    this computes the inverse.
*/

{
    if(m.nrows!=m.ncols) cout << "Error en inv: no es una matriz cuadrada" << endl;

    matrix minv(m);
    double **A, *solc, *iden;
    int n=m.ncols;

    if(n==1){
        for(int l=0;l<m.coef[0][0].nelem;l++){
            minv.coef[0][0].elem[l]=1.0/m.coef[0][0].elem[l];
        }
    }

    A=new double*[n];
    solc=new double[n];
    iden=new double[n];
    for (int i=0;i<n;i++){
        A[i]=new double[n];
    }

    for(int l=0;l<m.coef[0][0].nelem;l++){
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                A[i][j]=m.coef[i][j].elem[l].real;
            }
        }
        qrdcmp(A,n,n,tolqr);
        for (int j=0; j<n; j++){
            for (int i=0;i<n;i++) iden[i]=0.0;
            iden[j]=1.0;
            qrbksb(A,n,n,iden,solc,tolqr);
            for (int i=0;i<n;i++) minv.coef[i][j].elem[l]=solc[i];
        }
    }

    for (int i=0;i<n;i++){
        delete [] A[i];
    }
    delete [] A;
    delete [] solc;
    delete [] iden;

    return minv;
}

void qrdcmp(double **a , int m, int n, const double& tol)
{
    double *w,alf,bet,gam,aux,alfa;
    int i,j,k;
    double val2=2.0;
    double val1=1.0;
    
    w= new double[n];
    if (w==NULL)
    {
        cout << "qr. No hi ha espai pel vector de treball" << endl;
        return;
    }
    
    for (k=0;k<n;k++)
    {
        alf=0.0;
        for (i=k+1;i<m;i++) alf=alf+a[i][k]*a[i][k];
        gam=sqrt(a[k][k]*a[k][k]+alf);
        bet=( a[k][k]> 0.0 ? a[k][k]+gam : a[k][k]-gam);
        if (abs(bet)<tol) bet=1.0;
        alfa=-val2/(val1+alf/(bet*bet));
        w[k]=a[k][k];
        for (i=k+1;i<m;i++)
        {
            aux=a[i][k]/bet;
            w[k]=w[k]+aux*a[i][k];
            a[i][k]=aux;
        }
        w[k]=w[k]*alfa;
        for (j=k+1;j<n;j++)
        {
            w[j]=a[k][j];
            for (i=k+1;i<m;i++) w[j]=w[j]+a[i][j]*a[i][k];
            w[j]=w[j]*alfa;
        }
        for (j=k;j<n;j++) a[k][j]=a[k][j]+w[j];
        for (i=k+1;i<m;i++) for(j=k+1;j<n;j++) a[i][j]=a[i][j]+w[j]*a[i][k];
    }
    delete [] w;
}

int qrbksb(double **a ,int m, int n, double *b, double *x, const double& tol)
{
    int j,k;
    double bet,w,alfa;
    double val2=2.0;
    
    for (j=0;j<m;j++) x[j]=b[j];
    for (k=0;k<n;k++)
    {
        bet=1.0;
        for (j=k+1;j<m;j++) bet=bet+a[j][k]*a[j][k];
        alfa=-val2/bet;
        w=x[k];
        for (j=k+1;j<m;j++) w=w+x[j]*a[j][k];
        w=w*alfa;
        x[k]=x[k]+w;
        for (j=k+1;j<m;j++) x[j]=x[j]+w*a[j][k];
    }
    for (k=n-1;k>=0;k--)
    {
        w=0.0;
        
        for (j=k+1;j<n;j++) w=w+a[k][j]*x[j];
        x[k]=(x[k]-w)/a[k][k];
    }
    return(1);
}

#endif

