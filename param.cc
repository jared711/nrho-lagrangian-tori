/********************************************************************
------------------------------------------------
Alex Haro          http://www.maia.ub.es/~alex/
Alejandro Luque    http://www.maia.ub.es/~luque/
------------------------------------------------

Last version October 21, 2014
------------------------------------------------
Edited by Jared Blanchard (jaredb711@gmail.com) May 2024

Description of the routine:
------------------------------------------------
    This routine continues numerically an invariant torus of prefixed
    frequencies using the parameterization method. We follow the
    notation and implementation described in

    A. Haro, A. Luque: The parameterization method in KAM theory.
    Chapter 3 in "A. Haro. J.M. Mondelo, J.Ll. Figueras, A. Luque. and M.Canadell:
    The parameterization method for invariant manifolds: from rigorous
    results to effective computations". Applied Mathematica Sciences, Springer.

    We refer the user to the extense literature on the parameterization
    method and related algorithms. In particular to the work of
    R. de la Llave, A. González, A. Jorba, J. Villanueva, E. Fontich, Y. Sire, G. Huguet,
    R. C. Calleja, A. Celletti, etc

    The user has to provide the following subroutines for a given particular
    problem

    1) Evaluation of the studied map
       void map(complex *z, complex *fz, complex **Dfz, complex *depfz)

    2) Evaluation of the metric
       void gform_standard(complex *z, complex **Omegaz)

    3) Evaluation of the symplectic form

    4) Evaluation of a transversal field (only if we don't want to use a metric))
       void normal0_standard(matrix &N0, int *nn, int nelem)

    In this code, subroutines are provided for the case of the standard map
    and the Froeschle map.
********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <cstring>
#include "complex.h"
#include "grid.h"
#include "matrix.h"

extern "C"
{
#include "campvp.h"
#include "rtbphp.h"
#include "fluxvp.h"
#include "seccp.h"
#include "scread.h"
#include "utils.h"
}
// #include "utils/utils.h"
// #include "vfld/vfldsh.h"
// #include "vfld/rtbphsh.h"
// #include "stfWPeutils.h"
// #include "stfWPe_manifxflw.h"
// #include "maniutils.h"

using namespace std;

#define myreal double

// Question: Do I need to redefine these constants? YES
// #define DTOR (1)    // Dimension of the invariant torus
// #define DMAP (2)    // Dimension of phase space
// #define NPAR (0)    // Number of (constant) parameters in the map
// #define MNEW (10)   // Maximum number of Newton iterations
// #define MAXF (8192) // Maximum number of Fourier coefficients allowed

// CR3BP map
#define DTOR (2)    // Dimension of the invariant torus because I'm doing the generator of the full lagrangian torus
#define DMAP (4)    // Dimension of phase space (5/8/24, we've reduced the dimnsion to 4 because we're constraining C and y=0 with a Poincare map)
#define NPAR (2)    // Number of (constant) parameters in the map
#define MNEW (10)   // Maximum number of Newton iterations
#define MAXF (8192) // Maximum number of Fourier coefficients allowed

int position(int *nn, int *index, int ndim);
void indices(int pos, int *nn, int *index, int ndim);
void realloc_torus(matrix &paramR, matrix &paramF,
                   matrix &paramR0, matrix &paramF0,
                   int &nelem, int *nn, int *newnn);
int kam_torus(matrix &paramR, matrix &paramF, myreal *omega, myreal &error, int *nn, int nelem,
              int &tail0, int *tails, int Case,
              void (*map)(complex *, complex *, complex **, complex *),
              void (*sform)(complex *, complex **),
              void (*gform)(complex *, complex **),
              void (*normal0)(matrix &, int *, int));
void map_froeschle(complex *z, complex *fz, complex **Dfz, complex *depfz);
void sform_froeschle(complex *z, complex **Omegaz);
void gform_froeschle(complex *z, complex **Metricz);
void normal0_froeschle(matrix &N0, int *nn, int nelem);
void map_standard(complex *z, complex *fz, complex **Dfz, complex *depfz);
void sform_standard(complex *z, complex **Omegaz);
void gform_standard(complex *z, complex **Metricz);
void normal0_standard(matrix &N0, int *nn, int nelem);

/* functions created by Jared Blanchard May, 2024*/
void nu(complex *z, double *x, double (*Dnu)[4]);
void get_dp3(complex *z, double *dp3);
void map_CR3BP(complex *z, complex *fz, complex **Dfz, complex *depfz);
void sform_CR3BP(complex *z, complex **Omegaz);
void gform_CR3BP(complex *z, complex **Metricz);
void normal0_CR3BP(matrix &N0, int *nn, int nelem);
void wrtf(int n, int nv, double t, double x[], int aon, void *prm);
double get_p3(double q1, double q2, double q3, double p1, double p2, double mu);
void state2ham(double state[]);

/* General global variables */
myreal pi, pi2, val0, val05, val1, val2, val3, val4, val5, val10;
myreal global_twist; // Norm of inverse of torsion matrix T(theta)
myreal toltail;      // Tolerance on the size of the tails of the parameterization
myreal tolinva;      // Tolerance on the error of invariance
myreal tolinte;      // Tolerance on intermediate computation (e.g. matrix inverses)

/* Global variables of the problem */
myreal *lambda; // Fixed parameters or constants
myreal epsilon; // Continuation parameter

fstream file_torus, file_input;

#define N RTBPHP_N

int main(int argc, char *argv[])
{
    double auxdouble;
    char name[80];

    int *nn, *index, nelem, *newnn;
    complex *z, *fz, **Dfz, *depfz;
    myreal *omega, tol;
    myreal error, aux, deps, epsilon0;
    matrix paramR, paramF, paramR0, paramF0;
    int *tails, tail0;
    int power2;
    int conv, iter;
    nn = new int[DTOR];
    newnn = new int[DTOR];
    index = new int[DTOR];
    omega = new myreal[DTOR];
    tails = new int[DTOR];

    /**** START We set global variables ****/
    lambda = new myreal[NPAR];
    val0 = 0.0;
    val1 = 1.0;
    val2 = 2.0;
    val05 = val1 / val2;
    val3 = 3.0;
    val4 = 4.0;
    val5 = 5.0;
    val10 = 10.0;
    pi2 = 6.2831853071795864769252867665590057684;
    /**** END   We set global variables ****/

    cout << scientific;
    cout.precision(15);

    /**** START We read an invariant torus from the file ****/
    file_input.open(argv[1], ios::in);
    file_input << scientific;
    file_input >> toltail; // tolerance of the tail from eq. 4.104
    cout << "#toltail: " << toltail << endl;
    file_input >> tolinva; //
    cout << "#tolinva: " << tolinva << endl;
    file_input >> tolinte; // toerance of the error at the grid
    cout << "#tolinte: " << tolinte << endl;
    start_matrix(tolinte);
    start_grid(tolinte);
    for (int i = 0; i < DTOR; i++)
    {
        file_input >> omega[i];
        cout << "#omega[" << i << "]: " << omega[i] << endl;
    }
    int auxint;
    file_input >> epsilon;
    cout << "#epsilon: " << epsilon << endl;
    for (int j = 0; j < NPAR; j++)
    {
        file_input >> lambda[j];
        cout << "#lambda[" << j << "]: " << lambda[j] << endl;
    }
    cout << "lambda[0] should be H0" << endl;
    cout << "lambda[1] should be mu" << endl;
    for (int j = 0; j < DTOR; j++)
    {
        file_input >> nn[j]; // dimesions of the torus that you're trying to compute
        cout << "#nn[" << j << "]: " << nn[j] << endl;
    }
    file_input >> deps;
    cout << "#deps: " << deps << endl;
    file_input >> auxint;
    cout << "#auxint: " << auxint << endl;
    paramR = matrix(DMAP, 1, DTOR, nn); // Memory for the torus
    nelem = 1;                          // initialized as one, but becomes the number of elements in the grid
    for (int i = 0; i < DTOR; i++)
        nelem *= nn[i]; // Number of elements in the grid is the product of the components of nn
    if (auxint != 0)
    { // if auxint is not 0, then we read the torus from the file
        for (int l = 0; l < nelem; l++)
        {
            // this line seems worthless, maybe I'm supposed to use it to keep the rows of the torus separate
            // since it does this for the first DTOR elements, I need to put DTOR (2) elements in before each state
            // It might be better to just rewrite this to import a csv file
            for (int j = 0; j < DTOR; j++)
            {
                file_input >> auxint; // assign the next line into auxint (which is an integer)
                cout << "auxint: " << auxint << endl;
            }
            for (int j = 0; j < DMAP; j++)
            {
                paramR.coef[j][0].elem[l] = val0;
                file_input >> paramR.coef[j][0].elem[l].real;
                cout << "paramR.coef[" << j << "][0].elem[" << l << "]: " << paramR.coef[j][0].elem[l].real << endl;
            }
        }
    }
    else
    {
        for (int l = 0; l < nelem; l++)
        {
            indices(l, nn, index, DTOR);
            for (int j = 0; j < DTOR; j++)
                paramR.coef[j][0].elem[l] = val0;
            for (int j = 0; j < DTOR; j++)
                paramR.coef[j + DTOR][0].elem[l] = omega[j];
        }
    }
    file_input.close();
    /**** END   We read an invariant torus from the file ****/


    paramF = fft_F(paramR);
    // print out all components of paramR and paramF
    // for (int i = 0; i < DMAP; i++)
    // {
    //     for (int j = 0; j < nelem; j++)
    //     {
    //         cout << "paramR.coef[" << i << "][0].elem[" << j << "]: " << paramR.coef[i][0].elem[j].real << endl;
    //     }
    // }
    // for (int i = 0; i < DMAP; i++)
    // {
    //     for (int j = 0; j < nelem; j++)
    //     {
    //         cout << "paramF.coef[" << i << "][0].elem[" << j << "]: " << paramF.coef[i][0].elem[j].real << endl;
    //     }
    // }
    // return 0;

    /* If the initial torus is not invariant, we shold uncomment the following line: */
    conv = kam_torus(paramR,paramF,omega,error,nn,nelem,tail0,tails,2,map_CR3BP,sform_CR3BP,gform_CR3BP,normal0_CR3BP);

    conv = kam_torus(paramR,paramF,omega,error,nn,nelem,tail0,tails,2,map_CR3BP,sform_CR3BP,gform_CR3BP,normal0_CR3BP);

    // return 0;


    paramR0 = matrix(paramR);
    paramF0 = matrix(paramF);
    epsilon0 = epsilon;

    /* not using this anymore 5/8/24 */
    // double mu = 1.901109735892602E-7;
    // double t_value = 0;
    // double *t = &t_value; //
    // double h_value = 1E-1;
    // double *h = &h_value; // this is type casting to make a pointer to the number 1E-1
    // vfldvp_t rtbphp;
    // double x[6] = {1.0025904252283664E+0, 1.2103594028238294E-20, 4.8802791235671040E-3, -1.3013300877803077E-15, -5.4593457318909567E-3, -1.2951908078310448E-14};

    // // This was to check that seccp actually works correctly (commented out 7/1/24)
    // double x[6] = {-0.9975334497794613, 0, -0.00489434320310318, 1.383426648384139E-16, -1.002834494433839, -1.991818118201947E-16};
    // double mu = 1.901109735892602e-7;
    // lambda[0] = rtbphp_h(6, 0, &mu, x);
    // int ibck = 0; // If ibck==1, backward in time (forward if ==0)
    // int isiggrad = 1;
    // double tolJM = 1e-12;
    // double maxts = 13;
    // int ivb = 1;
    // int nsecss = 1;
    // int nsec = 1; // nsec1 is the number of passes through each section. Since we have only one section, I'm going to make it just an int
    // double cp[7] = {0, 0, 1, 0, 0, 0, 0};
    // double t = 0;
    // double h = fluxvp_pas0;
    // if (ibck)
    //     h = -h;
    // // FILE *fp = NULL;
    // FILE *fp = fopen("mypoint.txt", "w");
    // int indict = 0;
    // indict = seccp(6, 6 /*nv*/, 0 /*np*/, rtbphp /*camp*/, &mu /*prm*/, &t, x, &h, cp,
    //                nsec, isiggrad, tolJM, ivb, 0 /*idt*/, NULL /*dt*/, wrtf, fp, maxts);
    // cout << "indict: " << indict << endl;
    // return 0;
    // state2ham(x);

    // return 0;

    // cout << "hey"   << endl;
    // flowvp(1/*N*/,6/*n number of states?*/,0/*number of parameters?*/,&mu/*parameter*/,rtbphp/*vfld*/,
    //   t,x,h,1.0/*T*/,NULL/*forb*/,0/*ivb*/);
    // cout << "hey2"   << endl;

    // for (int i = 0; i < 6; i++)
    // {
    //    cout << "x[" << i << "]: " << x[i] << endl;
    // }
    // return 0;

    /*This is the debugging block, to make sure I'm reading everything in correctly*/

    // return 0;
    /* END DEBUG BLOCK*/

    /**** START Continuation with respect to epsilon ****/
    int fail = 0;
    do
    {
        paramR = paramR0;
        paramF = paramF0;
        epsilon = epsilon0 + deps;

        /**** START Newton method to correct the invariant torus ****/
        cout << "# We try to compute the torus for epsilon=" << epsilon << endl;
        iter = 0;
        do
        {
            cout << "# Iteration " << iter + 1 << " : " << endl;
            conv = kam_torus(paramR, paramF, omega, error, nn, nelem, tail0, tails, 1, map_CR3BP, sform_CR3BP, gform_CR3BP, normal0_CR3BP);
            // conv = kam_torus(paramR, paramF, omega, error, nn, nelem, tail0, tails, 1, map_standard, sform_standard, gform_standard, normal0_standard);
            // conv = kam_torus(paramR,paramF,omega,error,nn,nelem,tail0,tails,3,map_froeschle,sform_froeschle,gform_froeschle);
            if (tail0 == 1)
            {
                /* If the tail is too large, we refine the grid and restart the computation */
                if (nelem == MAXF)
                {
                    cout << "# We reched the maximum number of Fourier modes" << endl;
                    delete[] nn;
                    delete[] newnn;
                    delete[] index;
                    delete[] omega;
                    delete[] tails;
                    delete[] lambda;
                    return 0;
                }
                for (int i = 0; i < DTOR; i++)
                {
                    newnn[i] = nn[i];
                    if (tails[i] == 1)
                        newnn[i] = 2 * nn[i];
                }
                realloc_torus(paramR, paramF, paramR0, paramF0, nelem, nn, newnn);
            }
            iter++;
        } while (conv == 0 && iter < MNEW && tail0 == 0);
        /**** END   Newton method to correct the invariant torus ****/

        if (conv == 1)
        {
            cout.precision(15);
            /* If we converge, we update the last computed torus and store the results */
            paramR0 = paramR;
            paramF0 = paramF;
            epsilon0 = epsilon;

            /**** START We print the information related to the successful continuation step ****/
            cout << epsilon << " ";
            for (int j = 0; j < DTOR; j++)
                cout << nn[j] << " ";
            cout << global_twist << " ";
            for (int j = 0; j < DTOR; j++)
                cout << normsobo(paramF.coef[j][0], val2) << " ";
            cout << error << endl;
            /**** END   We print the information related to the successful continuation step ****/

            /**** START We save the computed invariant torus in a separated file ****/
            sprintf(name, "output_torus%.4lf", epsilon);
            file_torus.open(name, ios::out);
            file_torus << scientific;
            file_torus.precision(15);
            file_torus << toltail << endl;
            file_torus << tolinva << endl;
            file_torus << tolinte << endl;
            for (int j = 0; j < DTOR; j++)
                file_torus << omega[j] << endl;
            file_torus << epsilon << endl;
            for (int j = 0; j < NPAR; j++)
                file_torus << lambda[j] << endl;
            for (int j = 0; j < DTOR; j++)
                file_torus << nn[j] << endl;
            file_torus << deps << endl;
            file_torus << 1 << endl;
            for (int l = 0; l < nelem; l++)
            {
                indices(l, nn, index, DTOR);
                for (int j = 0; j < DTOR; j++)
                    file_torus << index[j] << " ";
                for (int j = 0; j < DMAP - 1; j++)
                    file_torus << paramR.coef[j][0].elem[l].real << " ";
                file_torus << paramR.coef[DMAP - 1][0].elem[l].real << endl;
            }
            file_torus.close();
            /**** END   We save the computed invariant torus in a separated file ****/
        }
        else
        {
            if (tail0 == 1)
                cout << "# We need more Fourier modes" << endl;
            else
            {
                cout << "# Newton method does not converge" << endl;
                deps = deps / val10;
                if (deps < 1e-5)
                    fail = 1;
            }
            cout << "###########################################" << endl;
        }
    } while (fail == 0);
    /**** END Continuation with respect to epsilon ****/

    delete[] nn;
    delete[] newnn;
    delete[] index;
    delete[] omega;
    delete[] tails;
    delete[] lambda;

    return 0;
}

void realloc_torus(matrix &paramR, matrix &paramF,
                   matrix &paramR0, matrix &paramF0,
                   int &nelem, int *nn, int *newnn)
{
    matrix newparamF;
    int *index, *indexserie;
    int pos;

    index = new int[DTOR];
    indexserie = new int[DTOR];

    newparamF = matrix(DMAP, 1, DTOR, newnn);

    for (int l = 0; l < nelem; l++)
    {
        indices(l, nn, index, DTOR);
        trigo_to_series(nn, index, indexserie, DTOR);
        series_to_trigo(newnn, indexserie, index, DTOR);
        pos = position(newnn, index, DTOR);
        for (int j = 0; j < DMAP; j++)
            newparamF.coef[j][0].elem[pos] = paramF.coef[j][0].elem[l];
    }
    paramF = newparamF;

    for (int l = 0; l < nelem; l++)
    {
        indices(l, nn, index, DTOR);
        trigo_to_series(nn, index, indexserie, DTOR);
        series_to_trigo(newnn, indexserie, index, DTOR);
        pos = position(newnn, index, DTOR);
        for (int j = 0; j < DMAP; j++)
            newparamF.coef[j][0].elem[pos] = paramF0.coef[j][0].elem[l];
    }
    paramF0 = newparamF;

    paramR = fft_B(paramF);
    paramR0 = fft_B(paramF0);

    nelem = 1;
    for (int i = 0; i < DTOR; i++)
    {
        nelem = nelem * newnn[i];
        nn[i] = newnn[i];
    }

    delete[] index;
    delete[] indexserie;
}

int kam_torus(matrix &paramR, matrix &paramF, myreal *omega, myreal &error, int *nn, int nelem,
              int &tail0, int *tails, int Case, /**/
              void (*map)(complex *, complex *, complex **, complex *),
              void (*sform)(complex *, complex **),
              void (*gform)(complex *, complex **),
              void (*normal0)(matrix &, int *, int))
{
    complex *z, *fz, **Dfz, **Omegaz, **Metricz, *depfz, auxc;
    myreal **twist0, **neweta0;
    myreal **invT, *solc, *iden, tolqr, **xiN0;
    myreal aux;
    int *index;
    myreal *tg;
    matrix FparamR(DMAP, 1, DTOR, nn);
    matrix ErrorR(DMAP, 1, DTOR, nn);
    matrix ErrorF(DMAP, 1, DTOR, nn);
    matrix DFKR(DMAP, DMAP, DTOR, nn);
    matrix OmegaKR(DMAP, DMAP, DTOR, nn);
    matrix OmegaKF(DMAP, DMAP, DTOR, nn);
    matrix OmegaKshiftF(DMAP, DMAP, DTOR, nn);
    matrix OmegaKshiftR(DMAP, DMAP, DTOR, nn);
    matrix MetricKR(DMAP, DMAP, DTOR, nn);
    matrix KshiftF(DMAP, 1, DTOR, nn);
    matrix KshiftR(DMAP, 1, DTOR, nn);
    matrix DparamF(DMAP, DTOR, DTOR, nn);
    matrix DparamR(DMAP, DTOR, DTOR, nn);
    matrix LR(DMAP, DTOR, DTOR, nn);
    matrix LF(DMAP, DTOR, DTOR, nn);
    matrix LshiftR(DMAP, DTOR, DTOR, nn);
    matrix LshiftF(DMAP, DTOR, DTOR, nn);
    matrix GR(DTOR, DTOR, DTOR, nn);
    matrix AR(DTOR, DTOR, DTOR, nn);
    matrix BR(DTOR, DTOR, DTOR, nn);
    matrix NR(DMAP, DTOR, DTOR, nn);
    matrix NF(DMAP, DTOR, DTOR, nn);
    matrix NshiftR(DMAP, DTOR, DTOR, nn);
    matrix NshiftF(DMAP, DTOR, DTOR, nn);
    matrix etaLR(DTOR, 1, DTOR, nn);
    matrix etaNR(DTOR, 1, DTOR, nn);
    matrix etaLF(DTOR, 1, DTOR, nn);
    matrix etaNF(DTOR, 1, DTOR, nn);
    matrix RetaNR(DTOR, 1, DTOR, nn);
    matrix RetaNF(DTOR, 1, DTOR, nn);
    matrix newetaR(DTOR, 1, DTOR, nn);
    matrix newetaF(DTOR, 1, DTOR, nn);
    matrix xiNR(DTOR, 1, DTOR, nn);
    matrix xiLR(DTOR, 1, DTOR, nn);
    matrix xiLF(DTOR, 1, DTOR, nn);
    matrix twistR(DTOR, DTOR, DTOR, nn);
    matrix newparamR(DMAP, 1, DTOR, nn);
    matrix newparamF(DMAP, 1, DTOR, nn);

    index = new int[DTOR];
    z = new complex[DMAP];
    fz = new complex[DMAP];
    depfz = new complex[DMAP];
    Dfz = new complex *[DMAP];
    Omegaz = new complex *[DMAP];
    Metricz = new complex *[DMAP];
    for (int i = 0; i < DMAP; i++)
    {
        Dfz[i] = new complex[DMAP];
        Omegaz[i] = new complex[DMAP];
        Metricz[i] = new complex[DMAP];
    }
    tg = new myreal[DTOR];
    twist0 = new myreal *[DTOR];
    neweta0 = new myreal *[DTOR];
    xiN0 = new myreal *[DTOR];
    for (int i = 0; i < DTOR; i++)
    {
        twist0[i] = new myreal[DTOR];
        neweta0[i] = new myreal[1];
        xiN0[i] = new myreal[1];
    }
    invT = new myreal *[DTOR];
    solc = new myreal[DTOR];
    iden = new myreal[DTOR];
    for (int i = 0; i < DTOR; i++)
    {
        invT[i] = new myreal[DTOR];
    }

    /*****************************************************************
    STEP 0 Evaluation of the tail
    *****************************************************************/
    tail(paramF, tg); // The tail is defined in equation 4.104 of the book Haro A., et al. (2016) The parameterization method for invariant manifolds: From rigorous results to effective computations

    tail0 = 0;
    cout.precision(3);
    cout << "#     - Size of the grid: ";
    for (int i = 0; i < DTOR; i++)
        cout << nn[i] << " ";
    cout << endl;
    cout << "#     - Tails of the parameterization: ";
    for (int i = 0; i < DTOR; i++)
    {                         // goes through all of tg and if a single one is larger than toltail, then we'll stop the computation (step 4 of algorithm 4.32)
        cout << tg[i] << " "; // A big tail means that the parameterization is not good and we need to fix something (e.g. doubling the number of Fourier Modes)
        if (tg[i] > toltail)
        {
            tails[i] = 1;
            tail0 = 1;
        }
        else
            tails[i] = 0;
    }
    cout << endl;
    if (tail0 == 1)
        return 0;
    clean(paramF);
    paramR = fft_B(paramF);

    /*****************************************************************
    STEP 1 Evaluation of the invariance error
    *****************************************************************/
    for (int l = 0; l < nelem; l++)
    {
        indices(l, nn, index, DTOR);
        for (int i = 0; i < DMAP; i++)
        {
            z[i] = paramR.coef[i][0].elem[l];
            // if (i < DTOR) // Jared wants to comment out these lines on 6/27/24
                // z[i] = z[i] + ((double)index[i]) / ((double)nn[i]);
        }
        (*map)(z, fz, Dfz, depfz);
        (*sform)(z, Omegaz);
        if (Case != 1)
            (*gform)(z, Metricz);
        for (int i = 0; i < DMAP; i++)
        {
            FparamR.coef[i][0].elem[l] = fz[i];
            for (int j = 0; j < DMAP; j++)
            {
                DFKR.coef[i][j].elem[l] = Dfz[i][j];
                OmegaKR.coef[i][j].elem[l] = Omegaz[i][j];
                if (Case != 1)
                    MetricKR.coef[i][j].elem[l] = Metricz[i][j];
            }
        }
    }

    KshiftF = shift(paramF, omega);
    KshiftR = fft_B(KshiftF);

    for (int l = 0; l < nelem; l++)
    {
        indices(l, nn, index, DTOR);
        for (int i = 0; i < DTOR; i++)
        {
            KshiftR.coef[i][0].elem[l] = KshiftR.coef[i][0].elem[l] + ((double)index[i]) / ((double)nn[i]) + omega[i];
        }
    }

    ErrorR = FparamR - KshiftR;

    ErrorF = fft_F(ErrorR);
    error = norm(ErrorF);
    cout << "#     - Error of invariance: ";
    cout << error << endl;

    if (error < tolinva)
    {
        cout << "#     - No correction is needed!" << endl;
        delete[] index;
        for (int i = 0; i < DMAP; i++)
        {
            delete[] Dfz[i];
            delete[] Omegaz[i];
            delete[] Metricz[i];
        }
        delete[] z;
        delete[] fz;
        delete[] depfz;
        delete[] Dfz;
        delete[] Omegaz;
        delete[] Metricz;
        delete[] tg;
        for (int i = 0; i < DTOR; i++)
        {
            delete[] twist0[i];
            delete[] neweta0[i];
            delete[] invT[i];
            delete[] xiN0[i];
        }
        delete[] twist0;
        delete[] neweta0;
        delete[] invT;
        delete[] solc;
        delete[] iden;
        delete[] xiN0;

        for (int i = 0; i < DTOR; i++)
            tails[i] = 0;
        return 1;
    }

    /*****************************************************************
    STEP 2 Construction of the symplectic frame
    Here I need to do something with creating L
    *****************************************************************/
    DparamF = diff(paramF);
    DparamR = fft_B(DparamF);

    LR = DparamR;
    // for (int i = 0; i < DTOR; i++)
    // LR.coef[i][i] = LR.coef[i][i] + val1; /*Alex 5/3 get rid of the val1*/

    if (Case == 1)
    { /*Cases in Book*/
        matrix N0R(DMAP, DTOR, DTOR, nn);

        (*normal0)(N0R, nn, nelem); // This is the only place where the normal form gets used
        GR = -trans(LR) * OmegaKR * N0R;
        BR = inv(GR);
        AR = -trans(BR) * trans(N0R) * OmegaKR * N0R * BR * val05;
        NR = LR * AR + N0R * BR;
    }
    else if (Case == 2)
    {
        GR = trans(LR) * MetricKR * LR;
        BR = inv(GR);
        AR = trans(BR) * trans(LR) * MetricKR * inv(OmegaKR) * MetricKR * LR * BR * val05;
        NR = LR * AR - inv(OmegaKR) * MetricKR * LR * BR;
    }
    else if (Case == 3)
    {
        GR = trans(LR) * MetricKR * LR;
        BR = inv(GR);
        NR = inv(MetricKR) * OmegaKR * LR * BR;
    }
    LF = fft_F(LR);
    NF = fft_F(NR);

    /*****************************************************************
    STEP 3 Computation of the correction on the symplectic frame
    *****************************************************************/
    LshiftF = shift(LF, omega);
    LshiftR = fft_B(LshiftF);
    NshiftF = shift(NF, omega);
    NshiftR = fft_B(NshiftF);
    OmegaKF = fft_F(OmegaKR);
    OmegaKshiftF = shift(OmegaKF, omega);
    OmegaKshiftR = fft_B(OmegaKshiftF);
    etaLR = -trans(NshiftR) * OmegaKshiftR * ErrorR;
    etaNR = trans(LshiftR) * OmegaKshiftR * ErrorR;
    twistR = trans(NshiftR) * OmegaKshiftR * DFKR * NR;
    etaNF = fft_F(etaNR);
    RetaNF = cohomological(etaNF, omega);
    RetaNR = fft_B(RetaNF);
    newetaR = etaLR - twistR * RetaNR;

    aver(twistR, twist0);
    tolqr = tolinte;
    qrdcmp(twist0, DTOR, DTOR, tolqr);
    global_twist = val0;
    for (int j = 0; j < DTOR; j++)
    {
        for (int i = 0; i < DTOR; i++)
            iden[i] = val0;
        iden[j] = val1;
        qrbksb(twist0, DTOR, DTOR, iden, solc, tolqr);
        for (int i = 0; i < DTOR; i++)
        {
            invT[i][j] = solc[i];
            global_twist = global_twist + invT[i][j];
        }
    }
    aver(newetaR, neweta0);
    cout << "#     - Norm inverse twist: ";
    cout << global_twist << endl;

    for (int i = 0; i < DTOR; i++)
    {
        for (int k = 0; k < 1; k++)
        {
            xiN0[i][k] = 0.0;
            for (int j = 0; j < DTOR; j++)
            {
                xiN0[i][k] = xiN0[i][k] + invT[i][j] * neweta0[j][k];
            }
        }
    }
    xiNR = RetaNR;

    for (int l = 0; l < nelem; l++)
    {
        for (int i = 0; i < DTOR; i++)
        {
            xiNR.coef[i][0].elem[l] = xiNR.coef[i][0].elem[l] + xiN0[i][0];
        }
    }

    newetaR = etaLR - twistR * xiNR;
    newetaF = fft_F(newetaR);
    xiLF = cohomological(newetaF, omega);
    xiLR = fft_B(xiLF);

    /*****************************************************************
    STEP 4 New parameterization
    *****************************************************************/
    newparamR = paramR + LR * xiLR + NR * xiNR;
    newparamF = fft_F(newparamR);

    delete[] index;
    for (int i = 0; i < DMAP; i++)
    {
        delete[] Dfz[i];
        delete[] Omegaz[i];
        delete[] Metricz[i];
    }
    for (int i = 0; i < DTOR; i++)
    {
        delete[] twist0[i];
        delete[] neweta0[i];
        delete[] invT[i];
        delete[] xiN0[i];
    }
    delete[] z;
    delete[] fz;
    delete[] depfz;
    delete[] Dfz;
    delete[] Omegaz;
    delete[] Metricz;
    delete[] tg;
    delete[] twist0;
    delete[] neweta0;
    delete[] invT;
    delete[] solc;
    delete[] iden;
    delete[] xiN0;

    paramR = newparamR;
    paramF = newparamF;

    /* This is just to show the size of the correction */
    newparamR = LR * xiLR + NR * xiNR;
    newparamF = fft_F(newparamR);
    aux = norm(newparamF);
    cout << "#     - Norm of the correction: ";
    cout << aux << endl;

    return 0;
}

void map_froeschle(complex *z, complex *fz, complex **Dfz, complex *depfz)
{
    fz[2] = z[2] + lambda[0] * sin(pi2 * z[0]) / pi2 + epsilon * sin(pi2 * (z[0] + z[1])) / pi2;
    fz[3] = z[3] + lambda[1] * sin(pi2 * z[1]) / pi2 + epsilon * sin(pi2 * (z[0] + z[1])) / pi2;
    fz[0] = z[0] + fz[2];
    fz[1] = z[1] + fz[3];

    Dfz[2][0] = lambda[0] * cos(pi2 * z[0]) + epsilon * cos(pi2 * (z[0] + z[1]));
    Dfz[2][1] = epsilon * cos(pi2 * (z[0] + z[1]));
    Dfz[2][2] = val1;
    Dfz[2][3] = val0;

    Dfz[3][0] = epsilon * cos(pi2 * (z[0] + z[1]));
    Dfz[3][1] = lambda[1] * cos(pi2 * z[1]) + epsilon * cos(pi2 * (z[0] + z[1]));
    Dfz[3][2] = val0;
    Dfz[3][3] = val1;

    Dfz[0][0] = val1 + Dfz[2][0];
    Dfz[0][1] = Dfz[2][1];
    Dfz[0][2] = Dfz[2][2];
    Dfz[0][3] = Dfz[2][3];

    Dfz[1][0] = Dfz[3][0];
    Dfz[1][1] = val1 + Dfz[3][1];
    Dfz[1][2] = Dfz[3][2];
    Dfz[1][3] = Dfz[3][3];

    depfz[0] = sin(pi2 * (z[0] + z[1])) / pi2;
    depfz[1] = sin(pi2 * (z[0] + z[1])) / pi2;
    depfz[2] = sin(pi2 * (z[0] + z[1])) / pi2;
    depfz[3] = sin(pi2 * (z[0] + z[1])) / pi2;
}

void sform_froeschle(complex *z, complex **Omegaz)
{
    Omegaz[0][0] = val0;
    Omegaz[0][1] = val0;
    Omegaz[0][2] = -val1;
    Omegaz[0][3] = val0;

    Omegaz[1][0] = val0;
    Omegaz[1][1] = val0;
    Omegaz[1][2] = val0;
    Omegaz[1][3] = -val1;

    Omegaz[2][0] = val1;
    Omegaz[2][1] = val0;
    Omegaz[2][2] = val0;
    Omegaz[2][3] = val0;

    Omegaz[3][0] = val0;
    Omegaz[3][1] = val1;
    Omegaz[3][2] = val0;
    Omegaz[3][3] = val0;
}

void gform_froeschle(complex *z, complex **Metricz)
{
    Metricz[0][0] = val1;
    Metricz[0][1] = val0;
    Metricz[0][2] = val0;
    Metricz[0][3] = val0;

    Metricz[1][0] = val0;
    Metricz[1][1] = val1;
    Metricz[1][2] = val0;
    Metricz[1][3] = val0;

    Metricz[2][0] = val0;
    Metricz[2][1] = val0;
    Metricz[2][2] = val1;
    Metricz[2][3] = val0;

    Metricz[3][0] = val0;
    Metricz[3][1] = val0;
    Metricz[3][2] = val0;
    Metricz[3][3] = val1;
}

void normal0_froeschle(matrix &N0, int *nn, int nelem)
{
    for (int l = 0; l < nelem; l++)
    {
        N0.coef[0][0].elem[l] = val0;
        N0.coef[1][0].elem[l] = val0;
        N0.coef[2][0].elem[l] = val1;
        N0.coef[3][0].elem[l] = val0;

        N0.coef[0][1].elem[l] = val0;
        N0.coef[1][1].elem[l] = val0;
        N0.coef[2][1].elem[l] = val0;
        N0.coef[3][1].elem[l] = val1;
    }
}

// double get_p3(double q1, double q2, double q3, double p1, double p2) //, double mu, double H)
void get_p3(double *x)
{
    // The Hamiltonian is defined in rtbphp.c as shown below
    // double rtbphp_h (int n, int np, void *prm, double x[]) {
    //    double mu=*((double *)prm), xmmu=X-mu, xmmup1=xmmu+1,
    // 	  r12=SQR(xmmu)+SQR(Y)+SQR(Z),
    // 	  r22=SQR(xmmup1)+SQR(Y)+SQR(Z),
    // 	  r1=sqrt(r12), r2=sqrt(r22),
    // 	  p1=(1-mu)/r1, p2=mu/r2;
    //    return .5*(SQR(PX)+SQR(PY)+SQR(PZ))+Y*PX-X*PY-p1-p2;
    // }

    // The Hamiltonian should be contained in lambda[0] and mu in lambda[1]
    double H = lambda[0];
    double mu = lambda[1];
    double q1 = x[0];   double q2 = x[1];   double q3 = x[2];
    double p1 = x[3];   double p2 = x[4];
    double xmmu = q1 - mu, xmmup1 = xmmu + 1,
           r12 = SQR(xmmu) + SQR(q2) + SQR(q3),
           r22 = SQR(xmmup1) + SQR(q2) + SQR(q3),
           r1 = sqrt(r12), r2 = sqrt(r22);
    x[5] = sqrt(2 * (H - q2 * p1 + q1 * p2 + (1 - mu) / r1 + mu / r2) - SQR(p1) - SQR(p2)); // computing pz from the remaining variables and the Hamiltonian
    // return p3;
}

void get_dp3(complex *z, double *dp3){
    double mu = lambda[1];
    double q1 = z[0].real;  double q2 = z[1].real;
    double p1 = z[2].real;  double p2 = z[4].real;
    double xmmu = q1 - mu, xmmup1 = xmmu + 1;
    double r12 = SQR(xmmu) + SQR(q2);
    double r22 = SQR(xmmup1) + SQR(q2);
    double r1 = sqrt(r12);
    double r2 = sqrt(r22);
    double r13 = r1 * r1 * r1;
    double r23 = r2 * r2 * r2;
    dp3[0] = 2 * ( p2 - (1 - mu) * (xmmu) / r13 - mu * (xmmup1) / r23); // dp3dq1
    dp3[1] = 2 * (-p1 - (1 - mu) *   (q2) / r13 - mu *     (q2) / r23); // dp3dq2
    dp3[2] = -2 * q2 - 2 * p1; // dp3dp1
    dp3[3]= 2 * q1 - 2 * p2; // dp3dp2
}

/* nu takes the 4D state and returns the 6D state on the Poincare Map along with the differential*/
void nu(complex *z, double *x, double (*Dnu)[4]){
    x[0] = z[0].real; // q1
    x[1] = z[1].real; // q2
    x[2] = 0;         // q3
    x[3] = z[2].real; // p1
    x[4] = z[3].real; // p2
    get_p3(x); // p3

    /*Compute the differential of the mapping from R4 to Sigma*/
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Dnu[i][j] = 0; // fill it with zeros
        }
    }
    // put in the ones
    Dnu[0][0] = 1;  Dnu[1][1] = 1;  Dnu[3][2] = 1;  Dnu[4][3] = 1;
    
    // get the differential of p3 with respect to the 4D state
    double dp3[4];
    get_dp3(z, dp3); // dp3dq1, dp3dq2, dp3dp1, dp3dp2
    Dnu[5][0] = dp3[0]; Dnu[5][1] = dp3[1]; Dnu[5][2] = dp3[2]; Dnu[5][3] = dp3[3];
}

void map_CR3BP(complex *z, complex *fz, complex **Dfz, complex *depfz)
{
    /* I'm doing a poincare map at a fixed energy level, so I have 4 DOF
    z \in R^4 = (x, z, px, pz)
    I know that y = 0 because we are in the x-z plane (the Poincare map I've defined)
    I can compute py from the Jacobi constant, which is conserved
    */

    //   I need to make a 6D state to propagate and get the poincare map back
    double mu = lambda[1];
    
    /* Map from 4D to 6D space on Poincare section */
    /* Todo function nu()*/
    double x[42];
    double Dnu[6][4];
    nu(z, x, Dnu); // Does the mapping from 4D (z) to 6D (x), the extra 36 states will be filled by the STM
    
    // Fill in the rest of the state with the STM
    for (int i = 6; i < 42; i++)
    {
        x[i] = 0;
    }
    // filling the rest of the array with a vectorized identity matrix
    x[6] = 1;   x[13] = 1;  x[20] = 1;  x[27] = 1;  x[34] = 1;  x[41] = 1;

    // int n = 6;
    // int nv = 42;   // Number of variables in the vector field??? Why is this different than n?, because it could include the STM and go up to 42
    // int np = 0;
    int ibck = 0; // If ibck==1, backward in time (forward if ==0)
    int isiggrad = 1; // needs to be 1 so that the pmap stops in the direction of cp
    double tolJM = 1e-12;
    double maxts = 13;
    // int ivb = 1;
    // int nsecss = 1;
    int nsec = 1; // nsec1 is the number of passes through each section. Since we have only one section, I'm going to make it just an int
    double cp[7] = {0, 0, 1, 0, 0, 0, 0}; // This is 7 dimensional because the last one is the constant term, in case the hyperplane is offset from the origin
    double t = 0;
    double h = fluxvp_pas0;
    if (ibck)
        h = -h;
    FILE *fp = NULL;
    // FILE *fp = fopen("mytraj.txt", "w");
    // FILE
    // setvbuf(fp,NULL,_IOLBF,1024); // Josep-Maria added this so that it would print to the mytraj.txt file correctly (but I just commented it on 7/1/24 because it was causing a segmentation fault after several trajectories)
    /*
     * n : dimension of the initial condition vector
     * nv : number of variables in the vector field
     * np : number of parameters in the vector field
     * camp : vector field
     * prm : parameters of the vector field
     * t : time
     * x : initial condition (gets updated with the final condition)
     * h : step
     * cp : coefficients of the Poincaré section
     * nsec : number of sections to cross
     * isiggrad : if 1, the gradient of the return map is computed
     * tol : tolerance
     * ivb : verbosity level???
     * idt : if 1, the differential of the return map is computed
     * dt : differential of the return map
     * wrtf : function to write the states to a file
     * wrtf_prm : parameters of the function to write the states to a file
     * maxts : maximum time to integrate
     */
    double Dtau[6]; // initialize the derivative of the time of flight with respect to the initial condition (gets filled in within seccp)
    seccp(6 /*n*/, 42 /*nv*/, 0 /*np*/, rtbphp /*camp*/, &mu /*prm*/, &t /*&t*/, x /*x*/, &h /*&h*/, cp /*psec hyperplane*/,
          1 /*nsec*/, isiggrad /*isiggrad*/, tolJM /*tol*/, 0 /*ivb*/, 1 /*idt*/, Dtau /*dt*/, wrtf /*write function*/, fp /*filename*/, 13 /*maxts*/);

    fz[0].real = x[0]; // x
    fz[0].imag = val0;
    fz[1].real = x[1]; // y
    fz[1].imag = val0;
    // x[2], or z, is not part of my state
    fz[2].real = x[3]; // px
    fz[2].imag = val0;
    fz[3].real = x[4]; // py
    fz[3].imag = val0;
    // p3, or x[5], is not part of my state
    /*We fill in DP with the STM, then we'll change it to be the actual differential of the Pmap by adding f_tau*Dtau */
    double DP[6][6];
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            DP[i][j] = x[6 + 6 * i + j]; // fill the rest of the matrix with the STM
        }
    }
    
    // Compute the Jacobian of the map
    // DP = PHI + f_tau*Dtau (this is the differential of the Poincare map)
    // DP = (I - (f_tau*Dsigma/Dsigma*f_tau))*PHI (This is another way of writing it, but I already get Dtau from the seccp function, so the line above is easier)
    double f_tau[6]; // initializes f_tau (the vector field at the final point)
    rtbphp(6 /*n*/, 0 /*np*/, &mu /*prm*/, 0 /*t*/, x /*x*/, f_tau /*dx/dt*/); // fills f_tau with the time derivatives of x (at the final time)

    // double Dsigma[6];
    // for (int i = 0; i < 6; i++)
    // {
    //     Dsigma[i] = cp[i]; // pull from the definition of cp (it's the normal vector to the hyperplane)
    // }

    // double denom = 0;
    // for (int i = 0; i < 6; i++)
    // {
    //     denom += Dsigma[i]*f_tau[i]; // This is actually done inside of seccp function
    // }

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            DP[i][j] = DP[i][j] + f_tau[i]*Dtau[j]; // we add the outer product
        }
    }

    /*Compute the differential of the mapping from Sigma to R4*/
    int DnuInv[4][6];
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            DnuInv[i][j] = 0; // fill it with zeros
        }
    }
    DnuInv[0][0] = 1;  DnuInv[1][1] = 1;  DnuInv[2][3] = 1;  DnuInv[3][4] = 1;
    
    /*Multiply DP and Dnu*/
    double DPDnu[6][4];
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            DPDnu[i][j] = 0;
            for (int k = 0; k < 6; k++)
            {
                DPDnu[i][j] += DP[i][k]*Dnu[k][j];
            }
        }
    }
    /*multiply DnuInv and DPDnu*/
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Dfz[i][j].real = 0;
            Dfz[i][j].imag = 0;
            for (int k = 0; k < 6; k++)
            {
                Dfz[i][j].real += DnuInv[i][k]*DPDnu[k][j];
            }
        }
    }

    // Todo: compute the derivative of the map with respect to the parameter
    for (int i = 0; i < 4; i++)
    {
        depfz[i].real = 0;
        depfz[i].imag = 0;
    }
}

void sform_CR3BP(complex *z, complex **Omegaz)
{
    // fill omega matrix with 2x2 negative identity in the top right corner, and normal identity matrix in the bottom left corner
    Omegaz[0][0] = val0;
    Omegaz[0][1] = val0;
    Omegaz[0][2] = -val1;
    Omegaz[0][3] = val0;

    Omegaz[1][0] = val0;
    Omegaz[1][1] = val0;
    Omegaz[1][2] = val0;
    Omegaz[1][3] = -val1;

    Omegaz[2][0] = val1;
    Omegaz[2][1] = val0;
    Omegaz[2][2] = val0;
    Omegaz[2][3] = val0;

    Omegaz[3][0] = val0;
    Omegaz[3][1] = val1;
    Omegaz[3][2] = val0;
    Omegaz[3][3] = val0;
}

void gform_CR3BP(complex *z, complex **Metricz)
{
    // Do I need to convert all of these to doubles? yes
    double q1, q2, p1, p2;
    q1 = z[0].real;
    q2 = z[1].real;
    p1 = z[2].real;
    p2 = z[3].real;
    double mu = lambda[1];
    double xmmu = q1 - mu, xmmup1 = xmmu + 1;
    double r12 = SQR(xmmu) + SQR(q2);
    double r22 = SQR(xmmup1) + SQR(q2);
    double r1 = sqrt(r12);
    double r2 = sqrt(r22);
    double r13 = r1 * r1 * r1;
    double r23 = r2 * r2 * r2;
    double dp3dq1 = 2 * (p2 - (1 - mu) * (xmmup1) / r13 - mu * (xmmup1) / r23);
    double dp3dq2 = 2 * (-p1 - (1 - mu) * (q2) / r13 - mu * (q2) / r23);
    double dp3dp1 = -2 * q2 - 2 * p1;
    double dp3dp2 = 2 * q1 - 2 * p2;
    // Just a big Identity matrix
    Metricz[0][0] = val1 + SQR(dp3dq1);
    Metricz[0][1] = val0 + dp3dq1 * dp3dq2;
    Metricz[0][2] = val0 + dp3dq1 * dp3dp1;
    Metricz[0][3] = val0 + dp3dq1 * dp3dp2;

    Metricz[1][0] = val0 + dp3dq2 * dp3dq1;
    Metricz[1][1] = val1 + SQR(dp3dq2);
    Metricz[1][2] = val0 + dp3dq2 * dp3dp1;
    Metricz[1][3] = val0 + dp3dq2 * dp3dp2;

    Metricz[2][0] = val0 + dp3dp1 * dp3dq1;
    Metricz[2][1] = val0 + dp3dp1 * dp3dq2;
    Metricz[2][2] = val1 + SQR(dp3dp1);
    Metricz[2][3] = val0 + dp3dp1 * dp3dp2;

    Metricz[3][0] = val0 + dp3dp2 * dp3dq1;
    Metricz[3][1] = val0 + dp3dp2 * dp3dq2;
    Metricz[3][2] = val0 + dp3dp2 * dp3dp1;
    Metricz[3][3] = val1 + SQR(dp3dp2);
}

void normal0_CR3BP(matrix &N0, int *nn, int nelem)
{
    for (int l = 0; l < nelem; l++)
    {
        N0.coef[0][0].elem[l] = val0;
        N0.coef[1][0].elem[l] = val0;
        N0.coef[2][0].elem[l] = val1;
        N0.coef[3][0].elem[l] = val0;

        N0.coef[0][1].elem[l] = val0;
        N0.coef[1][1].elem[l] = val0;
        N0.coef[2][1].elem[l] = val0;
        N0.coef[3][1].elem[l] = val1;
    }
}

void map_standard(complex *z, complex *fz, complex **Dfz, complex *depfz)
{
    fz[1] = z[1] + epsilon * sin(pi2 * z[0]) / pi2;
    fz[0] = z[0] + fz[1];

    Dfz[1][0] = epsilon * cos(pi2 * z[0]);
    Dfz[1][1] = val1;

    Dfz[0][0] = val1 + Dfz[1][0];
    Dfz[0][1] = Dfz[1][1];

    depfz[0] = sin(pi2 * z[0]) / pi2;
    depfz[1] = sin(pi2 * z[0]) / pi2;
}

void sform_standard(complex *z, complex **Omegaz)
{
    Omegaz[0][0] = val0;
    Omegaz[0][1] = -val1;

    Omegaz[1][0] = val1;
    Omegaz[1][1] = val0;
}

void gform_standard(complex *z, complex **Metricz)
{
    Metricz[0][0] = val1;
    Metricz[0][1] = val0;

    Metricz[1][0] = val0;
    Metricz[1][1] = val1;
}

void normal0_standard(matrix &N0, int *nn, int nelem)
{
    for (int l = 0; l < nelem; l++)
    {
        N0.coef[0][0].elem[l] = val0;
        N0.coef[1][0].elem[l] = val1;
    }
}

// int pmap()
// {
//     double mu, t, x[N], h, (*cp)[1 + N], tol, maxts;
//     int nsecss, i, *nsec, isiggrad, ivb, ibck;
//     FILE *fp;
//     /* Precisió de la integració numèrica */
//     fluxvp_pasmin = 1e-6;
//     fluxvp_pasmax = 1000;
//     fluxvp_tol = 1e-14;
//     fluxvp_pas0 = .01;
//     fluxvp_pasminfet = DBL_MAX;
//     fluxvp_pasmaxfet = 0;
//     fluxvp_maxit = 100000;
// /*
//  * Línia de comandes
//     double mu int ibck int isiggrad double tol char* fitxout double maxts int ivb int nsecss
//  * Afegir:
// nsec1 cp1[0..6] nsec2 cp2[0..6] ...\
// - Si ibck==1, endarrere en el temps (si ==0, endavant)\n\
// - Punts per stdin\n\
// - Per no especificar fitxout, s'ha de posar -\n\
// - Per cada punt d'entrada, escriu per stdout temps i punt a la secció\n\
//  */
// // #define FITXOUT argv[5]
// #define FITXOUT '-'

//     // echo -0.9975334497794613 0 -0.00489434320310318 1.383426648384139E-16 -1.002834494433839 -1.991818118201947E-16 | ./rtbp_seccp_main 1.901109735892602e-7=mu 0=ibck 0=isiggrad 1e-12=tol - 13=maxts 1=ivb 1=nsecss 1=nsecc1 0 0 1 0 0 0 0=cp[] > point.txt
//     // By convention, argv[0] is the command with which the program is invoked. argv[1] is the first command-line argument. The last argument from the command line is argv[argc - 1] , and argv[argc] is always NULL.
//     double mu = 1.901109735892602e-7;
//     int ibck = 0; // If ibck==1, backward in time (forward if ==0)
//     int isiggrad = 0;
//     double tol = 1e-12;
//     double maxts = 13;
//     int ivb = 1;
//     int nsecss = 1;
//     int nsec = 1; // nsec1 is the number of passes through each section. Since we have only one section, I'm going to make it just an int
//     double cp[7] = {0, 0, 1, 0, 0, 0, 0};

// #define NARGS 9
//     if (argc < NARGS || sscanf(argv[1], "%lf", &mu) != 1 || sscanf(argv[2], "%d", &ibck) != 1 || sscanf(argv[3], "%d", &isiggrad) != 1 || sscanf(argv[4], "%lf", &tol) != 1 || sscanf(argv[6], "%lf", &maxts) != 1 || sscanf(argv[7], "%d", &ivb) != 1 || sscanf(argv[8], "%d", &nsecss) != 1)
//     {
//         fprintf(stderr, "%s mu ibck isiggrad tol fitxout maxts ivb nsecss \
// nsec1 cp1[0..6] nsec2 cp2[0..6] ...\
// \n\
// - If ibck==1, backward in time (forward if ==0)\n\
// - Points through stdin\n\
// - Use - in order not to specify fitxout\n\
// - For every input point, writes trough stdout return time and point \
//   at the section\n\
// ",
//                 argv[0]);
//         return -1;
//     }
//     /* Fi línia de comandes */
//     /* Fitxer de sortida */
//     if (FITXOUT[0] != '-')
//     {
//         fp = fopen(FITXOUT, "w");
//         if (fp == NULL)
//         {
//             fprintf(stderr, "%s : error obrint fitxout %s !!\n", argv[0],
//                     FITXOUT);
//             return -1;
//         }
//     }
//     else
//         fp = NULL;
//     /* Hiperplans de secció */
//     if (argc < NARGS + nsecss * (2 + N))
//     {
//         fprintf(stderr, "%s : falten nseci o cpi[0..6] !!\n", argv[0]);
//         return -1;
//     }
//     /*This is in case you want another surface of section, it's going to run through the total number of inputs until it's filled in all the surfaces of section*/
//     cp = malloc(nsecss * (N + 1) * sizeof(double));
//     assert(cp != NULL);
//     nsec = malloc(nsecss * sizeof(int));
//     assert(nsec != NULL);
//     for (i = 0; i < nsecss; i++)
//         if (sscanf(argv[NARGS + i * (2 + N)], "%d", &nsec[i]) != 1 || sscanf(argv[1 + NARGS + i * (2 + N)], "%lf", &cp[i][0]) != 1 || sscanf(argv[2 + NARGS + i * (2 + N)], "%lf", &cp[i][1]) != 1 || sscanf(argv[3 + NARGS + i * (2 + N)], "%lf", &cp[i][2]) != 1 || sscanf(argv[4 + NARGS + i * (2 + N)], "%lf", &cp[i][3]) != 1 || sscanf(argv[5 + NARGS + i * (2 + N)], "%lf", &cp[i][4]) != 1 || sscanf(argv[6 + NARGS + i * (2 + N)], "%lf", &cp[i][5]) != 1 || sscanf(argv[7 + NARGS + i * (2 + N)], "%lf", &cp[i][6]) != 1)
//         {
//             fprintf(stderr, "%s : error llegint nsec%d o cp%d[0..6] !!\n",
//                     argv[0], i + 1, i + 1);
//             return -1;
//         }
//     while (scread(stdin, 1, "%lf", &x[0]) == 1)
//     {
//         for (i = 1; i < N; i++)
//             assert(scread(stdin, 1, "%lf", &x[i]) == 1);
//         t = 0;
//         h = fluxvp_pas0;
//         if (ibck)
//             h = -h;
//         for (i = 0; i < nsecss; i++)
//         {
//             if (seccp(N, N /*nv*/, 0 /*np*/, rtbphp /*camp*/, &mu /*prm*/, &t, x, &h, cp[i],
//                       nsec[i], isiggrad, tol, ivb, 0 /*idt*/, NULL /*dt*/, wrtf, fp, maxts)) // Jared 5/8/24 looks like this returns a -1 if it doesn't work, otherwise it returns a 0
//                 fprintf(stderr, "%s : problemes cridant seccp()!!\n", argv[0]);
//             else
//                 printf("%.16G %.16G %.16G %.16G %.16G %.16G %.16G\n",
//                        t, x[0], x[1], x[2], x[3], x[4], x[5]);
//         }
//         if (fp != NULL)
//             fprintf(fp, "\n\n");
//     }
//     return 0;
// }

void wrtf(int n, int nv, double t, double x[], int aon, void *prm)
{
    FILE *fp = (FILE *)prm;
    if (fp != NULL)
        fprintf(fp, "%.16G %.16G %.16G %.16G %.16G %.16G %.16G\n",
                t, x[0], x[1], x[2], x[3], x[4], x[5]);
}

void state2ham(double state[]) // now an unnecessary function
{
    double x, y, z, vx, vy, vz;
    x = state[0];
    y = state[1];
    z = state[2];
    vx = state[3];
    vy = state[4];
    vz = state[5];

    double q1 = x;
    double q2 = y;
    double q3 = z;
    double p1 = vx - y;
    double p2 = vy + x;
    double p3 = vz;

    state[0] = q1;
    state[1] = q2;
    state[2] = q3;
    state[3] = p1;
    state[4] = p2;
    state[5] = p3;
    for (int i = 0; i < 6; i++)
        cout << "state[" << i << "]: " << state[i] << endl;
}
