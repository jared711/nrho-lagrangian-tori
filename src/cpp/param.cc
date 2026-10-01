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
void map_twist(complex *z, complex *fz, complex **Dfz, complex *depfz);
void gform_identity(complex *z, complex **Metricz);
extern double twist_rho[2];
extern double twist_b[2][2];

/* functions created by Jared Blanchard May, 2024*/
void nu(complex *z, double *x, double (*Dnu)[4]);
void get_dp3(complex *z, double p3, double *dp3);
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
int map_failures = 0; // Number of failed Poincare map evaluations since last check
int free_frequency = 0; // 1: Newton corrects omega and fixes the average normal correction (--free-omega)
myreal tolfloor = 1e-10; // accept a stalled Newton iteration if its best error is below this (physical, --tol-floor)
myreal domega_cont[DTOR] = {0.0, 0.0}; // continuation direction in omega per unit epsilon (--domega)
myreal eps_max = 0.01;                  // continuation stops beyond this epsilon (--eps-max)
int max_cont_steps = 10;                // maximum number of continuation steps (--max-steps)
int lowpass_filter = 1;                 // zero the upper half of the DFT after each Newton step (--no-filter)

/* The torus is solved in local coordinates zeta, z = zcen + Mloc zeta, where z = (q1,q2,p1,p2)
   are the physical section coordinates and the columns of Mloc are the first-harmonic axes of
   the input torus, so that the torus is close to a product of round circles in zeta.
   map_CR3BP and gform_CR3BP take zeta; sform_CR3BP returns the constant form
   Omega_loc = Mloc^T Omega Mloc. The defaults (identity) give physical coordinates.
   zscale (the size of the input torus) is only used to convert tolerances and printed errors. */
myreal zcen[DMAP] = {0.0, 0.0, 0.0, 0.0};
myreal zscale = 1.0;
myreal Mloc[DMAP][DMAP] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
myreal Minv[DMAP][DMAP] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
myreal Omega_loc[DMAP][DMAP] = {{0, 0, -1, 0}, {0, 0, 0, -1}, {1, 0, 0, 0}, {0, 1, 0, 0}};
void to_physical(complex *zeta, complex *z);
int invert4(myreal A[DMAP][DMAP], myreal Ainv[DMAP][DMAP]);

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

    for (int i = 2; i < argc; i++)
        if (strcmp(argv[i], "--free-omega") == 0)
            free_frequency = 1;
        else if (strcmp(argv[i], "--tol-floor") == 0 && i + 1 < argc)
            tolfloor = atof(argv[++i]);
        else if (strcmp(argv[i], "--domega") == 0 && i + 2 < argc)
        {
            domega_cont[0] = atof(argv[++i]);
            domega_cont[1] = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--eps-max") == 0 && i + 1 < argc)
            eps_max = atof(argv[++i]);
        else if (strcmp(argv[i], "--max-steps") == 0 && i + 1 < argc)
            max_cont_steps = atoi(argv[++i]);
        else if (strcmp(argv[i], "--no-filter") == 0)
            lowpass_filter = 0;

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

    /* Orbit mode: iterate the Poincare map from z = (q1,q2,p1,p2) given on the command line.
       Usage: param input.csv --orbit niter q1 q2 p1 p2
       Prints each iterate; the first line also contains the Jacobian of the map at z. */
    if (argc > 7 && strcmp(argv[2], "--orbit") == 0)
    {
        int niter = atoi(argv[3]);
        complex zc[DMAP], f0[DMAP], dep[DMAP];
        complex *D0[DMAP];
        for (int i = 0; i < DMAP; i++)
        {
            D0[i] = new complex[DMAP];
            zc[i] = atof(argv[4 + i]);
        }
        cout.precision(16);
        for (int k = 0; k < niter; k++)
        {
            map_CR3BP(zc, f0, D0, dep);
            if (map_failures > 0)
            {
                cout << "# map failed at iterate " << k << endl;
                return 1;
            }
            if (k == 0)
            {
                cout << "# Dfz";
                for (int i = 0; i < DMAP; i++)
                    for (int j = 0; j < DMAP; j++)
                        cout << " " << D0[i][j].real;
                cout << endl;
            }
            for (int i = 0; i < DMAP; i++)
            {
                zc[i] = f0[i];
                cout << zc[i].real << (i < DMAP - 1 ? " " : "\n");
            }
        }
        return 0;
    }

    /* Validation mode: Newton iterations of kam_torus on the integrable twist map map_twist,
       starting from a perturbed exact torus. Uses the grid size and tolerances of the input file.
       Usage: param input.csv --test-twist [Case] */
    if (argc > 2 && strcmp(argv[2], "--test-twist") == 0)
    {
        int tcase = (argc > 3) ? atoi(argv[3]) : 2;
        double Istar[2] = {0.01, 0.02};
        for (int j = 0; j < 2; j++)
            omega[j] = twist_rho[j] + twist_b[j][0] * Istar[0] + twist_b[j][1] * Istar[1];
        cout << "# twist test: omega = " << omega[0] << " " << omega[1] << ", Case " << tcase << endl;
        cout << "# exact <T> = 2 pi b (up to the frame normalization)" << endl;
        for (int l = 0; l < nelem; l++)
        {
            indices(l, nn, index, DTOR);
            double th[2] = {(double)index[0] / nn[0], (double)index[1] / nn[1]};
            for (int j = 0; j < 2; j++)
            {
                double r = sqrt(2 * Istar[j] * 1.01); // 1% wrong action
                paramR.coef[j][0].elem[l] = r * cos(pi2 * th[j]);
                paramR.coef[j + 2][0].elem[l] = r * sin(pi2 * th[j]);
            }
            paramR.coef[0][0].elem[l] = paramR.coef[0][0].elem[l].real + 1e-4 * cos(pi2 * (th[0] + th[1]));
        }
        paramF = fft_F(paramR);
        for (int it = 0; it < 8; it++)
        {
            cout << "# Newton iteration " << it + 1 << endl;
            conv = kam_torus(paramR, paramF, omega, error, nn, nelem, tail0, tails, tcase,
                             map_twist, sform_CR3BP, gform_identity, normal0_CR3BP);
            if (conv != 0 || tail0 == 1)
                break;
        }
        return 0;
    }

    /* Evaluation mode: print F(K) and DF(K) at every grid point of the input torus
       (physical coordinates, grid order of indices()). One line per point:
       l F_0..F_3 DF_00 DF_01 .. DF_33 (row-major). Usage: param input.csv --eval */
    if (argc > 2 && strcmp(argv[2], "--eval") == 0)
    {
        complex zc[DMAP], f0[DMAP], dep[DMAP];
        complex *D0[DMAP];
        for (int i = 0; i < DMAP; i++)
            D0[i] = new complex[DMAP];
        cout.precision(17);
        for (int l = 0; l < nelem; l++)
        {
            for (int i = 0; i < DMAP; i++)
                zc[i] = paramR.coef[i][0].elem[l];
            map_CR3BP(zc, f0, D0, dep);
            cout << l;
            for (int i = 0; i < DMAP; i++)
                cout << " " << f0[i].real;
            for (int i = 0; i < DMAP; i++)
                for (int j = 0; j < DMAP; j++)
                    cout << " " << D0[i][j].real;
            cout << endl;
        }
        if (map_failures > 0)
            cout << "# map failures: " << map_failures << endl;
        return 0;
    }

    /* Verification mode: compare Dfz from map_CR3BP against central finite differences */
    if (argc > 2 && strcmp(argv[2], "--fdcheck") == 0)
    {
        complex zc[DMAP], zp[DMAP], zm[DMAP], f0[DMAP], fp[DMAP], fm[DMAP], dep[DMAP];
        complex *D0[DMAP], *Dtmp[DMAP];
        for (int i = 0; i < DMAP; i++)
        {
            D0[i] = new complex[DMAP];
            Dtmp[i] = new complex[DMAP];
            zc[i] = paramR.coef[i][0].elem[0];
        }
        map_CR3BP(zc, f0, D0, dep);
        double hfd = 1e-7, maxrel = 0;
        for (int j = 0; j < DMAP; j++)
        {
            for (int i = 0; i < DMAP; i++)
                zp[i] = zm[i] = zc[i];
            zp[j].real += hfd;
            zm[j].real -= hfd;
            map_CR3BP(zp, fp, Dtmp, dep);
            map_CR3BP(zm, fm, Dtmp, dep);
            for (int i = 0; i < DMAP; i++)
            {
                double fd = (fp[i].real - fm[i].real) / (2 * hfd);
                double an = D0[i][j].real;
                double rel = fabs(fd - an) / fmax(1.0, fabs(fd));
                maxrel = fmax(maxrel, rel);
                cout << "Dfz[" << i << "][" << j << "] analytic " << an << " fd " << fd << " diff " << rel << endl;
            }
        }
        cout << "# max |analytic - fd| / max(1,|fd|) = " << maxrel << endl;
        cout << "# map failures: " << map_failures << endl;
        return 0;
    }
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

    /* Change to local coordinates z = zcen + Mloc zeta, centered at the average of the torus,
       with Mloc built from its first harmonics. In zeta the torus is close to a product of round
       circles, so L^T L is nearly constant and the normal frame N is well resolved on the grid.
       The change is linear, so the form stays constant: Omega_loc = Mloc^T Omega Mloc. */
    for (int i = 0; i < DMAP; i++)
        zcen[i] = paramF.coef[i][0].elem[0].real;
    zscale = 0.0;
    for (int i = 0; i < DMAP; i++)
        for (int l = 0; l < nelem; l++)
            zscale = fmax(zscale, fabs(paramR.coef[i][0].elem[l].real - zcen[i]));
    {
        /* First harmonics: K(theta) ~ zcen + sum_j a_j cos(2 pi theta_j) + b_j sin(2 pi theta_j),
           with a_j = 2 Re K_{e_j}, b_j = -2 Im K_{e_j}. Columns of Mloc: a_1 a_2 b_1 b_2. */
        int pos1 = nn[1]; // index (1,0): position 1 * nn[1] + 0
        int pos2 = 1;     // index (0,1)
        for (int i = 0; i < DMAP; i++)
        {
            Mloc[i][0] = 2.0 * paramF.coef[i][0].elem[pos1].real;
            Mloc[i][1] = 2.0 * paramF.coef[i][0].elem[pos2].real;
            Mloc[i][2] = -2.0 * paramF.coef[i][0].elem[pos1].imag;
            Mloc[i][3] = -2.0 * paramF.coef[i][0].elem[pos2].imag;
        }
        /* Symplectic normalization of each pair (a_j, b_j): scale so that a_j^T Omega b_j = -1, the
           value of e_q^T Omega e_p for the standard form. Omega_loc is then close to the standard
           form. Afterwards all columns are multiplied by a common factor sym_scale = sqrt(mean action)
           so that the circles of the torus have radius ~1 in zeta and every quantity in kam_torus
           (L, N, T, eta) is O(1): the grid utilities use absolute tolerances (tolgrid in
           cohomological() and clean(), tolqr in inv()), which fail for a torus of size ~1e-5.
           The form is divided by sym_scale^2 accordingly (a constant multiple of a symplectic form
           is preserved by the same maps). */
        double sym_scale = 0.0;
        {
            double Om0[DMAP][DMAP] = {{0, 0, -1, 0}, {0, 0, 0, -1}, {1, 0, 0, 0}, {0, 1, 0, 0}};
            for (int j = 0; j < DTOR; j++)
            {
                double w = 0.0;
                for (int k = 0; k < DMAP; k++)
                    for (int m = 0; m < DMAP; m++)
                        w += Mloc[k][j] * Om0[k][m] * Mloc[m][j + DTOR];
                if (w == 0.0)
                    continue;
                sym_scale += fabs(w) / DTOR;
                double sa = 1.0 / sqrt(fabs(w)), sb = (w > 0) ? -sa : sa;
                for (int k = 0; k < DMAP; k++)
                {
                    Mloc[k][j] *= sa;
                    Mloc[k][j + DTOR] *= sb;
                }
            }
            sym_scale = sqrt(sym_scale);
            if (sym_scale > 0.0)
                for (int k = 0; k < DMAP; k++)
                    for (int j = 0; j < DMAP; j++)
                        Mloc[k][j] *= sym_scale;
            else
                sym_scale = 1.0;
        }
        if (!invert4(Mloc, Minv))
        {
            cout << "# First-harmonic matrix is singular; using isotropic scaling" << endl;
            for (int i = 0; i < DMAP; i++)
                for (int j = 0; j < DMAP; j++)
                {
                    Mloc[i][j] = (i == j) ? zscale : 0.0;
                    Minv[i][j] = (i == j) ? 1.0 / zscale : 0.0;
                }
        }
        double Om[DMAP][DMAP] = {{0, 0, -1, 0}, {0, 0, 0, -1}, {1, 0, 0, 0}, {0, 1, 0, 0}};
        for (int i = 0; i < DMAP; i++)
            for (int j = 0; j < DMAP; j++)
            {
                Omega_loc[i][j] = 0.0;
                for (int k = 0; k < DMAP; k++)
                    for (int m = 0; m < DMAP; m++)
                        Omega_loc[i][j] += Mloc[k][i] * Om[k][m] * Mloc[m][j];
                Omega_loc[i][j] /= sym_scale * sym_scale;
            }
        /* zscale converts local lengths to physical ones: |dz| <~ zscale |dzeta| */
        zscale = 0.0;
        for (int j = 0; j < DMAP; j++)
        {
            double c = 0.0;
            for (int i = 0; i < DMAP; i++)
                c += Mloc[i][j] * Mloc[i][j];
            zscale = fmax(zscale, sqrt(c));
        }
        cout << "# Omega_loc:";
        for (int i = 0; i < DMAP; i++)
            for (int j = 0; j < DMAP; j++)
                cout << " " << Omega_loc[i][j];
        cout << endl;
    }
    for (int l = 0; l < nelem; l++)
    {
        double dz[DMAP];
        for (int i = 0; i < DMAP; i++)
            dz[i] = paramR.coef[i][0].elem[l].real - zcen[i];
        for (int i = 0; i < DMAP; i++)
        {
            double v = 0.0;
            for (int j = 0; j < DMAP; j++)
                v += Minv[i][j] * dz[j];
            paramR.coef[i][0].elem[l] = v;
        }
    }
    paramF = fft_F(paramR);
    toltail /= zscale; // tolerances in the input file are in physical units
    tolinva /= zscale;
    cout << "# zcen: " << zcen[0] << " " << zcen[1] << " " << zcen[2] << " " << zcen[3] << endl;
    cout << "# zscale: " << zscale << endl;

    /* No Newton step here: the Newton loop below starts from the input torus. A step taken here
       (as in the original code) would be stored in paramR0 below, so every continuation step and
       every non-converged restart would start from the once-corrected torus instead of the input. */

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
    /* Continuation in the rotation vector at fixed energy (Haro & Mondelo 2021, arXiv:2101.07665,
       Sec. 3.5, Alg. 3.5.4): omega(epsilon) = omega_start + (epsilon - epsilon_start) * domega, with
       domega given by --domega (cycles per iterate per unit epsilon). Without --domega, epsilon does
       not enter the problem (the map does not depend on it). Predictor: secant through the last two
       converged tori. Step control as in Haro & Mondelo, Alg. 3.6.1, pp. 29-30: halve the step on
       failure, and set deps <- deps * n_des / n_it after a success (factor clamped to [0.5, 2]). */
    int fail = 0;
    int cont_step = 0;
    myreal omega_start[DTOR], epsilon_start = epsilon0, deps_start = deps, deps_last = 0.0;
    for (int i = 0; i < DTOR; i++)
        omega_start[i] = omega[i];
    matrix paramRm1, paramFm1; // previous converged torus, for the secant predictor
    int have_prev = 0;
    const int N_DES = 4;
    do
    {
        paramR = paramR0;
        paramF = paramF0;
        epsilon = (cont_step == 0) ? epsilon0 : epsilon0 + deps; // first converge the input torus
        for (int i = 0; i < DTOR; i++)
            omega[i] = omega_start[i] + (epsilon - epsilon_start) * domega_cont[i];
        if (have_prev && deps_last > 0.0)
        {
            paramR = paramR0 + (paramR0 - paramRm1) * (deps / deps_last);
            paramF = fft_F(paramR);
        }

        // Check continuation termination criteria
        if ((deps > 0 && epsilon > eps_max) || (deps < 0 && epsilon < eps_max)) {
            cout << "# Reached maximum epsilon = " << epsilon << " (limit: " << eps_max << ")" << endl;
            cout << "# Terminating continuation." << endl;
            fail = 1;
            break;
        }
        if (cont_step >= max_cont_steps) {
            cout << "# Reached maximum continuation steps = " << cont_step << " (limit: " << max_cont_steps << ")" << endl;
            cout << "# Terminating continuation." << endl;
            fail = 1;
            break;
        }
        cont_step++;

        /**** START Newton method to correct the invariant torus ****/
        cout << "# Continuation step " << cont_step << " / " << max_cont_steps << endl;
        cout << "# We try to compute the torus for epsilon=" << epsilon << " (limit: " << eps_max << ")" << endl;
        iter = 0;
        /* Stall detection: the error of the map evaluation sets a floor below which Newton cannot go;
           once there, further steps only amplify that noise through the small divisors. We keep the
           best torus seen and stop when a step fails to reduce the error by 10%. The best torus is accepted
           as converged if its error is below tolfloor (physical units, --tol-floor). */
        double best_error = DBL_MAX, prev_error = DBL_MAX;
        matrix bestR, bestF;
        myreal best_omega[DTOR];
        do
        {
            cout << "# Iteration " << iter + 1 << " : " << endl;
            matrix entryR = paramR, entryF = paramF;
            myreal entry_omega[DTOR];
            for (int i = 0; i < DTOR; i++)
                entry_omega[i] = omega[i];
            conv = kam_torus(paramR, paramF, omega, error, nn, nelem, tail0, tails, 2, map_CR3BP, sform_CR3BP, gform_identity, normal0_CR3BP); // Case 2: Case 1 (constant N0) is not transversal for tori around an elliptic point
            if (conv == 0 && tail0 == 0)
            {
                if (error < best_error)
                {
                    best_error = error;
                    bestR = entryR;
                    bestF = entryF;
                    for (int i = 0; i < DTOR; i++)
                        best_omega[i] = entry_omega[i];
                }
                if (error > 0.9 * prev_error)
                {
                    paramR = bestR;
                    paramF = bestF;
                    for (int i = 0; i < DTOR; i++)
                        omega[i] = best_omega[i];
                    error = best_error;
                    conv = (best_error * zscale < tolfloor) ? 1 : 0;
                    cout << "# Newton stalled; best error " << best_error * zscale << " (physical), "
                         << (conv ? "accepted" : "not accepted") << " (tol-floor " << tolfloor << ")" << endl;
                    if (conv == 0)
                        conv = -1;
                    break;
                }
                prev_error = error;
            }
            if (tail0 == 1)
                best_error = prev_error = DBL_MAX; // the grid changes: restart the stall detection
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
                have_prev = 0; // the previous torus is on the old grid
            }
            iter++;
        } while (conv == 0 && iter < MNEW && tail0 == 0);
        /**** END   Newton method to correct the invariant torus ****/

        if (conv == 1)
        {
            cout.precision(15);
            /* If we converge, we update the last computed torus and store the results */
            if (paramR0.coef[0][0].nelem == paramR.coef[0][0].nelem)
            {
                paramRm1 = paramR0;
                paramFm1 = paramF0;
                have_prev = 1;
            }
            else
                have_prev = 0;
            deps_last = deps;
            paramR0 = paramR;
            paramF0 = paramF;
            epsilon0 = epsilon;
            double factor = (double)N_DES / (double)(iter > 0 ? iter : 1);
            deps *= fmin(2.0, fmax(0.5, factor));

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
            sprintf(name, "output_torus%.6lf", epsilon);
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
                complex zeta[DMAP], zphys[DMAP];
                for (int j = 0; j < DMAP; j++)
                    zeta[j] = paramR.coef[j][0].elem[l];
                to_physical(zeta, zphys);
                for (int j = 0; j < DMAP - 1; j++)
                    file_torus << zphys[j].real << " ";
                file_torus << zphys[DMAP - 1].real << endl;
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
                cout << "# Newton method does not converge; halving the continuation step" << endl;
                deps = deps / val2;
                have_prev = (deps_last > 0.0) ? have_prev : 0;
                if (fabs(deps) < 1e-3 * fabs(deps_start))
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

static int nyquist_index(int *nn, int *index)
{
    for (int j = 0; j < DTOR; j++)
        if (2 * index[j] == nn[j])
            return 1;
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
        if (nyquist_index(nn, index)) // has no conjugate partner; copying it one-sided makes K complex
            continue;
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
        if (nyquist_index(nn, index)) // has no conjugate partner; copying it one-sided makes K complex
            continue;
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

    if (map_failures > 0)
    {
        cout << "#     - Poincare map failed at " << map_failures << " grid points" << endl;
        map_failures = 0;
        for (int i = 0; i < DTOR; i++)
            tails[i] = 0;
        return -1;
    }

    KshiftF = shift(paramF, omega);
    KshiftR = fft_B(KshiftF);

    /* The torus lives in Cartesian coordinates (q1,q2,p1,p2), so K is fully
       periodic: no theta+omega lift is added to the first DTOR components
       (that lift is only for tori in angle-action coordinates, e.g. Froeschle). */

    ErrorR = FparamR - KshiftR;

    ErrorF = fft_F(ErrorR);
    error = norm(ErrorF);
    cout << "#     - Error of invariance: ";
    cout << error * zscale << " (physical), " << error << " (local)" << endl;

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
    if (Case == 2)
    {
        /* N(theta + omega) from the frame formula applied pointwise to L(theta + omega), the metric
           and the form at K(theta + omega), instead of Fourier-shifting N: N contains
           B = (L^T G L)^-1, whose spectrum decays more slowly than that of L, so shifting it
           on the grid is less accurate than shifting L. */
        matrix MetricKshiftR = fft_B(shift(fft_F(MetricKR), omega));
        matrix GRs = trans(LshiftR) * MetricKshiftR * LshiftR;
        matrix BRs = inv(GRs);
        matrix ARs = trans(BRs) * trans(LshiftR) * MetricKshiftR * inv(OmegaKshiftR) * MetricKshiftR * LshiftR * BRs * val05;
        NshiftR = LshiftR * ARs - inv(OmegaKshiftR) * MetricKshiftR * LshiftR * BRs;
    }
    etaLR = -trans(NshiftR) * OmegaKshiftR * ErrorR;
    etaNR = trans(LshiftR) * OmegaKshiftR * ErrorR;
    twistR = trans(NshiftR) * OmegaKshiftR * DFKR * NR;
    etaNF = fft_F(etaNR);
    RetaNF = cohomological(etaNF, omega);
    RetaNR = fft_B(RetaNF);
    newetaR = etaLR - twistR * RetaNR;

    aver(twistR, twist0);
    cout << "#     - Average torsion <T>: [" << twist0[0][0] << " " << twist0[0][1] << "; "
         << twist0[1][0] << " " << twist0[1][1] << "], det " << twist0[0][0] * twist0[1][1] - twist0[0][1] * twist0[1][0] << endl;
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

    /* Free-frequency variant: fix the average of xi_N to zero (the "action" of the torus)
       and correct the frequency instead. With K(theta + omega + domega) in the invariance
       error, the tangent row of the reduced equation gains -domega, so its solvability
       condition <eta_L - T xi_N> + domega = 0 gives domega = -<eta_L - T R(eta_N)> = -neweta0.
       No inverse of the torsion <T> is needed. */
    myreal domega[DTOR];
    if (free_frequency)
    {
        for (int i = 0; i < DTOR; i++)
        {
            xiN0[i][0] = 0.0;
            domega[i] = -neweta0[i][0];
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
    if (free_frequency)
        for (int l = 0; l < nelem; l++)
            for (int i = 0; i < DTOR; i++)
                newetaR.coef[i][0].elem[l] = newetaR.coef[i][0].elem[l] + domega[i];
    newetaF = fft_F(newetaR);
    xiLF = cohomological(newetaF, omega);
    xiLR = fft_B(xiLF);

    /*****************************************************************
    STEP 4 New parameterization
    *****************************************************************/
    newparamR = paramR + LR * xiLR + NR * xiNR;
    newparamF = fft_F(newparamR);
    if (lowpass_filter)
    {
        /* Zero the upper half of the DFT coefficients (|k_j| >= nn[j]/4 in some direction) after
           each Newton step, as in Haro & Mondelo 2021 (arXiv:2101.07665, pp. 27-28), to prevent
           "divergence after apparent convergence": the relative error of the DFT coefficients
           grows with |k|, and the small divisors amplify it. Same band as Step 0 of Figueras,
           Haro & Luque 2017 (arXiv:1601.00084, Sec. 5.2). The tail test then measures how much
           the step tried to put in that band. */
        for (int l = 0; l < nelem; l++)
        {
            int idx[DTOR], ser[DTOR];
            indices(l, nn, idx, DTOR);
            trigo_to_series(nn, idx, ser, DTOR);
            int high = 0;
            for (int j = 0; j < DTOR; j++)
                if (4 * abs(ser[j]) >= nn[j])
                    high = 1;
            if (high)
                for (int i = 0; i < DMAP; i++)
                    newparamF.coef[i][0].elem[l] = 0.0;
        }
        newparamR = fft_B(newparamF);
    }

    /* Diagnostic: residual of the linearized equation DF(K) dK - dK(theta+omega) = -E */
    {
        matrix dKR = LR * xiLR + NR * xiNR;
        matrix dKshiftR = fft_B(shift(fft_F(dKR), omega));
        matrix linres = DFKR * dKR - dKshiftR + ErrorR;
        cout << "#     - Residual of the linearized equation: " << norm(fft_F(linres)) * zscale << " (physical)" << endl;
        matrix symp = trans(LR) * OmegaKR * LR;
        cout << "#     - Lagrangian defect |L^T Omega L|: " << norm(fft_F(symp)) << endl;
        matrix frame = trans(LR) * OmegaKR * NR;
        cout << "#     - Frame check |L^T Omega N + I| (should be 0): ";
        for (int l = 0; l < nelem; l++)
            for (int i = 0; i < DTOR; i++)
                frame.coef[i][i].elem[l] = frame.coef[i][i].elem[l] + val1;
        cout << norm(fft_F(frame)) << endl;
    }

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
    if (free_frequency)
    {
        for (int i = 0; i < DTOR; i++)
            omega[i] = omega[i] + domega[i];
        cout.precision(15);
        cout << "#     - Frequency correction: " << domega[0] << " " << domega[1] << ", new omega: " << omega[0] << " " << omega[1] << endl;
        cout.precision(3);
    }
    paramF = newparamF;

    /* This is just to show the size of the correction */
    newparamR = LR * xiLR + NR * xiNR;
    newparamF = fft_F(newparamR);
    aux = norm(newparamF);
    cout << "#     - Norm of the correction: ";
    cout << aux * zscale << " (physical), " << aux << " (local)" << endl;

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

/* Gradient of p3(q1,q2,p1,p2) on the section q3=0 at fixed energy H.
   d(p3^2) is computed first and then divided by 2*p3. */
void get_dp3(complex *z, double p3, double *dp3){
    double mu = lambda[1];
    double q1 = z[0].real;  double q2 = z[1].real;
    double p1 = z[2].real;  double p2 = z[3].real;
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
    for (int i = 0; i < 4; i++)
        dp3[i] /= 2 * p3;
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
    get_dp3(z, x[5], dp3); // dp3dq1, dp3dq2, dp3dp1, dp3dp2
    Dnu[5][0] = dp3[0]; Dnu[5][1] = dp3[1]; Dnu[5][2] = dp3[2]; Dnu[5][3] = dp3[3];
}

void to_physical(complex *zeta, complex *z)
{
    for (int i = 0; i < 4; i++)
    {
        double v = zcen[i];
        for (int j = 0; j < 4; j++)
            v += Mloc[i][j] * zeta[j].real;
        z[i] = complex(v, 0.0);
    }
}

/* Inverse of a 4x4 matrix by Gauss-Jordan elimination with partial pivoting. Returns 0 if singular. */
int invert4(myreal A[DMAP][DMAP], myreal Ainv[DMAP][DMAP])
{
    double a[DMAP][2 * DMAP];
    for (int i = 0; i < DMAP; i++)
        for (int j = 0; j < DMAP; j++)
        {
            a[i][j] = A[i][j];
            a[i][j + DMAP] = (i == j) ? 1.0 : 0.0;
        }
    for (int c = 0; c < DMAP; c++)
    {
        int p = c;
        for (int r = c + 1; r < DMAP; r++)
            if (fabs(a[r][c]) > fabs(a[p][c]))
                p = r;
        if (fabs(a[p][c]) < 1e-300)
            return 0;
        for (int j = 0; j < 2 * DMAP; j++)
        {
            double t = a[c][j];
            a[c][j] = a[p][j];
            a[p][j] = t;
        }
        double d = a[c][c];
        for (int j = 0; j < 2 * DMAP; j++)
            a[c][j] /= d;
        for (int r = 0; r < DMAP; r++)
            if (r != c)
            {
                double m = a[r][c];
                for (int j = 0; j < 2 * DMAP; j++)
                    a[r][j] -= m * a[c][j];
            }
    }
    for (int i = 0; i < DMAP; i++)
        for (int j = 0; j < DMAP; j++)
            Ainv[i][j] = a[i][j + DMAP];
    return 1;
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
    complex zphys[4];
    to_physical(z, zphys); // z is in local coordinates (see zcen, zscale)
    nu(zphys, x, Dnu); // Does the mapping from 4D (z) to 6D (x), the extra 36 states will be filled by the STM
    
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
    double tolJM = 1e-14; // section-crossing tolerance; 1e-12 left a map-noise floor of ~2e-11 (see commit message)
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
    int sec_ret = seccp(6 /*n*/, 42 /*nv*/, 0 /*np*/, rtbphp /*camp*/, &mu /*prm*/, &t /*&t*/, x /*x*/, &h /*&h*/, cp /*psec hyperplane*/,
          1 /*nsec*/, isiggrad /*isiggrad*/, tolJM /*tol*/, 0 /*ivb*/, 1 /*idt*/, Dtau /*dt*/, wrtf /*write function*/, fp /*filename*/, maxts /*maxts*/);
    if (sec_ret < 0 || x[5] != x[5]) // seccp failed, or p3 is NaN (point is off the energy surface)
        map_failures++;

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
    { // back to local coordinates: fz = Minv (F - zcen)
        double dz[4];
        for (int i = 0; i < 4; i++)
            dz[i] = fz[i].real - zcen[i];
        for (int i = 0; i < 4; i++)
        {
            fz[i].real = 0.0;
            for (int j = 0; j < 4; j++)
                fz[i].real += Minv[i][j] * dz[j];
        }
    }
    /*We fill in DP with the STM, then we'll change it to be the actual differential of the Pmap by adding f_tau*Dtau */
    double DP[6][6];
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            DP[i][j] = *vr1(6, 0, x, i, j); // STM is stored column-major: x[6 + 6*j + i] = Phi_ij
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
    { // Jacobian in local coordinates: Minv DF Mloc
        double DM[4][4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                DM[i][j] = 0.0;
                for (int k = 0; k < 4; k++)
                    DM[i][j] += Dfz[i][k].real * Mloc[k][j];
            }
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                double v = 0.0;
                for (int k = 0; k < 4; k++)
                    v += Minv[i][k] * DM[k][j];
                Dfz[i][j] = complex(v, 0.0);
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
    // Constant symplectic form in the torus coordinates: Omega_loc = Mloc^T Omega Mloc, where
    // Omega has the 2x2 negative identity in the top right corner and the identity in the bottom
    // left corner (Omega_loc = Omega in physical coordinates)
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            Omegaz[i][j] = Omega_loc[i][j];
}

void gform_CR3BP(complex *z, complex **Metricz)
{
    // Metric induced on the section by the Euclidean metric of R^6: I + grad(p3) grad(p3)^T
    double x[6], Dnu[6][4], dp3[4];
    complex zphys[4];
    to_physical(z, zphys);
    nu(zphys, x, Dnu);
    get_dp3(zphys, x[5], dp3);
    double dp3dq1 = dp3[0], dp3dq2 = dp3[1], dp3dp1 = dp3[2], dp3dp2 = dp3[3];
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

/* Integrable 4D symplectic twist map, used to validate kam_torus on Cartesian (non-lifted) tori.
   In each plane (q_j, p_j) it rotates by phi_j = 2 pi (rho_j + sum_k b_jk I_k), with
   I_k = (q_k^2 + p_k^2) / 2; it is the time-1 map of H = 2 pi rho.I + pi I^T b I.
   The torus of frequency omega is K(theta)_j = sqrt(2 I*_j) (cos 2 pi theta_j, sin 2 pi theta_j)
   with rho + b I* = omega, and its averaged torsion is known in closed form. */
double twist_rho[2] = {0.3435, 0.1733};
double twist_b[2][2] = {{1.0, 0.3}, {0.3, 0.5}};

void map_twist(complex *z, complex *fz, complex **Dfz, complex *depfz)
{
    double q[2] = {z[0].real, z[1].real}, p[2] = {z[2].real, z[3].real}, I[2], phi[2];
    for (int k = 0; k < 2; k++)
        I[k] = 0.5 * (q[k] * q[k] + p[k] * p[k]);
    for (int j = 0; j < 2; j++)
        phi[j] = pi2 * (twist_rho[j] + twist_b[j][0] * I[0] + twist_b[j][1] * I[1]);
    for (int i = 0; i < 4; i++)
        for (int k = 0; k < 4; k++)
            Dfz[i][k] = val0;
    for (int j = 0; j < 2; j++)
    {
        double c = cos(phi[j]), s = sin(phi[j]);
        fz[j] = c * q[j] - s * p[j];     // q_j'
        fz[j + 2] = s * q[j] + c * p[j]; // p_j'
        // rotation part
        Dfz[j][j] = c;
        Dfz[j][j + 2] = -s;
        Dfz[j + 2][j] = s;
        Dfz[j + 2][j + 2] = c;
        // dependence of phi_j on (q_k, p_k): dphi_j/dq_k = 2 pi b_jk q_k, dphi_j/dp_k = 2 pi b_jk p_k
        double dq = -s * q[j] - c * p[j], dp = c * q[j] - s * p[j]; // d(q_j', p_j')/dphi_j
        for (int k = 0; k < 2; k++)
        {
            double a = pi2 * twist_b[j][k];
            Dfz[j][k] = Dfz[j][k].real + dq * a * q[k];
            Dfz[j][k + 2] = Dfz[j][k + 2].real + dq * a * p[k];
            Dfz[j + 2][k] = Dfz[j + 2][k].real + dp * a * q[k];
            Dfz[j + 2][k + 2] = Dfz[j + 2][k + 2].real + dp * a * p[k];
        }
    }
    for (int i = 0; i < 4; i++)
        depfz[i] = val0;
}

void gform_identity(complex *z, complex **Metricz)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            Metricz[i][j] = (i == j) ? val1 : val0;
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
