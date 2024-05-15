/*-
 * \file campvp.c
 *
 * \brief Convenis de camp vectorial amb variacionals 1eres respecte
 * de paràmetres (no es pot estendre a ordre superior) i 1eres,2ones 
 * respecte de c.i. (es pot estendre a ordre superior).
 */

double *vrp(int n, double *x, int i, int j) {
   return x+n*(1+j)+i;
}

double *vr1(int n, int np, double *x, int i, int j) {
   return x+n*(1+np+j)+i;
}

/* Atenció!!! Se suposa j>=i !! */

double *vr2(int n, int np, double *x, int i, int j, int k) {
   return x+((1+np)*n+(n*n))+(i*(2*n-i+1)/2+j-i)*n+k;
}

