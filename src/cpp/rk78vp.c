#include <stdio.h>
#include <string.h>
#include <math.h>
#include "campvp.h"
#include "utils.h"
#include "vbprintf.h"

/*!
 * \brief Versió rk78() amb suport per variacionals 1eres i 2ones,
 * variacionals respecte de paràmetres, i també per a camp
 * vectorial no evaluable (petar amb l'esfera de convergència).
 * \c nn és la "dimensió espaial", i es fa servir per a calcular
 * el nombre d'equacions diferencials com a <tt>n=nv+np*nn</tt>.
 * D'aquí que \c nv inclou equacions originals i variacionals (42 per
 * al RTBP amb variacionals), i np és el nombre de paràmetres resp.
 * del qual es volen trobar variacionals (1 per al RTBP, si es vol
 * trobar derivada resp mu).
 */

int rk78vp (int nn, int nv, int np, void *prm, campvp_t camp,
      double *t, double x[], double *h, double hmin, double hmax, double tol,
      double *eerr, int ivb) {
/* Coeficients RK (estatics) */
   static double alfa[13] = {
	      0.,     2./27.,      1./9.,      1./6.,     5./12.,
	      .5,      5./6.,      1./6.,      2./3.,      1./3.,
	      1.,         0.,         1.
   };
   static double beta[78] = {
	              2./27.,     1./36.,     1./12.,     1./24.,
	      0.,      1./8.,     5./12.,         0.,   -25./16.,
	 25./16.,        .05,         0.,         0.,        .25,
	      .2,  -25./108.,         0.,         0.,  125./108.,
	-65./27., 2.*(125./108.), 31./300.,       0.,         0.,
	      0.,   61./225.,     -2./9.,   13./900.,         2.,
	      0.,         0.,    -53./6.,   704./45.,   -107./9.,
	 67./90.,         3.,  -91./108.,         0.,         0.,
	23./108., -976./135.,   311./54.,   -19./60.,     17./6.,
	 -1./12., 2383./4100.,        0.,         0., -341./164.,
      4496./1025., -301./82., 2133./4100.,   45./82.,   45./164.,
	 18./41.,    3./205.,         0.,         0.,         0.,
	      0.,    -6./41.,   -3./205.,    -3./41.,     3./41.,
	  6./41.,         0., -1777./4100.,       0.,         0.,
      -341./164., 4496./1025., -289./82., 2193./4100.,   51./82.,
	33./164.,    12./41.,         0.,          1.
   };
   static double  c[11] = {
	41./840.,         0.,         0.,         0.,         0.,
	34./105.,     9./35.,     9./35.,    9./280.,    9./280.,
	41./840.
   };
   static double cp[13] = {
	      0.,         0.,         0.,         0.,         0.,
	34./105.,     9./35.,     9./35.,    9./280.,    9./280.,
	      0.,   41./840.,   41./840.
   };
/* Variables locals */
   int ib, j, k, l, n=nv+np*nn, iret;
   double tt, bet, d, dd, e3, wksp[15*n], *r=wksp, *b=wksp+13*n,
	  *f=wksp+14*n;
/* Bucle en el pas */
   do {
   /* Fem RK7 -> b[] i RK8 -> f[]
    * Promig per coordenades de RK7-RK8 -> d
    * Norma sub-1 de RK7-RK8 -> dd */
      ib=0;
      for (j=0; j<13; j++) {
         memcpy(b, x, n*sizeof(double));
         tt=(*t)+alfa[j]*(*h);
	 for (k=0; k<j; k++, ib++) {
            bet=beta[ib]*(*h);
            for (l=0; l<n; l++)
               b[l]+=bet*r[n*k+l];
         }
	 iret=camp(nv,np,prm,tt,b/*x*/,r+n*j/*f*/);
	 if (iret) return iret;
      }
      dd=d=0;
      for (l=0; l<n; l++) {
         b[l]=f[l]=x[l];
         for (k=0; k<11; k++) {
            bet=(*h)*r[k*n+l];
            b[l]+=bet*c[k];
            f[l]+=bet*cp[k];
         }
         f[l]+=(*h)*(cp[11]*r[11*n+l]+cp[12]*r[12*n+l]);
         d+=fabs(f[l]-b[l]);
         dd+=fabs(f[l]);
      }
      d/=n;	/* No es la norma sub-1, sino el promig dels errors	*/
		/*    comesos a cada coordenada. dd si que es la norma	*/
		/*    sub-1, pero no s'usa per mesurar errors sino per	*/
		/*    relativitzar la tolerancia.			*/
      e3=tol*(1.+.01*dd); /* Es una tolerancia absoluta per valors    */
			  /* petits de les coordenades, relativa      */
			  /* "retardada dos digits" per valors grans. */
   /* Pleguem si O.K. o pas minim. */
      if (d<e3 || fabs(*h)<=hmin) break;
   /* Corregim pas per a reintegrar. */
      (*h)*=.9*pow(e3/d,.125);
      if (fabs((*h))<hmin) (*h)=SIGNE(*h)*hmin;
   /* Torno a fer RK78 */
   } while (1);
/* Guardem temps final. */
   (*t)+=(*h);
/* Guardem punt final. */
   memcpy(x, f, n*sizeof(double));
/* Fem correccio de pas */
   if (d<e3) { double tmp=e3/256; if (d<tmp) d=tmp; }
				/* Si l'error ha estat molt petit	  */
				/*    no volem que el nou pas es dispari. */
   (*h)*=.9*pow(e3/d,.125);  /* Correccio Fehlberg (Stoer 2a ed (7.2.5.16), */
			     /*   noti's que NO es la (7.2.5.17)).          */
   if (fabs(*h)<hmin) {      /* Fem que estigui dins */
      (*h)=SIGNE(*h)*hmin;
      vbprintf(1,ivb,"rk78vp():: t %G : ajusto a pasmin %G !!\n", *t, hmin);
   }
   else if (fabs(*h)>hmax) { /* els limits permesos. */
      (*h)=SIGNE(*h)*hmax;
      vbprintf(2,ivb,"rk78vp():: t %G : ajusto a pasmax %G\n", *t, hmax);
   }
/* Tornem l'estimació de l'error */
   if (eerr!=NULL) *eerr=d;
   return 0;
}
