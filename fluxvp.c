#include <stdio.h>
#include <math.h>
#include <float.h>
#include "campvp.h"
#include "rk78vp.h"
#include "vbprintf.h"
// #include "utils.h"

/*!
 * \file fluxvp.c
 * \brief Rutina fluxvp().
 */

// /* Globals (coherents amb tol rel 1e-14 a la SOI de Vesta, en km i s) */
// /* ATENCIÓ! Crec que els càlculs del mamotreto estan fet amb les
//  * toleràncies de Vesta!! */
// double fluxvp_pasmin=1e-4, fluxvp_pasmax=1000,
//        fluxvp_tol=4e-10, fluxvp_pas0=.01,
//        fluxvp_pasminfet=DBL_MAX, fluxvp_pasmaxfet=0;
// int fluxvp_maxit=100000;
/* Globals (RTBP) */
double fluxvp_pasmin=1e-6, fluxvp_pasmax=.5,
       fluxvp_tol=4e-14, fluxvp_pas0=.01,
       fluxvp_pasminfet=DBL_MAX, fluxvp_pasmaxfet=0;
int fluxvp_maxit=100000;

/*!
 * \brief Rutina flux temps T. Necessita un camp d'acord amb els
 * convenis de <tt>utils.h</tt> (<tt>typedef campvp_t</tt>).
 */
int fluxvp (int n, int nv, int np, void *prm, campvp_t camp,
      double *t, double x[], double *h, double T, FILE *forb, int ivb) {
   int s, it, iret, i;
   double t0, h2, tf, pmin;
   if (T*(*h)<0) *h=-(*h); T=fabs(T);
   s=(*h>=0) ? 1: -1; it=0; t0=*t;
   if (forb!=NULL) {
      fprintf(forb, "%.16G", *t);
      for (i=0; i<n; i++) fprintf(forb, " %.16G", x[i]);
      fprintf(forb, "\n");
   }
/* Integro fins que, amb el pas recomanat, em passo. */
   while (s*(*t+*h-t0)<T) {
      if (it++ >= fluxvp_maxit)  {
	 vbprintf(0,ivb,
	       "flux(): sobrepassat fluxvp_maxit (%d) integrant!!\n",
	       fluxvp_maxit);
	 return -1;
      }
      iret=rk78vp(n,nv,np,prm,camp,t,x,h,fluxvp_pasmin,fluxvp_pasmax,
	    fluxvp_tol,NULL,ivb);
      if (iret) {
	 vbprintf(2,ivb,
	       "fluxvp(): problemes cridant rk78vp() (vol lliure)!!\n");
	 return -iret;
      }
      if (forb!=NULL) {
	 fprintf(forb, "%.16G", *t);
	 for (i=0; i<n; i++) fprintf(forb, " %.16G", x[i]);
	 fprintf(forb, "\n");
      }
      if (fabs(*h)<fluxvp_pasminfet) fluxvp_pasminfet=fabs(*h);
      if (fabs(*h)>fluxvp_pasmaxfet) fluxvp_pasmaxfet=fabs(*h);
   }
/* Ajusto */
   tf=t0+s*T;
   while (fabs(*t-tf)>fluxvp_tol*(1.+fabs(*t)/100)) {
      if (it++ >= fluxvp_maxit) {
	 vbprintf(0,ivb,"flux(): sobrepassat MAXITFLUX (%d) ajustant!!!\n",
	       fluxvp_maxit);
	 return -2;
      }
      h2=tf-*t;
      pmin=fabs(h2); if (fluxvp_pasmin<pmin) pmin=fluxvp_pasmin;
      iret=rk78vp(n,nv,np,prm,camp,t,x,&h2,pmin/*pasmin*/,fluxvp_pasmax,
	    fluxvp_tol,NULL,ivb);
      if (iret) {
	 vbprintf(2,ivb,
	       "fluxvp(): problemes cridant rk78vp() (ajustament)!!!\n");
	 return iret;
      }
   }
   if (forb!=NULL) {
      fprintf(forb, "%.16G", *t);
      for (i=0; i<n; i++) fprintf(forb, " %.16G", x[i]);
      fprintf(forb, "\n");
   }
   return 0;
}

