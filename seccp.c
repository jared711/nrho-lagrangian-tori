#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "campvp.h"
#include "rk78vp.h"
#include "fluxvp.h"
#include "seccp.h"
#include "vbprintf.h"
#include "utils.h"

/* Globals */
int seccp_maxit=10;
double seccp_factolst=1e4; /* Factor pel qual es multiplica la tolerància
			      per considerar que s'ha sortit de la secció,
			      i considerar no convergència un error */

/* Rutina auxiliar */
#ifdef __GNUC__
inline
#endif
static double seccp_aval (int n, double cp[], double x[]) {
   double ret=0; int i;
   for (i=0; i<n; i++) ret+=cp[i]*x[i];
   ret+=cp[n];
   return ret;
}


/*!
 * \brief Porta una c.i. (vector de variables dependents amb el
 * conveni de campvp.h) fins a una secció de Poincaré. Torna
 * opcionalment la diferencial de l'aplicació de retorn. Permet fer un
 * nombre de seccions arbitrari (paràmetre \c nsec). Pot considerar
 * seccP travessar de negatiu a positiu o de qualsevol manera. Permet
 * bolcar a un fitxer els punts que va obtenint a mida que fa passos
 * d'integració numèrica.
 *
 * El \c camp ha de ser autònom, of course.
 *
 * Usa les variables globals <tt>fluxvp_pasmin,fluxvp_pasmax,fluxvp_tol</tt>,
 * definides a \c fluxvp.c
 *
 * \todo:
 * - acabar de documentar!!
 */
int seccp (int n, int nv, int np, campvp_t camp, void *prm,
      double *t, double x[], double *h, double cp[], int nsec, int isiggrad,
      double tol, int ivb, int idt, double dt[],
      seccp_wrtf_t wrtf, void *wrtf_prm, double maxts) {
   double g, dg, hb, tt, g0, f[n], t0=*t;
   int i, j, it, itfl, vret=0;
/*   assert(!idt || nv>=42);  només valid per 3dof!! */
   itfl=0; i=0;
/* Integro fins al nsec-esim tal (tenint restriccions en compte) */
   while (i<nsec) {
   /* Integro fins sotir de la secció */
      while (fabs(seccp_aval(n,cp,x))<seccp_factolst*tol) {
	 if (wrtf!=NULL) wrtf(n,nv,*t,x,0/*aon*/,wrtf_prm);
	 rk78vp(n,nv,np,prm,camp,t,x,h,fluxvp_pasmin,fluxvp_pasmax,fluxvp_tol,
	       NULL/*eerr*/,ivb);
      }
   /* M'apunto g un cop fora (en part: sé que és diferent de zero) */
      g0=seccp_aval(n,cp,x);
   /* Ara integro fins tallar */
      do {
	 if (wrtf!=NULL) wrtf(n,nv,*t,x,0/*aon*/,wrtf_prm);
	 rk78vp(n,nv,np,prm,camp,t,x,h,fluxvp_pasmin,fluxvp_pasmax,fluxvp_tol,
	       NULL/*eerr*/,ivb);
	 if (fabs(*t-t0)>maxts) {
	    vbprintf(0,ivb,"seccp() : superat maxts %G !!\n", maxts);
	    return -2;
	 }
	 if (++itfl>fluxvp_maxit) {
	    vbprintf(0,ivb,"seccp() : superat fluxvp_maxit %d !!\n",
		  fluxvp_maxit);
	    return -3;
	 }
	 g=seccp_aval(n,cp,x);
	 vbprintf(4,ivb,"seccp(): itfl %d g %G\n", itfl, g);
      } while (g*g0>0);
      if (!isiggrad || g0<0) i++;
   }
/* Guardo el pas abans de fer Newton */
   hb=*h;
/* Ara faig Newton per caure a la secció */
   it=0;
   while (fabs(g=seccp_aval(n,cp,x))>=tol && it<seccp_maxit) {
      camp(n/*nv*/,0/*np*/,prm,0./*t*/,x,f);
      dg=0; for (i=0; i<n; i++) dg+=cp[i]*f[i];
      tt=-g/dg;
      vbprintf(3,ivb,"seccp()::REF: it %d g %G corr %G\n", it, g, tt);
      *h=SIGNE(tt)*fabs(*h);
   /* ATENCIO! Per coherència, hauria de cridar fluxvpw() (dins fluxvp.c) */
      fluxvp(n,nv,np,prm,camp,t,x,h,tt,NULL/*forb*/,ivb);
      it++;
   }
   *h=hb; /* Recupero el pas d'abans de fer Newton */
   vbprintf(3,ivb,"seccp()::ref: it %d g %G\n", it, g);
   if (wrtf!=NULL) wrtf(n,nv,*t,x,1/*aon*/,wrtf_prm);
/* Em queixo si no me n'he sortit */
   if (fabs(g)>=tol) {
      vbprintf(2,ivb,"%s(): no convergència en %d its (g %G)!!\n",
	    __func__, it, g);
      vret=-1;
      if (fabs(g)>=seccp_factolst*tol) {
	 vbprintf(0,ivb,"%s(): error >= %G : no cola!!!\n",
	    __func__, seccp_factolst*tol);
	 return -4;
      }
   }
/* Faig la diferencial del temps de retorn, si me la demanen */
   if (idt) {
      double denom;
   /* Trobo la diferencial de l'aplicació temps de retorn */
      camp(n/*nv*/,0/*np*/,prm,0./*t*/,x,f);
      denom=0; for (i=0; i<n; i++) denom+=cp[i]*f[i];
      for (j=0; j<n; j++) {
	 dt[j]=0;
	 for (i=0; i<n; i++)
	    dt[j]-=cp[i]*(*vr1(n,np,x,i,j));
	 dt[j]/=denom;
      }
   }
/* All done! */
   return vret;
}

