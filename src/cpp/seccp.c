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
int seccp_maxit = 10;
double seccp_factolst = 1e4; /* Factor pel qual es multiplica la tolerància
                 per considerar que s'ha sortit de la secció,
                 i considerar no convergència un error */

/* Rutina auxiliar */
#ifdef __GNUC__
inline
#endif
    static double
    seccp_aval(int n, double cp[], double x[])
{
   double ret = 0;
   int i;
   for (i = 0; i < n; i++)
      ret += cp[i] * x[i];
   ret += cp[n];
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
 * \brief Brings an initial condition (vector of dependent variables with the
 * convention of campvp.h) to a Poincaré section. Optionally returns the
 * differential of the return application. It allows an arbitrary number of
 * sections ( \c nsec parameter). It can consider seccP crossing from negative
 * to positive or in any way. It allows dumping the points it obtains as it
 * takes numerical integration steps to a file.
 *
 * El \c camp ha de ser autònom, of course.
 * The vector field \c camp must be autonomous, of course.
 *
 * Usa les variables globals <tt>fluxvp_pasmin,fluxvp_pasmax,fluxvp_tol</tt>,
 * definides a \c fluxvp.c
 * Uses the global variables <tt>fluxvp_pasmin,fluxvp_pasmax,fluxvp_tol</tt>,
 * defined in \c fluxvp.c
 *
 * \todo:
 * - acabar de documentar!!
 * n is the number of dimensions
 * nv is the number of variables accounting for the variationals
 * np is the number of parameters
 * camp is the vector field
 * prm is a pointer to the parameters (in my case, something pointing to mu) The vector field function deals with it, so it should know what to expect
 * t is the time
 * x is the state vector of variables including the variationals
 * h is the step (gets modified) puts the last step size in there
 * cp is the vector of coefficients (coefficients of the plane)
 * nsec is the number of sections (1,2,3 cuts)
 * isiggrad is a flag indicating whether you want to record the cuts in the direction of the vector cp ()
 * tol is the tolerance
 * ivb is the verbosity level
 * idt is the flag for the differential of the return application
 * dt is the differential of the return application
 * wrtf is the function to write the points to a file
 * wrtf_prm is the name of the file
 * maxts is the maximum time step
 * \return 0 if everything went OK
 * \return -1 if there was no convergence
 * \return -2 if the maximum time step was exceeded
 * \return -3 if the maximum number of iterations was exceeded
 * \return -4 if the error was greater than the tolerance
 * \return -5 if the number of variables is less than 42
 * \return -6 if the number of variables is less than 3
 * \return -7 if the number of variables is less than 2
 * \return -8 if the number of variables is less than 1
 *
 */
int seccp(int n, int nv, int np, campvp_t camp, void *prm,
          double *t, double x[], double *h, double cp[], int nsec, int isiggrad,
          double tol, int ivb, int idt, double dt[],
          seccp_wrtf_t wrtf, void *wrtf_prm, double maxts)
{
   double g, dg, hb, tt, g0, f[n], t0 = *t;
   int i, j, it, itfl, vret = 0;
   /*   assert(!idt || nv>=42);  només valid per 3dof!! */
   itfl = 0;
   i = 0;
   /* Integro fins al nsec-esim tal (tenint restriccions en compte) */
   while (i < nsec)
   {
      /* Integro fins sotir de la secció */
      while (fabs(seccp_aval(n, cp, x)) < seccp_factolst * tol)
      {
         if (wrtf != NULL)
            wrtf(n, nv, *t, x, 0 /*aon*/, wrtf_prm);
         rk78vp(n, nv, np, prm, camp, t, x, h, fluxvp_pasmin, fluxvp_pasmax, fluxvp_tol,
                NULL /*eerr*/, ivb);
      }
      /* M'apunto g un cop fora (en part: sé que és diferent de zero) */
      /*I point to g once out (partially: I know it's non-zero)*/
      g0 = seccp_aval(n, cp, x);
      /* Ara integro fins tallar */
      // Now integrate until cut
      do
      {
         if (wrtf != NULL)
            wrtf(n, nv, *t, x, 0 /*aon*/, wrtf_prm);
         rk78vp(n, nv, np, prm, camp, t, x, h, fluxvp_pasmin, fluxvp_pasmax, fluxvp_tol,
                NULL /*eerr*/, ivb);
         if (fabs(*t - t0) > maxts)
         {
            vbprintf(0, ivb, "seccp() : superat maxts %G !!\n", maxts);
            return -2;
         }
         if (++itfl > fluxvp_maxit)
         {
            vbprintf(0, ivb, "seccp() : superat fluxvp_maxit %d !!\n",
                     fluxvp_maxit);
            return -3;
         }
         g = seccp_aval(n, cp, x);
         vbprintf(4, ivb, "seccp(): itfl %d g %G\n", itfl, g);
      } while (g * g0 > 0);
      if (!isiggrad || g0 < 0)
         i++;
   }
   /* Guardo el pas abans de fer Newton */
   // Save the step before doing Newton
   hb = *h;
   /* Ara faig Newton per caure a la secció */
   // Now do Newton to reach the section
   it = 0;
   while (fabs(g = seccp_aval(n, cp, x)) >= tol && it < seccp_maxit)
   {
      camp(n /*nv*/, 0 /*np*/, prm, 0. /*t*/, x, f); /*fill f with the value of the vector field*/
      dg = 0;
      for (i = 0; i < n; i++)
         dg += cp[i] * f[i];
      tt = -g / dg;
      vbprintf(3, ivb, "seccp()::REF: it %d g %G corr %G\n", it, g, tt);
      *h = SIGNE(tt) * fabs(*h);
      /* ATENCIO! Per coherència, hauria de cridar fluxvpw() (dins fluxvp.c) */
      fluxvp(n, nv, np, prm, camp, t, x, h, tt, NULL /*forb*/, ivb);
      it++;
   }
   *h = hb; /* Recupero el pas d'abans de fer Newton */
   vbprintf(3, ivb, "seccp()::ref: it %d g %G\n", it, g);
   if (wrtf != NULL)
      wrtf(n, nv, *t, x, 1 /*aon*/, wrtf_prm);
   /* Em queixo si no me n'he sortit */
   // Complain if I didn't get out
   if (fabs(g) >= tol)
   {
      vbprintf(2, ivb, "%s(): no convergència en %d its (g %G)!!\n",
               __func__, it, g);
      vret = -1;
      if (fabs(g) >= seccp_factolst * tol)
      {
         vbprintf(0, ivb, "%s(): error >= %G : no cola!!!\n",
                  __func__, seccp_factolst * tol);
         return -4;
      }
   }
   /* Faig la diferencial del temps de retorn, si me la demanen */
   /* Do the differential of the return time, if they request it*/
   if (idt)
   {
      double denom;
      /* Trobo la diferencial de l'aplicació temps de retorn */
      // I find the differential of the return time application
      camp(n /*nv*/, 0 /*np*/, prm, 0. /*t*/, x, f);
      denom = 0;
      for (i = 0; i < n; i++)
         denom += cp[i] * f[i];
      for (j = 0; j < n; j++)
      {
         dt[j] = 0;
         for (i = 0; i < n; i++)
            /* Just not sure why the following line decrements rather than increments (answer the negative sign is in the definition of the time derivative*/
            dt[j] -= cp[i] * (*vr1(n, np, x, i, j)); /*This is taking the contents of the large vector x that has been augmented with the state transition matrix*/
         dt[j] /= denom;
      }
   }
   /* All done! */
   return vret;
}
