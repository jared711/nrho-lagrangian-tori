extern double fluxvp_pasmin, fluxvp_pasmax, fluxvp_tol, fluxvp_pas0,
       fluxvp_pasminfet, fluxvp_pasmaxfet;
extern int fluxvp_maxit;
typedef struct {
   double pasmin, pasmax,
          tol, pas0,
          pasminfet, pasmaxfet;
   int maxit;
} flxprm_t;
int fluxvp (int n, int nv, int np, void *prm, campvp_t camp,
      double *t, double x[], double *h, double T, FILE *forb, int ivb);
