extern int seccp_maxit;
extern double seccp_factolst;

/* aon==0 mentre es vola, ==1 mentre es refina per caure a la secció. */
typedef void (*seccp_wrtf_t)(int n, int nv, double t, double x[], int aon,
      void *prm);
int seccp (int n, int nv, int np, campvp_t camp, void *prm,
      double *t, double x[], double *h, double cp[], int nsec, int isiggrad,
      double tol, int ivb, int idt, double dt[],
      seccp_wrtf_t wrtf, void *wrtf_prm, double maxts);

