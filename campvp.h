typedef int (*campvp_t) (int n, int np, void *prm,
   double t, double x[], double f[]);
typedef double (*hamvp_t) (int n, int np, void *prm, double x[]);

double *vrp(int n, double *x, int i);
double *vr1(int n, int np, double *x, int i, int j);
double *vr2(int n, int np, double *x, int i, int j, int k); /* Ull: j>=i !! */
double *vr1u(int n, int np, double *x, int i);
double *vr1v(int n, int np, double *x, int i);
double *vr2uv(int n, int np, double *x, int i);
