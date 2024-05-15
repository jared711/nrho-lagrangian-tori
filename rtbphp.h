/* - rtbphp.c */
#define RTBPHP_N 6
#define RTBPHP_NV1 42
#define RTBPHP_NV2 168
#define RTBPHP_MU_TL 1.215058560962404e-2
#define RTBPHP_MU_ST 3.040357143e-06
#define RTBPHP_MU_SJ 9.53875e-04
int rtbphp (int nv, int np, void *prm, double t, double x[], double f[]);
double rtbphp_h (int n, int np, void *prm, double x[]);
double rtbpxli (double xmu, int li);

