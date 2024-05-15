/*
 * Copyright (C) 2016 Àlex Haro, Josep-Maria Mondelo
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include "campvp.h"
#include "rtbphp.h"
#include "fluxvp.h"
#include "utils.h"

#define N RTBPHP_N
#define NV1 RTBPHP_NV1
#define NV2 RTBPHP_NV2

#define X x[0]
#define Y x[1]
#define Z x[2]
#define PX x[3]
#define PY x[4]
#define PZ x[5]
#define MU (*((double *)prm))

int rtbphp (int nv, int np, void *prm, double t, double x[], double f[]) {
   double r1x, r2x, r1x2, r2x2, r1ypz, r12, r22, r13, r23, p13, p23,
          p123, p15, p25, p125, p15x, p25x, p125x, dxpx, dypx, dzpx, dypy,
          dzpy, dzpz;
   int j;
   f[0]=PX+Y; f[1]=PY-X; f[2]=PZ;
   r1x=X-MU; r2x=r1x+1; r1x2=r1x*r1x; r2x2=r2x*r2x;
   r1ypz=Y*Y+Z*Z;
   r12=r1x2+r1ypz; r22=r2x2+r1ypz;
   r13=r12*sqrt(r12); r23=r22*sqrt(r22);
   p13=(1-MU)/r13; p23=MU/r23;
   f[3]=PY-(p13*r1x+p23*r2x);
   p123=p13+p23;
   f[4]=-PX-Y*p123;
   f[5]=-Z*p123;
   if (nv>=NV1 || np>0) {
      p15=p13/r12; p25=p23/r22;
      p125=p15+p25;
      p15x=p15*r1x; p25x=p25*r2x;
      p125x=p15x+p25x;
      dxpx=-p123+3*(p15x*r1x+p25x*r2x);
      dzpx=dypx=3*p125x;
      dypx*=Y; dzpx*=Z;
      dzpy=dypy=3*Y*p125;
      dypy=Y*dypy-p123;
      dzpy*=Z;
      dzpz=-p123+3*Z*Z*p125;
   }
   if (nv>=NV1) {
      for (j=0; j<N; j++) {
	 *vr1(N,np,f,0,j)=*vr1(N,np,x,1,j)+*vr1(N,np,x,3,j);
	 *vr1(N,np,f,1,j)=-*vr1(N,np,x,0,j)+*vr1(N,np,x,4,j);
	 *vr1(N,np,f,2,j)=*vr1(N,np,x,5,j);
	 *vr1(N,np,f,3,j)=
	       dxpx*(*vr1(N,np,x,0,j))
			   +dypx*(*vr1(N,np,x,1,j))+dzpx*(*vr1(N,np,x,2,j))
	       +*vr1(N,np,x,4,j);
	 *vr1(N,np,f,4,j)=
	       dypx*(*vr1(N,np,x,0,j))
			   +dypy*(*vr1(N,np,x,1,j))+dzpy*(*vr1(N,np,x,2,j))
	       -*vr1(N,np,x,3,j);
	 *vr1(N,np,f,5,j)=
	       dzpx*(*vr1(N,np,x,0,j))
			   +dzpy*(*vr1(N,np,x,1,j))+dzpz*(*vr1(N,np,x,2,j));
      }
   }
   assert(np==0 && nv<NV2);
   return 0;
}

double rtbphp_h (int n, int np, void *prm, double x[]) {
   double mu=*((double *)prm), xmmu=X-mu, xmmup1=xmmu+1,
	  r12=SQR(xmmu)+SQR(Y)+SQR(Z),
	  r22=SQR(xmmup1)+SQR(Y)+SQR(Z),
	  r1=sqrt(r12), r2=sqrt(r22),
	  p1=(1-mu)/r1, p2=mu/r2;
   return .5*(SQR(PX)+SQR(PY)+SQR(PZ))+Y*PX-X*PY-p1-p2;
}

#undef MU

