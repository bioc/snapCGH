/* Most of the code for the optimization was written by Karl Gegenfurtner
and the original version is available from:
http://archives.math.utk.edu/software/msdos/numerical.analysis/praxis/.html */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <R_ext/Applic.h>
#include <stdlib.h>
#include <time.h>
#include "machine.h"

//#ifdef MSDOS
//static double e[N];	/* to save stack space */
//#endif

/* control parameters */
double tol = SQREPSILON,
       scbd = 1.0,
       step = 1.0;
int    ktm = 1,
       prin = 0,
       maxfun = 0,
       illc = 0;
       
/* some global variables */
static int i, j, k, k2, nl, nf, kl, kt;
static double s, sl, dn, dmin,
       fx, f1, lds, ldt, sf, df1,
       qf1, qd0, qd1, qa, qb, qc,
       m2, m4, small, vsmall, large, 
       vlarge, ldfac, t2;
static double d[N], y[N], z[N],
       q0[N], q1[N], v[N][N];

/* these will be set by praxis to point to it's arguments */
static int n;
double *x;
Exts *ext;
double (*fun)();

/* these will be set by praxis to the global control parameters */
static double h, macheps, t;

float
rand1()	/* return random no between 0 and 1 */
{
  float x;
  //  return (double)(rand()%(8192*2))/(double)(8192*2);
  x = (float) rand()/RAND_MAX;
  return(x);
}

void sort()		/* d and v in descending order */
{
   int k, i, j;
   double s;

   for (i=0; i<n-1; i++) {
       k = i; s = d[i];
       for (j=i+1; j<n; j++) {
           if (d[j] > s) {
	      k = j;
	      s = d[j];
	   }
       }
       if (k > i) {
	  d[k] = d[i];
	  d[i] = s;
	  for (j=0; j<n; j++) {
	      s = v[j][i];
	      v[j][i] = v[j][k];
	      v[j][k] = s;
	  }
       }
   }
}


/* singular value decomposition */

void minfit(int n, double eps, double tol, double ab[N][N], double q[N])
     //int n;
     //double eps, tol, ab[N][N], q[N];
{
   int l, kt, l2, i, j, k;
   double c, f, g, h, s, x, y, z;

   double e[N];		/* plenty of stack on a vax */


   /* householder's reduction to bidiagonal form */
   l = 0;
   x = g = 0.0;
   for (i=0; i<n; i++) {
       e[i] = g; s = 0.0; l = i+1;
       for (j=i; j<n; j++)
	   s += ab[j][i] * ab[j][i];
       if (s < tol) {
	  g = 0.0;
       }
       else {
	  f = ab[i][i];
          if (f < 0.0) 
	     g = sqrt(s);
	  else
	     g = -sqrt(s);
	  h = f*g - s; ab[i][i] = f - g;
	  for (j=l; j<n; j++) {
	      f = 0.0;
	      for (k=i; k<n; k++)
		  f += ab[k][i] * ab[k][j];
	      f /= h;
	      for (k=i; k<n; k++)
		  ab[k][j] += f * ab[k][i];
	  }
       }
       q[i] = g; s = 0.0;
       if (i < n)
	  for (j=l; j<n; j++)
	      s += ab[i][j] * ab[i][j];
       if (s < tol) {
	  g = 0.0;
       }
       else {
	  f = ab[i][i+1];
	  if (f < 0.0)
	     g = sqrt(s);
	  else 
	     g = - sqrt(s);
	  h = f*g - s; ab[i][i+1] = f - g;
	  for (j=l; j<n; j++)
	      e[j] = ab[i][j]/h;
	  for (j=l; j<n; j++) {
	      s = 0;
	      for (k=l; k<n; k++) s += ab[j][k]*ab[i][k];
	      for (k=l; k<n; k++) ab[j][k] += s * e[k];
	  }
       }
       y = fabs(q[i]) + fabs(e[i]);
       if (y > x) x = y;
   }
   /* accumulation of right hand transformations */
   for (i=n-1; i >= 0; i--) {
       if (g != 0.0) {
          h = ab[i][i+1]*g;
	  for (j=l; j<n; j++) ab[j][i] = ab[i][j] / h;
	  for (j=l; j<n; j++) {
              s = 0.0;
	      for (k=l; k<n; k++) s += ab[i][k] * ab[k][j];
	      for (k=l; k<n; k++) ab[k][j] += s * ab[k][i];
	  }
       }
       for (j=l; j<n; j++)
           ab[i][j] = ab[j][i] = 0.0;
       ab[i][i] = 1.0; g = e[i]; l = i;
   }
   /* diagonalization to bidiagonal form */
   eps *= x;
   for (k=n-1; k>= 0; k--) {
       kt = 0;
TestFsplitting:
       if (++kt > 30) {
          e[k] = 0.0;
	  fprintf(stderr, "\n+++ qr failed\n");
       }
       for (l2=k; l2>=0; l2--) {
           l = l2;
	   if (fabs(e[l]) <= eps)
	      goto TestFconvergence;
	   if (fabs(q[l-1]) <= eps)
   	      break;	/* goto Cancellation; */
       }
TestFconvergence:
       z = q[k];
       if (l == k)
          goto Convergence;
       /* shift from bottom 2x2 minor */
       x = q[l]; y = q[k-l]; g = e[k-1]; h = e[k];
       f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2.0*h*y);
       g = sqrt(f*f+1.0);
       if (f <= 0.0)
          f = ((x-z)*(x+z) + h*(y/(f-g)-h))/x;
       else
          f = ((x-z)*(x+z) + h*(y/(f+g)-h))/x;
       /* next qr transformation */
       s = c = 1.0;
       for (i=l+1; i<=k; i++) {
           g = e[i]; y = q[i]; h = s*g; g *= c;
	   if (fabs(f) < fabs(h)) {
	      double fh = f/h;
	      z = fabs(h) * sqrt(1.0 + fh*fh);
	   }
	   else {
	      double hf = h/f;
	      z = (f!=0.0 ? fabs(f)*sqrt(1.0+hf*hf) : 0.0);
	   }
	   e[i-1] = z;
	   if (z == 0.0) 
 	      f = z = 1.0;
	   c = f/z; s = h/z;
	   f = x*c + g*s; g = - x*s + g*c; h = y*s;
	   y *= c;
	   for (j=0; j<n; j++) {
	       x = ab[j][i-1]; z = ab[j][i];
	       ab[j][i-1] = x*c + z*s;
	       ab[j][i] = - x*s + z*c;
	   }
	   if (fabs(f) < fabs(h)) {
	      double fh = f/h;
	      z = fabs(h) * sqrt(1.0 + fh*fh);
	   }
	   else {
	      double hf = h/f;
	      z = (f!=0.0 ? fabs(f)*sqrt(1.0+hf*hf) : 0.0);
	   }
           q[i-1] = z;
	   if (z == 0.0) z = f = 1.0;
	   c = f/z; s = h/z;
	   f = c*g + s*y; x = - s*g + c*y;
       }
       e[l] = 0.0; e[k] = f; q[k] = x;
       goto TestFsplitting;
Convergence:
       if (z < 0.0) {
          q[k] = - z;
	  for (j=0; j<n; j++) ab[j][k] = - ab[j][k];
       }
   }
}


void print()		/* print a line of traces */
{
   printf("\n");
   printf("... chi square reduced to ... %20.10e\n", fx);
   printf("... after %u function calls ...\n", nf);
   printf("... including %u linear searches ...\n", nl);
}

void matprint(s, v, n)
char *s;
double v[N][N];
{
   int k, i;
   
   printf("%s\n", s);
   for (k=0; k<n; k++) {
       for (i=0; i<n; i++) {
           printf("%20.10e ", v[k][i]);
       }
       printf("\n");
   }
}



#ifdef MSDOS
static double tflin[N];
#endif

double
flin(l, j)
double l;
{
   int i;
#ifndef MSDOS
   double tflin[N];
#endif   

   if (j != -1) {		/* linear search */
      for (i=0; i<n; i++)
          tflin[i] = x[i] + l *v[i][j];
   }
   else {			/* search along parabolic space curve */
      qa = l*(l-qd1)/(qd0*(qd0+qd1));
      qb = (l+qd0)*(qd1-l)/(qd0*qd1);
      qc = l*(l+qd0)/(qd1*(qd0+qd1));
      for (i=0; i<n; i++)
          tflin[i] = qa*q0[i]+qb*x[i]+qc*q1[i];
   }
   nf++;
   return (*fun)(n, tflin, ext);
}

void min(j, nits, d2, x1, f1, fk)
double *d2, *x1, f1;
{
   int k, i, dz;
   double x2, xm, f0, f2, fm, d1, t2,
          s, sf1, sx1;
   sf1 = f1; sx1 = *x1;
   k = 0; xm = 0.0; fm = f0 = fx; dz = *d2 < macheps;
   /* find step size */
   s = 0;
   for (i=0; i<n; i++) s += x[i]*x[i];
   s = sqrt(s);
   if (dz)
      t2 = m4*sqrt(fabs(fx)/dmin + s*ldt) + m2*ldt;
   else
      t2 = m4*sqrt(fabs(fx)/(*d2) + s*ldt) + m2*ldt;
   s = s*m4 + t;
   if (dz && t2 > s) t2 = s;
   if (t2 < small) t2 = small;
   if (t2 > 0.01*h) t2 = 0.01 * h;
   if (fk && f1 <= fm) {
      xm = *x1;
      fm = f1;
   }

   if (!fk || fabs(*x1) < t2) {
      *x1 = (*x1 > 0 ? t2 : -t2);
      f1 = flin(*x1, j);
   }

   if (f1 <= fm) {
      xm = *x1;
      fm = f1;
   }
   
next:
   if (dz) {
      x2 = (f0 < f1 ? -(*x1) : 2*(*x1));
      f2 = flin(x2, j);
      if (f2 <= fm) {
         xm = x2;
	 fm = f2;
      }
      *d2 = (x2*(f1-f0) - (*x1)*(f2-f0))/((*x1)*x2*((*x1)-x2));
   }
   d1 = (f1-f0)/(*x1) - *x1**d2; dz = 1;

   if (*d2 <= small) {
      x2 = (d1 < 0 ? h : -h);
   }
   else {
      x2 = - 0.5*d1/(*d2);
   }
   if (fabs(x2) > h)
      x2 = (x2 > 0 ? h : -h);
try:
   f2 = flin(x2, j);
   if ((k < nits) && (f2 > f0)) {
      k++;
      if ((f0 < f1) && (*x1*x2 > 0.0))
         goto next;
      x2 *= 0.5;
      goto try;
   }
   nl++;
   if (f2 > fm) x2 = xm; else fm = f2;
   if (fabs(x2*(x2-*x1)) > small) {
      *d2 = (x2*(f1-f0) - *x1*(fm-f0))/(*x1*x2*(*x1-x2));
   }
   else {
      if (k > 0) *d2 = 0;
   }
   if (*d2 <= small) *d2 = small;
   *x1 = x2; fx = fm;
   if (sf1 < fx) {
      fx = sf1;
      *x1 = sx1;
   }
   if (j != -1)
      for (i=0; i<n; i++)
          x[i] += (*x1)*v[i][j];
}

void quad()	/* look for a minimum along the curve q0, q1, q2	*/
{
   int i;
   double l, s;

   s = fx; fx = qf1; qf1 = s; qd1 = 0.0;
   for (i=0; i<n; i++) {
       s = x[i]; l = q1[i]; x[i] = l; q1[i] = s;
       qd1 = qd1 + (s-l)*(s-l);
   }
   s = 0.0; qd1 = sqrt(qd1); l = qd1;
   if (qd0>0.0 && qd1>0.0 &&nl>=3*n*n) {
      min(-1, 2, &s, &l, qf1, 1);
      qa = l*(l-qd1)/(qd0*(qd0+qd1));
      qb = (l+qd0)*(qd1-l)/(qd0*qd1);
      qc = l*(l+qd0)/(qd1*(qd0+qd1));
   }
   else {
      fx = qf1; qa = qb = 0.0; qc = 1.0;
   }
   qd0 = qd1;
   for (i=0; i<n; i++) {
       s = q0[i]; q0[i] = x[i];
       x[i] = qa*s + qb*x[i] + qc*q1[i];
   }
}

double praxis(optimfn fun1, double *par, int n1, Exts *ext1){

   /* init global extern variables and parameters */
   macheps = EPSILON; h = step; t = tol;
   n = n1; x = par; fun = fun1; ext = ext1;
   small = macheps*macheps; vsmall = small*small;
   large = 1.0/small; vlarge = 1.0/vsmall;
   m2 = sqrt(macheps); m4 = sqrt(m2);
   ldfac = (illc ? 0.1 : 0.01);
   nl = kt = 0; nf = 1; fx = (*fun)(n, x, ext); qf1 = fx;
   t2 = small + fabs(t); t = t2; dmin = small;

   if (h < 100.0*t) h = 100.0*t;
   ldt = h;
   for (i=0; i<n; i++) for (j=0; j<n; j++)
       v[i][j] = (i == j ? 1.0 : 0.0);
   d[0] = 0.0; qd0 = 0.0;
   for (i=0; i<n; i++) q1[i] = x[i];
   if (prin > 1) {
      printf("\n------------- enter function praxis -----------\n");
      printf("... current parameter settings ...\n");
      printf("... scaling ... %20.10e\n", scbd);
      printf("...   tol   ... %20.10e\n", t);
      printf("... maxstep ... %20.10e\n", h);
      printf("...   illc  ... %20u\n", illc);
      printf("...   ktm   ... %20u\n", ktm);
      printf("... maxfun  ... %20u\n", maxfun);
   }
   if (prin) print();
  
mloop:
   sf = d[0];
   s = d[0] = 0.0;

   /* minimize along first direction */
   min(0, 2, &d[0], &s, fx, 0);

   if (s <= 0.0)
      for (i=0; i < n; i++)
          v[i][0] = -v[i][0];
   if ((sf <= (0.9 * d[0])) || ((0.9 * sf) >= d[0]))
      for (i=1; i<n; i++)
          d[i] = 0.0;
   for (k=1; k<n; k++) {
       for (i=0; i<n; i++)
           y[i] = x[i];
       sf = fx;
       illc = illc || (kt > 0);

next:
       kl = k;
       df1 = 0.0;
       if (illc) {        /* random step to get off resolution valley */
          for (i=0; i<n; i++) {
              z[i] = (0.1 * ldt + t2 * pow(10.0,(double)kt)) * (rand1() - 0.5);
              s = z[i];
              for (j=0; j < n; j++)
                  x[j] += s * v[j][i];
  	  }
          fx = (*fun)(n, x, ext);
          nf++;
       }
       /* minimize along non-conjugate directions */ 
       for (k2=k; k2<n; k2++) {  
           sl = fx;
           s = 0.0;
           min(k2, 2, &d[k2], &s, fx, 0);
           if (illc) {
	      double szk = s + z[k2];
              s = d[k2] * szk*szk;
	   }
           else 
	      s = sl - fx;
           if (df1 < s) {
              df1 = s;
              kl = k2;
           }
       }
       if (!illc && (df1 < fabs(100.0 * macheps * fx))) {
          illc = 1;
          goto next;
       }
       /* minimize along conjugate directions */ 
       for (k2=0; k2<=k-1; k2++) {
           s = 0.0;
           min(k2, 2, &d[k2], &s, fx, 0);
       }
       f1 = fx;
       fx = sf;
       lds = 0.0;
       for (i=0; i<n; i++) {
           sl = x[i];
           x[i] = y[i];
           y[i] = sl - y[i];
           sl = y[i];
           lds = lds + sl*sl;
       }
       lds = sqrt(lds);
       if (lds > small) {
          for (i=kl-1; i>=k; i--) {
              for (j=0; j < n; j++)
                  v[j][i+1] = v[j][i];
                  d[i+1] = d[i];
              }
              d[k] = 0.0;
              for (i=0; i < n; i++)
                  v[i][k] = y[i] / lds;
              min(k, 4, &d[k], &lds, f1, 1);
              if (lds <= 0.0) {
                 lds = -lds;
                 for (i=0; i<n; i++)
                     v[i][k] = -v[i][k];
              }
       }
       ldt = ldfac * ldt;
       if (ldt < lds)
          ldt = lds;
       if (prin > 1)
          print();
       t2 = 0.0;
       for (i=0; i<n; i++)
           t2 += x[i]*x[i];
       t2 = m2 * sqrt(t2) + t;
       if (ldt > (0.5 * t2))
          kt = 0;
       else 
	  kt++;
       if (kt > ktm)
          goto fret; 
   }
   /*  try quadratic extrapolation in case    */
   /*  we are stuck in a curved valley        */

   quad();
   dn = 0.0;
   for (i=0; i<n; i++) {
       d[i] = 1.0 / sqrt(d[i]);
       if (dn < d[i])
          dn = d[i];
   }
   if (prin > 2)
      matprint("\n... New Matrix of Directions ...",v,n);
   for (j=0; j<n; j++) {
       s = d[j] / dn;
       for (i=0; i < n; i++)
           v[i][j] *= s;
   }
   if (scbd > 1.0) {       /* scale axis to reduce condition number */
      s = vlarge;
      for (i=0; i<n; i++) {
          sl = 0.0;
          for (j=0; j < n; j++)
              sl += v[i][j]*v[i][j];
          z[i] = sqrt(sl);
          if (z[i] < m4)
             z[i] = m4;
          if (s > z[i])
             s = z[i];
      }
      for (i=0; i<n; i++) {
          sl = s / z[i];
          z[i] = 1.0 / sl;
          if (z[i] > scbd) {
             sl = 1.0 / scbd;
             z[i] = scbd;
          }
      }
   }
   for (i=1; i<n; i++)
       for (j=0; j<=i-1; j++) {
           s = v[i][j];
           v[i][j] = v[j][i];
           v[j][i] = s;
       }
   minfit(n, macheps, vsmall, v, d);
   if (scbd > 1.0) {
      for (i=0; i<n; i++) {
          s = z[i];
          for (j=0; j<n; j++)
              v[i][j] *= s;
      }
      for (i=0; i<n; i++) {
          s = 0.0;
          for (j=0; j<n; j++)
              s += v[j][i]*v[j][i];
          s = sqrt(s);
          d[i] *= s;
          s = 1.0 / s;
          for (j=0; j<n; j++)
              v[j][i] *= s;
      }
   }
   for (i=0; i<n; i++) {
       if ((dn * d[i]) > large)
          d[i] = vsmall;
       else if ((dn * d[i]) < small)
          d[i] = vlarge;
       else 
          d[i] = pow(dn * d[i],-2.0);
   }
   sort();               /* the new eigenvalues and eigenvectors */
   dmin = d[n-1];
   if (dmin < small)
      dmin = small;
   illc = (m2 * d[0]) > dmin;

   if ((maxfun > 0) && (nl > maxfun)) {
      if (prin)
	 printf("\n... maximum number of function calls reached ...\n");
      goto fret;
   }
   goto mloop; 	 /* back to main loop */

fret:
   if (prin > 0) {
         printf("\n... ChiSq reduced to %20.10e ...\n", fx);
	 printf("... after %20u fucntion calls.\n", nf);
   }
   
   return(fx);
}


double fr_one (int n, double *par, void *ex)
{
  Exts *ext;
  double *data;
  ext = ex;

  int nrow1, i;
  nrow1 = ext->nrow;
  data = ext->data;

  double S, mu, Sigma, output;
  double alpha[nrow1], alphahat[nrow1], emis_prob[nrow1], denom[nrow1];

  mu = par[0];
  Sigma = par[1];

  S = exp(Sigma);

  output = 0;

  for(i=0; i < nrow1; i++){
    emis_prob[i] = dnorm(data[i], mu, S, 0);
  }
  
  alpha[0] = emis_prob[0];
  alphahat[0] = alpha[0]/(alpha[0]);

  for(i=1; i < nrow1; i++){
   alpha[i] = alphahat[i-1]*emis_prob[i];
   alphahat[i] = alpha[i]/(alpha[i]);
  }
 
  for(i=0; i < nrow1; i++){
    denom[i] = alpha[i];
    output = output+(log(denom[i]));
  }
   return(-1*(output));
}

double fr_two (int n, double *par, void *ex){

  Exts *ext;
  double *data, *covars;
  ext = ex;
  int nrow1, j, k, t, i;

  nrow1 = ext->nrow;
  data = ext->data;
  covars = ext->covars1;

  double prior, eta, zeta, omega, output, pr1, p1, p2, rate1;
  double S[2], mu[2], Sigma[2];
  double gammaA[2][2], gammaB[2][2], gammaC[2][2], alpha[2][nrow1], alphahat[2][nrow1], emis_prob[2][nrow1];
  double denom[nrow1];

  //initialize output
    output = 0;

  mu[0] = par[0];
  mu[1] = par[1];
  Sigma[0] = par[2];
  Sigma[1] = par[3];
  prior = par[4];
  eta = par[5];
  zeta = par[6];
  omega = par[7];
  S[0] = exp(Sigma[0]);
  //Change this later to take the VAR.FIXED argument
  S[1] = exp(Sigma[1]);
  
  if (prior < 150) {
    pr1 = exp(prior)/(1+exp(prior)); 
      }
  else {
    pr1 = 1;
      }
  
  if (eta < 150) {
    p1 = exp(eta)/(1+exp(eta));
      }
  else {
    p1 = 1;
      }

  if (zeta < 150) {
    p2 = exp(zeta)/(1+exp(zeta)); 
      }
  else {
    p2 = 1;
      }

  if (omega < 50) {
    rate1 = exp(omega); 
      }
  else {
    rate1 = exp(50);
      }

  gammaA[0][0] = (1-p1); gammaA[0][1] = p1; gammaA[1][0] = p2; gammaA[1][1] = (1-p2);
  //  gammaA[0][0] = 1; gammaA[0][1] = 2; gammaA[1][0] = 3; gammaA[1][1] = 4;
  gammaB[0][0] = p1; gammaB[0][1] = -p1; gammaB[1][0] = -p2; gammaB[1][1] = p2;
  
   //Rprintf("%f \t %f\n%f \t %f\n", gammaA[0][0],gammaA[0][1],gammaA[1][0],gammaA[1][1]);

  for(j = 0; j < nrow1; j++){
    for(k = 0; k < 2; k++){
      if(S[k] > 0.001){
	emis_prob[k][j] = dnorm(data[j], mu[k], S[k], 0);
      }
      else {
	if(data[j] >= mu[k]*(0.999) && data[j] <= mu[k]*(1.001)){
	  emis_prob[k][j] = 1;
	}
	else {
	  emis_prob[k][j] = 0;
	}
      }
    }
  }

  alpha[0][0] = pr1*(emis_prob[0][0]);
  alpha[1][0] = (1 - pr1)*emis_prob[1][0];

  alphahat[0][0] = alpha[0][0]/(alpha[0][0]+alpha[1][0]);
  alphahat[1][0] = alpha[1][0]/(alpha[0][0]+alpha[1][0]);

  

  for(t = 1; t < nrow1; t++){
    gammaC[0][0] = gammaA[0][0] + ((exp(-(pow(covars[t-1],rate1))))*gammaB[0][0]);
    gammaC[0][1] = gammaA[0][1] + ((exp(-(pow(covars[t-1],rate1))))*gammaB[0][1]);
    gammaC[1][0] = gammaA[1][0] + ((exp(-(pow(covars[t-1],rate1))))*gammaB[1][0]);
    gammaC[1][1] = gammaA[1][1] + ((exp(-(pow(covars[t-1],rate1))))*gammaB[1][1]);

    //    Rprintf("%f \t %f\n%f \t %f\n", gammaC[0][0],gammaC[0][1],gammaC[1][0],gammaC[1][1]);

    for(i = 0; i < 2; i++){
      alpha[i][t] = (alphahat[0][t-1]*gammaC[0][i] + alphahat[1][t-1]*gammaC[1][i])*emis_prob[i][t];
      //          Rprintf("alpha %d %d is %f\n", i, t, alpha[i][t]);
	}
    for(i = 0; i < 2; i++){
      alphahat[i][t] = alpha[i][t]/(alpha[0][t]+alpha[1][t]);
    }
  }

  for(i=0; i < nrow1; i++){
    denom[i] = alpha[0][i]+alpha[1][i];
    //    Rprintf("%f\n", denom[i]);
    output = output+(log(denom[i]));
  }

  return(-1*(output));
}

double fr_three(int n, double *par, void *ex){

  Exts *ext;
  double *data, *covars;
  ext = ex;

  int nrow1, j, k, t, i, m;
  nrow1 = ext->nrow;
  data = ext->data;
  covars = ext->covars1;

  double prior1, prior2, eta, zeta, theta, beta, gamma, xi, omega, output, pr1, pr2, p1, p2, p3, p4, p5, p6, rate1;
  double S[3], mu[3], Sigma[3];
  double gammaA[3][3], gammaB[3][3], gammaC[3][3], alpha[3][nrow1], alphahat[3][nrow1], emis_prob[3][nrow1];
  double denom[nrow1];

  //initialize output
    output = 0;
    
  mu[0] = par[0];
  mu[1] = par[1];
  mu[2] = par[2];
  Sigma[0] = par[3];
  Sigma[1] = par[4];
  Sigma[2] = par[5];
  prior1 = par[6];
  prior2 = par[7];
  eta = par[8];
  zeta = par[9];
  theta = par[10];
  beta = par[11];
  gamma = par[12];
  xi = par[13];
  omega = par[14];
  S[0] = exp(Sigma[0]);
  //Change this later to take the VAR.FIXED argument
  S[1] = exp(Sigma[1]);
  S[2] = exp(Sigma[2]);
  
  if (prior1 < 150) {
    pr1 = exp(prior1)/(1+exp(prior1)); 
      }
  else {
    pr1 = 1;
      }

   if (prior2 < 150) {
    pr2 = (1 - pr1)*exp(prior2)/(1+exp(prior2)); 
      }
  else {
    pr2 = (1- pr1);
      }
  
  if (eta < 150) {
    p1 = exp(eta)/(1+exp(eta));
      }
  else {
    p1 = 1;
      }

  if (zeta < 150) {
    p2 = (1 - p1)*exp(zeta)/(1+exp(zeta)); 
      }
  else {
    p2 = (1 - p1);
      }

  if (theta < 150) {
    p3 = exp(theta)/(1+exp(theta));
      }
  else {
    p3 = 1;
      }

  if (beta < 150) {
    p4 = (1 - p3)*exp(beta)/(1+exp(beta)); 
      }
  else {
    p4 = (1 - p3);
      }

  if (gamma < 150) {
    p5 = exp(gamma)/(1+exp(gamma));
      }
  else {
    p5 = 1;
      }

  if (xi < 150) {
    p6 = (1 - p5)*exp(xi)/(1+exp(xi)); 
      }
  else {
    p6 = (1 - p5);
      }

  if (omega < 50) {
    rate1 = exp(omega); 
      }
  else {
    rate1 = exp(50);
      }

  gammaA[0][0] = (1-p1-p2); gammaA[0][1] = p1; gammaA[0][2] = p2; gammaA[1][0] = p3; gammaA[1][1] = (1-p3-p4); gammaA[1][2] = p4; gammaA[2][0] = p5; gammaA[2][1] = p6; gammaA[2][2] = (1-p5-p6);
 
  gammaB[0][0] = (p1+p2); gammaB[0][1] = -p1; gammaB[0][2] = -p2; gammaB[1][0] = -p3; gammaB[1][1] = (p3+p4); gammaB[1][2] = -p4; gammaB[2][0] = -p5; gammaB[2][1] = -p6; gammaB[2][2] = (p5+p6);
  
  //  Rprintf("%f \t %f\n%f \t %f\n", gammaA[0][0],gammaA[0][1],gammaA[1][0],gammaA[1][1]);

  for(j = 0; j < nrow1; j++){
    for(k = 0; k < 3; k++){
      if(S[k] > 0.001){
	emis_prob[k][j] = dnorm(data[j], mu[k], S[k], 0);
      }
      else {
	if(data[j] >= mu[k]*(0.999) && data[j] <= mu[k]*(1.001)){
	  emis_prob[k][j] = 1;
	}
	else {
	  emis_prob[k][j] = 0;
	}
      }
    }
  }

  alpha[0][0] = pr1*(emis_prob[0][0]);
  alpha[1][0] = pr2*(emis_prob[1][0]);
  alpha[2][0] = (1-pr1-pr2)*(emis_prob[2][0]);

  alphahat[0][0] = alpha[0][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]);
  alphahat[1][0] = alpha[1][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]);
  alphahat[2][0] = alpha[2][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]);

  for(t = 1; t < nrow1; t++){

    for(j = 0; j < 3; j++){
      for(m = 0; m < 3; m++){
	gammaC[j][m] = gammaA[j][m] + ((exp(-(pow(covars[t-1],rate1))))*gammaB[j][m]);
      }
    }
   
    //    Rprintf("GammaC:\n%f \t %f \t %f\n%f \t %f \t %f\n%f \t %f \t %f\n", gammaC[0][0],gammaC[0][1],gammaC[0][2],gammaC[1][0],gammaC[1][1],gammaC[1][2],gammaC[2][0],gammaC[2][1],gammaC[2][2]);

    for(i = 0; i < 3; i++){
      alpha[i][t] = (alphahat[0][t-1]*gammaC[0][i] + alphahat[1][t-1]*gammaC[1][i] + alphahat[2][t-1]*gammaC[2][i])*emis_prob[i][t];
      //  Rprintf("alpha %d %d is %f\n", i, t, alpha[i][t]);
	}
    for(i = 0; i < 3; i++){
      alphahat[i][t] = alpha[i][t]/(alpha[0][t]+alpha[1][t]+alpha[2][t]);
    }
  }

  

  for(i=0; i < nrow1; i++){
    denom[i] = alpha[0][i]+alpha[1][i]+alpha[2][i];
    //    Rprintf("%f\n", denom[i]);
    output = output+(log(denom[i]));
  }

  return (-1*(output));
}

double fr_four(int n, double *par, void *ex){

  Exts *ext;
  double *data, *covars;
  ext = ex;

  int nrow1, j, k, t, i, m;
  nrow1 = ext->nrow;
  data = ext->data;
  covars = ext->covars1;

  double prior1, prior2, prior3, eta, zeta, nu, theta, beta, phi, gamma, delt, epsilon, lambda, rho, xi, omega, output, pr1, pr2, pr3, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, rate1;
  double S[4], mu[4], Sigma[4];
  double gammaA[4][4], gammaB[4][4], gammaC[4][4], alpha[4][nrow1], alphahat[4][nrow1], emis_prob[4][nrow1];
  double denom[nrow1];

  //initialize output
    output = 0;
    
  mu[0] = par[0];
  mu[1] = par[1];
  mu[2] = par[2];
  mu[3] = par[3];
  Sigma[0] = par[4];
  Sigma[1] = par[5];
  Sigma[2] = par[6];
  Sigma[3] = par[7];
  prior1 = par[8];
  prior2 = par[9];
  prior3 = par[10];
  eta = par[11];
  zeta = par[12];
  nu = par[13];
  theta = par[14];
  beta = par[15];
  phi = par[16];
  gamma = par[17];
  delt = par[18];
  epsilon = par[19];
  lambda = par[20];
  rho = par[21];
  xi = par[22];
  omega = par[23];
  S[0] = exp(Sigma[0]);
  //Change this later to take the VAR.FIXED argument
  S[1] = exp(Sigma[1]);
  S[2] = exp(Sigma[2]);
  S[3] = exp(Sigma[3]);
  
  if (prior1 < 150) {
    pr1 = exp(prior1)/(1+exp(prior1)); 
      }
  else {
    pr1 = 1;
      }

   if (prior2 < 150) {
    pr2 = (1 - pr1)*exp(prior2)/(1+exp(prior2)); 
      }
  else {
    pr2 = (1- pr1);
      }

   if (prior3 < 150) {
    pr3 = (1 -pr1-pr2)*exp(prior3)/(1+exp(prior3)); 
      }
  else {
    pr3 = (1-pr1-pr2);
      }
  
  if (eta < 150) {
    p1 = exp(eta)/(1+exp(eta));
      }
  else {
    p1 = 1;
      }

  if (zeta < 150) {
    p2 = (1 - p1)*exp(zeta)/(1+exp(zeta)); 
      }
  else {
    p2 = (1 - p1);
      }

  if (nu < 150) {
    p3 = (1 -p1-p2)*exp(nu)/(1+exp(nu)); 
      }
  else {
    p3 = (1 - p1 - p2);
      }

  if (theta < 150) {
    p4 = exp(theta)/(1+exp(theta));
      }
  else {
    p4 = 1;
      }

  if (beta < 150) {
    p5 = (1 - p4)*exp(beta)/(1+exp(beta)); 
      }
  else {
    p5 = (1 - p4);
      }

  if (phi < 150) {
    p6 = (1 - p4 - p5)*exp(phi)/(1+exp(phi)); 
      }
  else {
    p6 = (1 - p4 - p5);
      }

  if (gamma < 150) {
    p7 = exp(gamma)/(1+exp(gamma));
      }
  else {
    p7 = 1;
      }

  if (delt < 150) {
    p8 = (1 - p7)*exp(delt)/(1+exp(delt)); 
      }
  else {
    p8 = (1 - p7);
      }

  if (epsilon < 150) {
    p9 = (1 - p7 - p8)*exp(epsilon)/(1+exp(epsilon)); 
      }
  else {
    p9 = (1 - p7 - p8);
      }

  if (lambda < 150) {
    p10 = exp(lambda)/(1+exp(lambda));
      }
  else {
    p10 = 1;
      }

  if (rho < 150) {
    p11 = (1 - p10)*exp(rho)/(1+exp(rho)); 
      }
  else {
    p11 = (1 - p10);
      }

  if (xi < 150) {
    p12 = (1 - p10 - p11)*exp(xi)/(1+exp(xi)); 
      }
  else {
    p12 = (1 - p10 - p11);
      }

  if (omega < 50) {
    rate1 = exp(omega); 
      }
  else {
    rate1 = exp(50);
      }

  gammaA[0][0] = (1-p1-p2-p3); gammaA[0][1] = p1; gammaA[0][2] = p2; gammaA[0][3] = p3; gammaA[1][0] = p4; gammaA[1][1] = (1-p4-p5-p6); gammaA[1][2] = p5; gammaA[1][3] = p6; gammaA[2][0] = p7; gammaA[2][1] = p8; gammaA[2][2] = (1-p7-p8-p9); gammaA[2][3] = p9; gammaA[3][0] = p10; gammaA[3][1] = p11; gammaA[3][2] = p12; gammaA[3][3] = (1-p10-p11-p12);
 
  gammaB[0][0] = (p1+p2+p3); gammaB[0][1] = -p1; gammaB[0][2] = -p2; gammaB[0][3] = -p3; gammaB[1][0] = -p4; gammaB[1][1] = (p4+p5+p6); gammaB[1][2] = -p5; gammaB[1][3] = -p6; gammaB[2][0] = -p7; gammaB[2][1] = -p8; gammaB[2][2] = (p7+p8+p9); gammaB[2][3] = -p9; gammaB[3][0] = -p10; gammaB[3][1] = -p11; gammaB[3][2] = -p12; gammaB[3][3] = (p10+p11+p12);
  


  for(j = 0; j < nrow1; j++){
    for(k = 0; k < 4; k++){
      if(S[k] > 0.001){
	emis_prob[k][j] = dnorm(data[j], mu[k], S[k], 0);
      }
      else {
	if(data[j] >= mu[k]*(0.999) && data[j] <= mu[k]*(1.001)){
	  emis_prob[k][j] = 1;
	}
	else {
	  emis_prob[k][j] = 0;
	}
      }
    }
  }

  alpha[0][0] = pr1*(emis_prob[0][0]);
  alpha[1][0] = pr2*(emis_prob[1][0]);
  alpha[2][0] = pr2*(emis_prob[2][0]);
  alpha[3][0] = (1-pr1-pr2-pr3)*(emis_prob[3][0]);

  alphahat[0][0] = alpha[0][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]+alpha[3][0]);
  alphahat[1][0] = alpha[1][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]+alpha[3][0]);
  alphahat[2][0] = alpha[2][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]+alpha[3][0]);
  alphahat[3][0] = alpha[3][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]+alpha[3][0]);

  for(t = 1; t < nrow1; t++){
    for(j = 0; j < 4; j++){
      for(m = 0; m < 4; m++){
	gammaC[j][m] = gammaA[j][m] + ((exp(-(pow(covars[t-1],rate1))))*gammaB[j][m]);
      }
    }

    for(i = 0; i < 4; i++){
      alpha[i][t] = (alphahat[0][t-1]*gammaC[0][i] + alphahat[1][t-1]*gammaC[1][i] + alphahat[2][t-1]*gammaC[2][i] + alphahat[3][t-1]*gammaC[3][i])*emis_prob[i][t];
	}
    for(i = 0; i < 4; i++){
      alphahat[i][t] = alpha[i][t]/(alpha[0][t]+alpha[1][t]+alpha[2][t]+alpha[3][t]);
    }
  }

  for(i=0; i < nrow1; i++){
    denom[i] = alpha[0][i]+alpha[1][i]+alpha[2][i]+alpha[3][i];
    //    Rprintf("%f\n", denom[i]);
    output = output+(log(denom[i]));
  }

  return (-1*(output));
}

double fr_five(int n, double *par, void *ex){

  Exts *ext;
  double *data, *covars;
  ext = ex;

  int nrow1, j, k, t, i, m;
  nrow1 = ext->nrow;
  data = ext->data;
  covars = ext->covars1;

  double prior1, prior2, prior3, prior4, eta, zeta, nu, omikron, theta, beta, phi, kappa, gamma, delt, epsilon, tau, lambda, rho, xi, iota, chi, upsilon, psi, aux, omega, output, pr1, pr2, pr3, pr4, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, rate1;
  double S[5], mu[5], Sigma[5];
  double gammaA[5][5], gammaB[5][5], gammaC[5][5], alpha[5][nrow1], alphahat[5][nrow1], emis_prob[5][nrow1];
  double denom[nrow1];

  //initialize output
    output = 0;
    
  mu[0] = par[0];
  mu[1] = par[1];
  mu[2] = par[2];
  mu[3] = par[3];
  mu[4] = par[4];
  Sigma[0] = par[5];
  Sigma[1] = par[6];
  Sigma[2] = par[7];
  Sigma[3] = par[8];
  Sigma[4] = par[9];
  prior1 = par[10];
  prior2 = par[11];
  prior3 = par[12];
  prior4 = par[13];
  eta = par[14];
  zeta = par[15];
  nu = par[16];
  omikron = par[17];
  theta = par[18];
  beta = par[19];
  phi = par[20];
  kappa = par[21];
  gamma = par[22];
  delt = par[23];
  epsilon = par[24];
  tau = par[25];
  lambda = par[26];
  rho = par[27];
  xi = par[28];
  iota = par[29];
  chi = par[30];
  upsilon = par[31];
  psi = par[32];
  aux = par[33];
  omega = par[34];
  S[0] = exp(Sigma[0]);
  //Change this later to take the VAR.FIXED argument
  S[1] = exp(Sigma[1]);
  S[2] = exp(Sigma[2]);
  S[3] = exp(Sigma[3]);
  S[4] = exp(Sigma[4]);
  
  if (prior1 < 150) {
    pr1 = exp(prior1)/(1+exp(prior1)); 
      }
  else {
    pr1 = 1;
      }

   if (prior2 < 150) {
    pr2 = (1 - pr1)*exp(prior2)/(1+exp(prior2)); 
      }
  else {
    pr2 = (1- pr1);
      }

   if (prior3 < 150) {
    pr3 = (1 -pr1-pr2)*exp(prior3)/(1+exp(prior3)); 
      }
  else {
    pr3 = (1-pr1-pr2);
      }

   if (prior4 < 150) {
    pr4 = (1 -pr1-pr2-pr3)*exp(prior4)/(1+exp(prior4)); 
      }
  else {
    pr4 = (1-pr1-pr2-pr3);
      }
  
  if (eta < 150) {
    p1 = exp(eta)/(1+exp(eta));
      }
  else {
    p1 = 1;
      }

  if (zeta < 150) {
    p2 = (1 - p1)*exp(zeta)/(1+exp(zeta)); 
      }
  else {
    p2 = (1 - p1);
      }

  if (nu < 150) {
    p3 = (1 -p1-p2)*exp(nu)/(1+exp(nu)); 
      }
  else {
    p3 = (1 - p1 - p2);
      }

  if (omikron < 150) {
    p4 = (1 -p1-p2-p3)*exp(omikron)/(1+exp(omikron)); 
      }
  else {
    p4 = (1 - p1 - p2 - p3);
      }

  if (theta < 150) {
    p5 = exp(theta)/(1+exp(theta));
      }
  else {
    p5 = 1;
      }

  if (beta < 150) {
    p6 = (1 - p5)*exp(beta)/(1+exp(beta)); 
      }
  else {
    p6 = (1 - p5);
      }

  if (phi < 150) {
    p7 = (1 - p5 - p6)*exp(phi)/(1+exp(phi)); 
      }
  else {
    p7 = (1 - p5 - p6);
      }

  if (kappa < 150) {
    p8 = (1 - p5 - p6 - p7)*exp(kappa)/(1+exp(kappa)); 
      }
  else {
    p8 = (1 - p5 - p6 - p7);
      }

  if (gamma < 150) {
    p9 = exp(gamma)/(1+exp(gamma));
      }
  else {
    p9 = 1;
      }

  if (delt < 150) {
    p10 = (1 - p9)*exp(delt)/(1+exp(delt)); 
      }
  else {
    p10 = (1 - p9);
      }

  if (epsilon < 150) {
    p11 = (1 - p9 - p10)*exp(epsilon)/(1+exp(epsilon)); 
      }
  else {
    p11 = (1 - p9 - p10);
      }

  if (tau < 150) {
    p12 = (1 - p9 - p10 - p11)*exp(tau)/(1+exp(tau)); 
      }
  else {
    p12 = (1 - p9 - p10 - p11);
      }

  if (lambda < 150) {
    p13 = exp(lambda)/(1+exp(lambda));
      }
  else {
    p13 = 1;
      }

  if (rho < 150) {
    p14 = (1 - p13)*exp(rho)/(1+exp(rho)); 
      }
  else {
    p14 = (1 - p13);
      }

  if (xi < 150) {
    p15 = (1 - p13 - p14)*exp(xi)/(1+exp(xi)); 
      }
  else {
    p15 = (1 - p13 - p14);
      }

  if (iota < 150) {
    p16 = (1 - p13 - p14 - p15)*exp(iota)/(1+exp(iota)); 
      }
  else {
    p16 = (1 - p13 - p14 - p15);
      }

  if (chi < 150) {
    p17 = exp(chi)/(1+exp(chi));
      }
  else {
    p17 = 1;
      }

  if (upsilon < 150) {
    p18 = (1 - p17)*exp(upsilon)/(1+exp(upsilon)); 
      }
  else {
    p18 = (1 - p17);
      }

  if (psi < 150) {
    p19 = (1 - p17 - p18)*exp(psi)/(1+exp(psi)); 
      }
  else {
    p19 = (1 - p17 - p18);
      }

  if (aux < 150) {
    p20 = (1 - p17 - p18 - p19)*exp(aux)/(1+exp(aux)); 
      }
  else {
    p20 = (1 - p17 - p18 - p19);
      }

  if (omega < 50) {
    rate1 = exp(omega); 
      }
  else {
    rate1 = exp(50);
      }

  gammaA[0][0] = (1-p1-p2-p3-p4); gammaA[0][1] = p1; gammaA[0][2] = p2; gammaA[0][3] = p3; gammaA[0][4] = p4;  gammaA[1][0] = p5; gammaA[1][1] = (1-p5-p6-p7-p8); gammaA[1][2] = p6; gammaA[1][3] = p7; gammaA[1][4] = p8; gammaA[2][0] = p9; gammaA[2][1] = p10; gammaA[2][2] = (1-p9-p10-p11-p12); gammaA[2][3] = p11; gammaA[2][4] = p12; gammaA[3][0] = p13; gammaA[3][1] = p14; gammaA[3][2] = p15; gammaA[3][3] = (1-p13-p14-p15-p16); gammaA[3][4] = p16; gammaA[4][0] = p17; gammaA[4][1] = p18; gammaA[4][2] = p19; gammaA[4][3] = p20; gammaA[4][4] = (1-p17-p18-p19-p20);
 
  gammaB[0][0] = (p1+p2+p3+p4); gammaB[0][1] = -p1; gammaB[0][2] = -p2; gammaB[0][3] = -p3; gammaB[0][4] = -p4;  gammaB[1][0] = -p5; gammaB[1][1] = (p5+p6+p7+p8); gammaB[1][2] = -p6; gammaB[1][3] = -p7; gammaB[1][4] = -p8; gammaB[2][0] = -p9; gammaB[2][1] = -p10; gammaB[2][2] = (p9+p10+p11+p12); gammaB[2][3] = -p11; gammaB[2][4] = -p12; gammaB[3][0] = -p13; gammaB[3][1] = -p14; gammaB[3][2] = -p15; gammaB[3][3] = (p13+p14+p15+p16); gammaB[3][4] = -p16; gammaB[4][0] = -p17; gammaB[4][1] = -p18; gammaB[4][2] = -p19; gammaB[4][3] = -p20; gammaB[4][4] = (p17+p18+p19+p20);
  


  for(j = 0; j < nrow1; j++){
    for(k = 0; k < 5; k++){
      if(S[k] > 0.001){
	emis_prob[k][j] = dnorm(data[j], mu[k], S[k], 0);
      }
      else {
	if(data[j] >= mu[k]*(0.999) && data[j] <= mu[k]*(1.001)){
	  emis_prob[k][j] = 1;
	}
	else {
	  emis_prob[k][j] = 0;
	}
      }
    }
  }

  alpha[0][0] = pr1*(emis_prob[0][0]);
  alpha[1][0] = pr2*(emis_prob[1][0]);
  alpha[2][0] = pr2*(emis_prob[2][0]);
  alpha[3][0] = pr4*(emis_prob[3][0]);
  alpha[4][0] = (1-pr1-pr2-pr3-pr4)*(emis_prob[4][0]);

  alphahat[0][0] = alpha[0][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]+alpha[3][0]+alpha[4][0]);
  alphahat[1][0] = alpha[1][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]+alpha[3][0]+alpha[4][0]);
  alphahat[2][0] = alpha[2][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]+alpha[3][0]+alpha[4][0]);
  alphahat[3][0] = alpha[3][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]+alpha[3][0]+alpha[4][0]);
  alphahat[4][0] = alpha[4][0]/(alpha[0][0]+alpha[1][0]+alpha[2][0]+alpha[3][0]+alpha[4][0]);

  for(t = 1; t < nrow1; t++){
     for(j = 0; j < 5; j++){
      for(m = 0; m < 5; m++){
	gammaC[j][m] = gammaA[j][m] + ((exp(-(pow(covars[t-1],rate1))))*gammaB[j][m]);
      }
    }


    for(i = 0; i < 5; i++){
      alpha[i][t] = (alphahat[0][t-1]*gammaC[0][i] + alphahat[1][t-1]*gammaC[1][i] + alphahat[2][t-1]*gammaC[2][i] + alphahat[3][t-1]*gammaC[3][i] + alphahat[4][t-1]*gammaC[4][i])*emis_prob[i][t];
	}
    for(i = 0; i < 5; i++){
      alphahat[i][t] = alpha[i][t]/(alpha[0][t]+alpha[1][t]+alpha[2][t]+alpha[3][t]+alpha[4][t]);
    }
  }

  for(i=0; i < nrow1; i++){
    denom[i] = alpha[0][i]+alpha[1][i]+alpha[2][i]+alpha[3][i]+alpha[4][i];
    output = output+(log(denom[i]));
  }

  return (-1*(output));
}


void one_state_praxis(int *nrow, double *xin, double *data, double *result){
  Exts *ext;
  ext = Calloc(1, Exts);

  ext->data = data;
  ext->nrow = *nrow;

  srand( (unsigned)time( NULL ) );
  *result = praxis(fr_one, xin, 2, ext);
}

void two_states_praxis(int *nrow, double *xin, double *data, double *covars, double *result){
  
  Exts *ext;
  ext = Calloc(1, Exts);

  ext->data = data;
  ext->nrow = *nrow;
  ext->covars1 = covars;

  srand( (unsigned)time( NULL ) );
  *result = praxis(fr_two, xin, 8, ext);
}

void three_states_praxis(int *nrow, double *xin, double *data, double *covars, double *result){

  Exts *ext;
  ext = Calloc(1, Exts);

  ext->data = data;
  ext->nrow = *nrow;
  ext->covars1 = covars;

  srand( (unsigned)time( NULL ) );
  *result = praxis(fr_three, xin, 15, ext);
}

void four_states_praxis(int *nrow, double *xin, double *data, double *covars, double *result){

  Exts *ext;
  ext = Calloc(1, Exts);

  ext->data = data;
  ext->nrow = *nrow;
  ext->covars1 = covars;

  srand( (unsigned)time( NULL ) );
  *result = praxis(fr_four, xin, 24, ext);
}

void five_states_praxis(int *nrow, double *xin, double *data, double *covars, double *result){

  Exts *ext;
  ext = Calloc(1, Exts);

  ext->data = data;
  ext->nrow = *nrow;
  ext->covars1 = covars;

  srand( (unsigned)time( NULL ) );
  *result = praxis(fr_five, xin, 35, ext);
}
