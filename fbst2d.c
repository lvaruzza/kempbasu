// Copyright (C) 2007 Leonardo Varuzza <varuzza@gmail.com>
//
// This program is free software; you can redistribute it and/or modify it
// under the term of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// _________________

 #include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <math.h>
#include "params.h"

inline double square(double x) {
  return x*x;
}
     
double g (double *k, size_t dim, void *_params) {
  double *params=(double*)_params;

  double ll=params[2]*(log(Na)+log(k[0]))+params[3]*(log(Nb)+log(k[1]))+(-params[2]-params[3])*log(Na*k[0]+Nb*k[1]);
  return ll;
}


double f (double *k, size_t dim, void *params) {
  double likenorm=((double*)params)[0];
  double ll=g(k,dim,params);
  return exp(ll-likenorm);
}

double t (double *k, size_t dim, void *params) {
  double val=g(k,dim,params);
  double cutoff=((double*)params)[1];
  double likenorm=((double*)params)[0];
  //printf("%lf %lf\n",val,cutoff);

  if (val > cutoff) {
    return exp(val-likenorm);
  } else {
    return 0;
  }
}
     
void display_results (char *title, double result, double error) {
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
}
    
const dim = 2;


void integrate(gsl_monte_function *G,gsl_rng *r,double *_res,double *_err) {
  double res, err;
     
  double xl[] = { 0, 0,};
  double xu[] = { 1.0,1.0 };
     
  size_t calls = 500000;

     
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
     
  gsl_monte_vegas_integrate (G, xl, xu, dim, 10000, r, s,
			       &res, &err);
  display_results ("vegas warm-up", res, err);
  
  printf ("converging...\n");
     
  do {
    gsl_monte_vegas_integrate (G, xl, xu, dim, calls/5, r, s,
				 &res, &err);
    printf ("result = % .6f sigma = % .6f "
	    "chisq/dof = %.1f\n", res, err, s->chisq);
  } while (fabs (s->chisq - 1.0) > 0.5);
  
  display_results ("vegas final", res, err);
  
  gsl_monte_vegas_free (s);

  *_res=res;
  *_err=err;
}

 
int main (int argc,char **argv) {
  const gsl_rng_type *RT;
  gsl_rng *r;

  double a = 144;
  double b = 221;

  if (argc >= 2) {
    a=atoi(argv[1]);
    b=atoi(argv[2]);
  }
  printf("%lf %lf\n",a,b);

  double params[]={0,0,a,b};
     
  gsl_monte_function G = { &f, dim, &params };
  gsl_monte_function T = { &t, dim, &params };

  double cutoff=a*log(Na)+b*log(Nb)-(a+b)*log(Na+Nb);

  printf("cutoff=%lf\n",cutoff);
  double phat[2];

  phat[0]=a*1.0/Na;
  phat[1]=b*1.0/Nb;

  double likenorm = g(phat,dim,params);
  params[0] = likenorm;

  printf("likenorm = %lf\n",likenorm);

  gsl_rng_env_setup ();
     
  RT = gsl_rng_default;
  r = gsl_rng_alloc (RT);

  double g,g_err;

  integrate(&G,r,&g,&g_err);

  printf("null = %lf\n",g);

  params[1]=cutoff;
  T.params=params;
  double t,t_err;

  integrate(&T,r,&t,&t_err);

  double err=t/g*sqrt(square(g_err/g)+square(t_err/t));

  printf("t/g = %lf\n",t/g);
  printf("ev = %lf +- %lf\n",1.0-t/g,err);

  gsl_rng_free (r);

  return 0;
}
