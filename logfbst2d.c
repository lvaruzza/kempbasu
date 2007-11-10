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
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <math.h>
     
/* Computation of the integral,
     
   I = int (dx dy dz)/(2pi)^3  1/(1-cos(x)cos(y)cos(z))
     
   over (-pi,-pi,-pi) to (+pi, +pi, +pi).  The exact answer
   is Gamma(1/4)^4/(4 pi^3).  This example is taken from
   C.Itzykson, J.M.Drouffe, "Statistical Field Theory -
   Volume 1", Section 1.1, p21, which cites the original
   paper M.L.Glasser, I.J.Zucker, Proc.Natl.Acad.Sci.USA 74
   1800 (1977) */
     
/* For simplicity we compute the integral over the region 
   (0,0,0) -> (pi,pi,pi) and multiply by 8 */
     
//double exact = 1.3932039296856768591842462603255;

const Na=5929;
const Nb=7460;

const a = 144;
const b = 221;
     
double g (double *k, size_t dim, void *params) {
  return a*(log(Na)+log(k[0]))+b*(log(Nb)+log(k[1]))+(-a-b)*log(Na*k[0]+Nb*k[1]);
}

double t (double *k, size_t dim, void *params) {
  double val=g(k,dim,params);
  double cutoff=((double*)params)[0];
  //printf("%lf %lf\n",val,cutoff);

  if (val > cutoff) {
    return val;
  } else {
    return 0;
  }
}
     
void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
}
    
const dim = 2;


double integrate(gsl_monte_function *G,gsl_rng *r) {
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
  return res;
}

 
int main (void)
{
  const gsl_rng_type *RT;
  gsl_rng *r;

  double params[]={0};
     
  gsl_monte_function G = { &g, dim, 0 };
  gsl_monte_function T = { &t, dim, 0 };

  double cutoff=a*log(Na)+b*log(Nb)-(a+b)*log(Na+Nb);

  printf("cutoff=%lf\n",cutoff);
  double p[]={0.5,0.5};

  printf("g(0.5,0.5) = %lf\n",g(p,2,NULL));
  params[0]=cutoff;
  T.params=params;

  gsl_rng_env_setup ();
     
  RT = gsl_rng_default;
  r = gsl_rng_alloc (RT);

  double g=integrate(&G,r);

  double t=integrate(&T,r);

  printf("t/g = %lf\n",t/g);

  gsl_rng_free (r);

  return 0;
}
