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

#include "kempbasu.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <math.h>

MCConfig _MC_DEFAULTS={ 
  100*1000,   // warmup
  1000*1000  // calls
};

MCConfig *MC_DEFAULTS=&_MC_DEFAULTS;


void display_results (char *title, double result, double error) {
  fprintf(stderr,"%s ==================\n", title);
  fprintf(stderr,"result = %.4g\n", result);
  fprintf(stderr,"sigma  = %.4g\n", error);
}
    

void mc_integrate(gsl_monte_function *G,
		  gsl_rng *r,
		  double *_res,double *_err,
		  MCConfig *config) {

  size_t i;

  double *xl = (double*)malloc(G->dim*sizeof(double));
  double *xu = (double*)malloc(G->dim*sizeof(double));


  for (i=0;i<G->dim;i++) {
    xl[i]=0.00001;
    xu[i]=1.0;
  }


  free(xl);
  free(xu);
}

void mc_integrate1(gsl_monte_function *G,
		   gsl_rng *r,
		   double *xl,
		   double *xu,
		   double *_res,double *_err,
		   MCConfig *config) {

  double res, err;
  size_t i;

  size_t calls = config->calls;

     
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (G->dim);
     
  gsl_monte_vegas_integrate (G, xl, xu, G->dim, config->warmup, r, s,
			       &res, &err);
  display_results ("vegas warm-up", res, err);
  
  fprintf(stderr,"converging...\n");
  int tries=0;

  do {
    gsl_monte_vegas_integrate (G, xl, xu, G->dim, calls, r, s,
				 &res, &err);
    fprintf(stderr,"result = %.4g sigma = %.4g "
	    "chisq/dof = %.1f\n", res, err, s->chisq);
  } while (fabs(s->chisq - 1.0) > 0.5 && (tries++ < 10));
  

  display_results ("vegas final", res, err);
  
  gsl_monte_vegas_free (s);

  *_res=res;
  *_err=err;
}


void mc_dual_integrate(gsl_monte_function *G,
		       gsl_monte_function *F,
		       gsl_rng *r,
		       double *_res1,double *_err1,
		       double *_res2,double *_err2,
		       MCConfig *config) {

  double res, err;
  size_t i;

  double *xl = (double*)malloc(G->dim*sizeof(double));
  double *xu = (double*)malloc(G->dim*sizeof(double));


  for (i=0;i<G->dim;i++) {
    xl[i]=0.00001;
    xu[i]=1.0;
  }

  size_t calls = config->calls;

     
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (G->dim);
     
  gsl_monte_vegas_integrate (G, xl, xu, G->dim, config->warmup, r, s,
			       &res, &err);
  display_results ("vegas warm-up", res, err);
  
  fprintf(stderr,"converging...\n");
  int tries=0;

  do {
    gsl_monte_vegas_integrate (G, xl, xu, G->dim, calls, r, s,
				 &res, &err);
    fprintf(stderr,"A: result = %.4g sigma = %.4g "
	    "chisq/dof = %.1f\n", res, err, s->chisq);
  } while (fabs(s->chisq - 1.0) > 0.5 && (tries++ < 10));
  

  display_results ("A: vegas final", res, err);

  *_res1=res;
  *_err1=err;

  s->stage=1;
  do {
    gsl_monte_vegas_integrate (F, xl, xu, G->dim, calls, r, s,
				 &res, &err);
    fprintf(stderr,"B: result = %.4g sigma = %.4g "
	    "chisq/dof = %.1f\n", res, err, s->chisq);
  } while (fabs(s->chisq - 1.0) > 0.5 && (tries++ < 10));
  

  display_results ("B: vegas final", res, err);
  

  *_res2=res;
  *_err2=err;


  gsl_monte_vegas_free (s);

  free(xl);
  free(xu);
}
