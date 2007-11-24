#ifndef FBST_H
#define FBST_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector_uint.h>

/*
 * FBST
 *
 */


typedef struct _FBSTConfig {
  MCConfig *mc_config;
} FBSTConfig;

extern FBSTConfig *FBST_DEFAULTS;

void fbst_poisson(gsl_rng *r,
		  gsl_vector_uint *x,
		  gsl_vector_uint *sums,
		  double alpha,
		  double beta,
		  double *_ev,double *_err,
		  FBSTConfig *config);


void fbst_loginorm(gsl_rng *r,
		   gsl_vector_uint *x,
		   gsl_vector_uint *sums,
		   double alpha,
		  double beta,
		   double *_ev,double *_err,
		   FBSTConfig *config);

#endif
