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


#ifndef SIDE_H
#define SIDE_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector_uint.h>
#include <gsl/gsl_monte.h>

/*
 * Pvalue
 *
 */

#define SMALL_P 0.0001

typedef struct _PvalueConfig {
  unsigned int runs;
} PvalueConfig;

extern PvalueConfig *PVALUE_DEFAULTS;


void pvalue(gsl_rng * r,
	    gsl_vector_uint *x,
	    gsl_vector_uint *sums,
	    double *_pvalue,
	    double *_alpha,
	    double *_alpha_score,
	    double *_beta_score,
	    PvalueConfig *config);


void pvalue_all_3d(gsl_rng * r,
		   gsl_vector_uint *x,
		   gsl_vector_uint *sums);


/*
 * Monte Carlo
 *
 */

typedef struct _MCConfig {
  unsigned int warmup;
  unsigned int calls;
} MCConfig;

extern MCConfig *MC_DEFAULTS;
extern MCConfig *MC_DEFAULTS_P_VALUE;

void mc_integrate(gsl_monte_function *G,
		  gsl_rng *r,
		  double *_res,double *_err,
		  MCConfig *config);

/*
 * FBST
 *
 */


typedef struct _FBSTConfig {
  MCConfig *mc_config;
} FBSTConfig;

extern FBSTConfig *FBST_DEFAULTS;

void fbst(gsl_rng *r,
	  gsl_vector_uint *x,
	  gsl_vector_uint *sums,
	  double *_ev,double *_err,
	  FBSTConfig *config);

/*
 *
 * P-value
 *
 */



typedef struct _P_valueConfig {
  MCConfig *mc_config;
  unsigned int runs;
} P_ValueConfig;

extern P_ValueConfig *P_VALUE_DEFAULTS;

void P_value(gsl_rng *r,
	     gsl_vector_uint *x,
	     gsl_vector_uint *sums,
	     double *_ev,
	     P_ValueConfig *config);

#endif
