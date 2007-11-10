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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "kempbasu.h"

#define ELT gsl_vector_uint_get

#define DUAL_INTEGRATE

FBSTConfig _FBST_DEFAULTS={NULL};
FBSTConfig *FBST_DEFAULTS=&_FBST_DEFAULTS;

inline double square(double x) {
  return x*x;
}
    
/*
 * log likelihood = Sum x_i ( log(N_i) + log(p_i) ) - ( Sum x_i) * log((Sum N_i * p_i))
 */



typedef struct {
  double cutoff;
  double likenorm;
  gsl_vector_uint *x;
  gsl_vector_uint *sums;
  unsigned int k;
} Params;

static double g (double *p, size_t dim, void *_params) {
  Params *params=(Params*)_params;

  double l1=0.0;
  double l2=0.0;
  size_t i;

  //a*(log(Na)+log(k[0]))+b*(log(Nb)+log(k[1]))+(-a-b)*log(Na*k[0]+Nb*k[1]);
  
  for(i=0;i<dim;i++) {
    unsigned int x=ELT(params->x,i);
    unsigned int sum=ELT(params->sums,i);

    l1 += x *(log(sum) + log(p[i]));
    l2 += sum*p[i]; 
    //fprintf(stderr,"%i %i %g %g\n",x,sum,l1,l2);
  }

  return l1 - params->k * log(l2);
}


static double f (double *p, size_t dim, void *_params) {
  Params *params=(Params*)_params;
  double ll=g(p,dim,params);
  return exp(ll-params->likenorm);
}

static double t (double *p, size_t dim, void *_params) {
  Params *params=(Params*)_params;
  double val=g(p,dim,params);

  if (val > params->cutoff) {
    return exp(val-params->likenorm);
  } else {
    return 0;
  }
}
     
void fbst(gsl_rng *r,
	  gsl_vector_uint *x,
	  gsl_vector_uint *sums,
	  double *_ev,double *_err,
	  FBSTConfig *config) {

  fprintf(stderr,"FBST\n");

  assert(x->size == sums->size);

  if (config->mc_config==NULL) {
    config->mc_config=MC_DEFAULTS;
  }

  size_t dim = x->size;

  fprintf(stderr,"dim = %zi\n",dim);

  Params *params=(Params*)malloc(sizeof(Params));
   
  gsl_monte_function G = { &f, dim, params };
  gsl_monte_function T = { &t, dim, params };

  size_t i;

  // a*log(Na)+b*log(Nb)-(a+b)*log(Na+Nb);

  double cutoff1=0.0;
  unsigned int k=0;
  unsigned int NN=0;

  for(i=0;i<dim;i++) {
    k  += ELT(x,i);
    NN += ELT(sums,i);
    cutoff1+=ELT(x,i)*log(ELT(sums,i));
  }

  // Verificar cutoff = NaN
  fprintf(stderr,"cutoff1 = %.5g\n",cutoff1);

  double cutoff=cutoff1-k*log(NN);
  fprintf(stderr,"cutoff  = %.5g\n",cutoff);

  params->cutoff=cutoff;
  params->k=k;
  params->x=x;

  params->sums=sums;

  double *phat=(double*)malloc(dim*sizeof(double));

  for(i=0;i<dim;i++) {
    // Add 1 to eliminet phat[i] = 0
    if (ELT(x,i) != 0) {
      phat[i]=ELT(x,i)*1.0/ELT(sums,i);
    } else {
      phat[i]=SMALL_P;
    }
  }


  double likenorm = g(phat,dim,params);
  params->likenorm = likenorm;

  fprintf(stderr,"likenorm = %.5g\n",likenorm);

  if (isnan(likenorm)) {
    likenorm = 0;
  }

  double g,g_err;
  double t,t_err;

#ifdef DUAL_INTEGRATE
  mc_dual_integrate(&G,&T,r,&g,&g_err,&t,&t_err,config->mc_config);
#else
  mc_integrate(&G,r,&g,&g_err,config->mc_config);

  fprintf(stderr,"null = %.5g\n",g);

  mc_integrate(&T,r,&t,&t_err,config->mc_config);
#endif


  *_err=t/g*sqrt(square(g_err/g)+square(t_err/t));
  *_ev=1.0-t/g;

  fprintf(stderr,"evidence = %.5g\n",*_ev);

  free(phat);
  free(params);
}
