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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "kempbasu.h"
#define ELT gsl_vector_uint_get

P_ValueConfig _P_VALUE_DEFAULTS={NULL,100};
P_ValueConfig *P_VALUE_DEFAULTS=&_P_VALUE_DEFAULTS;

MCConfig _MC_DEFAULTS_PVALUE={ 
  5000,   // warmup
  100*1000  // calls
};

MCConfig *MC_DEFAULTS_P_VALUE=&_MC_DEFAULTS_PVALUE;

typedef struct {
  double cutoff;
  double likenorm;
  gsl_vector_uint *x;
  gsl_vector_uint *sums;
  unsigned int k;
} Params;

static double g (double *p, size_t dim, void *_params) {
  Params *params=(Params*)_params;
  
  size_t i;
  double l1=0.0;
  double l2=0.0;

  for(i=0;i<dim;i++) {
    // xi*(log(Ni)+log(Pi))
    l1+=ELT(params->x,i)*(log(ELT(params->sums,i))+log(p[i]));

    // Ni*Pi
    l2+=ELT(params->sums,i)*p[i];
  }

  //printf("%lf %lf\n",l1,l2);

  return l1-params->k*log(l2);
}


static double f (double *p, size_t dim, void *_params) {
  Params *params=(Params*)_params;
  double ll=g(p,dim,params);
  return exp(ll-params->likenorm);
}

double FB(gsl_rng *r,size_t dims,gsl_vector_uint *x,
	  gsl_vector_uint *sums,unsigned int NN,MCConfig *config) {

  
  size_t i;
  unsigned int k=0;
  double *p=malloc(sizeof(double)*dims);

  for(i=0;i<dims;i++) {
    k+=ELT(x,i);
  }

  for(i=0;i<dims;i++) {
    unsigned int xx=ELT(x,i);
    if (xx==0) {
      p[i]=0.0001;
    } else {
      p[i]=xx*1.0/k;    
    }
  }

  Params params;

  params.sums=sums;
  params.x=x;
  params.k=k;

  /* Normalizing factor, to avoid numeric undeflow */

  params.likenorm = g(p,dims,&params);

  printf(" likenorm = %lf\n",params.likenorm);


  gsl_monte_function G = { &f, dims, &params };

  double g;
  double g_err;

  mc_integrate(&G,r,&g,&g_err,config);

  printf("g=%g error=%lf\n",g,g_err);
  free(p);

  double h=0.0;

  for(i=0;i<dims;i++) {
    h+=ELT(x,i)*log(ELT(sums,i));
  }
  h-=k*log(NN);

  printf("h=%lf\n",h);

  double fb=exp(h+params.likenorm)/g;

  printf("fb=%g\n",fb);


  return fb;
}

void P_value(gsl_rng *r,
	     gsl_vector_uint *x,
	     gsl_vector_uint *sums,
	     double *_pvalue,
	     P_ValueConfig *config) {


  printf("P-value\n");

  assert(x->size == sums->size);

  if (config->mc_config==NULL) {
    config->mc_config=MC_DEFAULTS_P_VALUE;
  }


  size_t i;
  size_t dim=x->size;
  unsigned int runs=config->runs;
  unsigned int k=0;
  unsigned int NN=0;

  for(i=0;i<dim;i++) {
    k+=ELT(x,i);
    NN+=ELT(sums,i);
  }

  double *p=(double*)malloc(sizeof(double)*dim);

  for(i=0;i<dim;i++) {
    p[i]=ELT(sums,i)*1.0/NN;
  }

  double cutoff = FB(r,dim,x,sums,NN,config->mc_config);
  printf("cutoff = %f\n",cutoff);


  double fb;
  unsigned int positives=0;

  //unsigned int *n=(unsigned int*)malloc(sizeof(unsigned int)*dim);

  gsl_vector_uint *n=gsl_vector_uint_alloc(dim);

  for (i = 0; i < runs; i++)  {
    gsl_ran_multinomial (r, dim, k, p, n->data);
    fb=FB(r,dim,n,sums,NN,config->mc_config);
    printf("fb=%g cutoff=%g\n",fb,cutoff);
    if (fb < cutoff) {
      positives++;
      printf("=======> postitive: %i\n",positives);
    } else { 
      printf("<======= negative\n");
    }

    //printf("%i: %i %i\t%f\t%i\n",i,n[0],n[1],rv,);    
  }
  double pvalue=positives*1.0/runs;

  printf("PVALUE: %i %i %f",positives,runs,pvalue);

  *_pvalue=pvalue;

  printf ("\n");

  free(p);
  //free(n);
  gsl_vector_uint_free(n);
}

