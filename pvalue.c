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

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector_uint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "kempbasu.h"

#define ELT gsl_vector_uint_get
#define SET_ELT gsl_vector_uint_set

PvalueConfig _PVALUE_DEFAULTS={ 
  1000000 // RUNS (1M)
};

PvalueConfig *PVALUE_DEFAULTS=&_PVALUE_DEFAULTS;



double logRV(size_t dims,gsl_vector_uint *x,gsl_vector_uint *sums,unsigned int NN) {
  unsigned int n=0;
  double r=0;

  size_t i;

  for(i=0;i<dims;i++) {
    n+=gsl_vector_uint_get(x,i);
  }

  // x*(log(Na)-log(x))+y*(log(Nb)-log(y))+n*(log(n)-log(NN));

  for(i=0;i<dims;i++) {
    unsigned int xi=gsl_vector_uint_get(x,i);
    if (xi != 0) {
      r+=xi*(log(ELT(sums,i))-log(xi));
    } else {
      r+=SMALL_P*(log(ELT(sums,i))-log(SMALL_P));
    }
  }

  if (n==0)
    return r;
  else 
    return r+n*(log(n)-log(NN));
}


void pvalue(gsl_rng * r,
	    gsl_vector_uint *x,
	    gsl_vector_uint *sums,
	    double *_pvalue,
	    double *_alpha,
	    double *_alpha_score,
	    double *_beta_score,
	    PvalueConfig *config) {

  fprintf(stderr,"pvalue\n");
  assert(x->size == sums->size);

  size_t i;
  size_t dim=x->size;
  unsigned int runs=config->runs;
  unsigned int N=0;
  unsigned int NN=0;

  for(i=0;i<dim;i++) {
    N+=ELT(x,i);
    NN+=ELT(sums,i);
  }
  printf("N=%i\n",N);

  double *p=(double*)malloc(sizeof(double)*dim);

  for(i=0;i<dim;i++) {
    p[i]=ELT(sums,i)*1.0/NN;
  }

  double cutoff = logRV(dim,x,sums,NN);
  fprintf(stderr,"cutoff = %f\n",cutoff);

  double rv;
  unsigned int positives=0;

  //unsigned int *n=(unsigned int*)malloc(sizeof(unsigned int)*dim);

  gsl_vector_uint *n=gsl_vector_uint_alloc(dim);

  for (i = 0; i < runs; i++)  {
    gsl_ran_multinomial (r, dim, N, p, n->data);
    rv=logRV(dim,n,sums,NN);
    if (rv <= cutoff) {
      positives++;
      /*fprintf(multi,"%i %i %i\t%i\n",n->data[0],n->data[1],n->data[2],
	n->data[0]+n->data[1]+n->data[2]);    */
    }
  }

  double pvalue=positives*1.0/runs;

  fprintf(stderr,"%i %i %f",positives,runs,pvalue);

  *_pvalue=pvalue;
  pvalue_alpha_beta(N,4,1,_alpha,_beta_score);

  *_alpha_score=1.0-pvalue/(*_alpha);

  fprintf (stderr,"\n");

  free(p);
  //free(n);
  gsl_vector_uint_free(n);
}


void pvalue_all_3d(gsl_rng * r,
		   gsl_vector_uint *x,
		   gsl_vector_uint *sums) {

  assert(x->size == sums->size);

  size_t i,j,k;

  size_t dim=x->size;
  unsigned int N=0;
  unsigned int NN=0;

  for(i=0;i<dim;i++) {
    N+=ELT(x,i);
    NN+=ELT(sums,i);
  }

  double *p=(double*)malloc(sizeof(double)*dim);

  for(i=0;i<dim;i++) {
    p[i]=((double)ELT(sums,i))/NN;
  }

  double cutoff = logRV(dim,x,sums,NN);
  fprintf(stderr,"cutoff = %f\n",cutoff);

  double rv;
  unsigned int positives=0;

  gsl_vector_uint *n=gsl_vector_uint_alloc(dim);


  FILE *graph_pos,*graph_neg;

  char buff[256];
  sprintf(buff,"pvalue_graph_pos-%i-%i-%i.dat",x->data[0],x->data[1],x->data[2]);

  if ( (graph_pos=fopen(buff,"w+"))==NULL) {
    fprintf(stderr,"ERROR: Can't open pvalue_graph_pos.dat");
    exit(-1);
  }

  sprintf(buff,"pvalue_graph_neg-%i-%i-%i.dat",x->data[0],x->data[1],x->data[2]);
  if ( (graph_neg=fopen(buff,"w+"))==NULL) {
    fprintf(stderr,"ERROR: Can't open pvalue_graph_neg.dat");
    exit(-1);
  }

  double pr;
  double pv=0;
  unsigned int positivos=0;
  unsigned int total=0;

  pr=gsl_ran_multinomial_pdf(dim,p,x->data);

  printf("pr(x) = %f\n",pr);

  for (i=0; i <= N;i+=1) {
    SET_ELT(n,0,i);
    for(j=0; i+j<= N;j+=1) { 
      SET_ELT(n,1,j);
      SET_ELT(n,2,N-(i+j));	

      rv=logRV(dim,n,sums,NN);
      pr=gsl_ran_multinomial_pdf(dim,p,n->data);
      //printf("%f %f\n",rv,cutoff);
	
      if (rv <= cutoff) { 
	  pv+=pr;
	  positives++;
	  fprintf(graph_pos,"%u %u %f\n",i,j,pr);
      } else {
	fprintf(graph_neg,"%u %u %f\n",i,j,pr);
      }
      total++;
    }
  }
  printf("pos = %u total = %u: ratio = %f\n",positives,total,((double)positives)/total);
  printf("pvalue2 = %f\n",pv);
  fclose(graph_pos);
  fclose(graph_neg);

  /*for (i = 0; i < runs; i++)  {
    gsl_ran_multinomial (r, dim, N, p, n->data);
    rv=logRV(dim,n,sums,NN);
    if (rv <= cutoff) 
      positives++;
    //fprintf(stderr,"%i: %i %i\t%f\t%i\n",i,n[0],n[1],rv,);    
    }*/

  free(p);
  //free(n);
  gsl_vector_uint_free(n);
}

