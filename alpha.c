#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*
 * Alpha is the cutoff for the pvalue
 *
 */
void pvalue_alpha_beta(unsigned int n,unsigned int u,unsigned int v,
		       double *_alpha,double *_beta_score) {
  unsigned int i,j,k;
  unsigned int middle=n/2+n%2;
  unsigned int rest=n%2;
  double error0=u+v;
  double alpha0=1.0;
  double beta0=1.0;
  double beta_score=1.0;
  unsigned int k0=1;

  double alpha,beta,error;

  //printf("n = %i middle = %i\n",n,middle);

  for(j=0;j<n/2;j++) {
    int a=middle-j-rest;
    int b=middle+j;
    //printf("%i: ",j);
    double sum=0.0;
    for(k=a;k<=b;k++) {
      double p=gsl_ran_binomial_pdf(k,0.5,n);
      //printf("[%i,%f] ",k,p);
      sum+=p;
    }
    //printf("\n");
    alpha=1-sum;
    k=(b-a+1);
    beta=(b-a+1)*1.0/(n+1);
    //printf("\n   alpha = %lf + beta = %lf = %lf\n",alpha,beta,alpha+beta);    
    error=u*alpha+v*beta;
    beta_score=1.0*k0/k;
    k0=k;

    if (error > error0) {
      break;
    }
    //printf("%u %f\n",j,error);
    error0=error;
    alpha0=alpha;
    beta0=beta;
  }
  if (_alpha != NULL)
    *_alpha=alpha0;

  if (_beta_score != NULL)
    *_beta_score=beta_score;
}
