#include <kempbasu.h>
#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix_uint.h>
#include <gsl/gsl_vector_uint.h>

unsigned int _sums[]={35076,73933,79020};
unsigned int _x[]={63,85,90};

PvalueConfig test_conf= { 1000*1000*1 };

int main(int argc,char **argv) {
  int i;

  const gsl_rng_type *RT;
  gsl_rng *r;
  gsl_rng_env_setup();
     
  RT = gsl_rng_default;
  r = gsl_rng_alloc (RT);
  
  double pv;
  double alpha;
  double beta_score,alpha_score;

  gsl_vector_uint *sums=gsl_vector_uint_alloc(3);
  gsl_vector_uint *x=gsl_vector_uint_alloc(3);
 
  for(i=0;i<3;i++) {
    gsl_vector_uint_set(sums,i,_sums[i]);
    gsl_vector_uint_set(x,i,atoi(argv[i+1]));
    printf("%u %u\n",gsl_vector_uint_get(x,i),gsl_vector_uint_get(sums,i));
  }

  pvalue(r,
	 x,sums,
	 &pv,
	 &alpha,
	 &alpha_score,
	 &beta_score,
	 &test_conf);

  printf("pvalue = %f alpha = %f\n",pv,alpha);


  printf("%u\n",x->size);

  pvalue_all_3d(r,x,sums);

  return 0;
}
