#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_uint.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#include <assert.h>
#include <stdio.h>

#include "kempbasu.h"
#include "fbst.h"

#define ELT gsl_vector_uint_get
#define ELTd gsl_vector_get

double normal_prod(double *x,gsl_vector *means,gsl_vector *vars) {
  size_t j;
  double prod=1.0;
  for(j=0;j<means->size;j++) {
    prod*=gsl_ran_gaussian_pdf(x[j]-ELTd(means,j),ELTd(vars,j));
  }
  return prod;
}

double gsl_vector_sum(gsl_vector *v) {
  size_t i;
  double sum=0.0;
  for(i=0;i<v->size;i++) {
    sum+=gsl_vector_get(v,i);
  }
  return sum;
}

// Calulate p(j)=Prod_(i!=j) sigma_i
void mutual_prod(gsl_vector *x,gsl_vector *prods) {
  assert(x->size == prods->size);
  size_t i,j;
  size_t n=x->size;
  double prod;

  printf("mutial prod: ");
  for(i=0;i<n;i++) {
    prod=1.0;
    for(j=0;j<n;j++) {
      if (i!=j) prod*=gsl_vector_get(x,j);
    }
    gsl_vector_set(prods,i,prod);
    printf("%le ",i,prod);
  }
  printf("\n");
}

double normal_null_maximum(gsl_vector *means,gsl_vector *vars) {
  assert(means->size == vars->size);

  gsl_vector *prods=gsl_vector_alloc(means->size);

  mutual_prod(vars,prods);
  double a,b;
  b=gsl_vector_sum(prods);
  gsl_vector_mul(prods,means);
  a=gsl_vector_sum(prods);
  printf("a=%lg b=%lg\n",a,b);
  gsl_vector_free(prods);

  return a/b;
}

typedef struct {
  double cutoff;
  gsl_vector *means;
  gsl_vector *vars;
} Params;

static double tangent_region_fun (double *p, size_t dim, void *_params) {
  Params *par=(Params*)_params;
  double val=normal_prod(p,par->means,par->vars);

  //return val;

  if (val > par->cutoff) {
    return val;
  } else {
    return 0;
  }
}


static double normal_fun (double *p, size_t dim, void *_params) {
  Params *par=(Params*)_params;
  double val=normal_prod(p,par->means,par->vars);

  return val;
}

static double max_fun (double x, void *_params) {
  Params *par=(Params*)_params;
  
  gsl_vector *p=gsl_vector_alloc(par->means->size);
  gsl_vector_set_all(p,x);

  double val=-normal_prod(p->data,par->means,par->vars);

  gsl_vector_free(p);

  return val;
}

double iteractive_max(Params *params,double a, double b) {
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;
  int status;
  int iter = 0, max_iter = 100;
  double m;

  printf("a=%lg b=%lg\n",a,b);

  printf("f(a)=%lg f(b)=%lg\n",max_fun(a,params),max_fun(b,params));

  F.function = &max_fun;
  F.params = (void*)params;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);

  printf ("using %s method\n",
          gsl_min_fminimizer_name (s));

  fflush(stdout);

  printf("Start Interactions\n");

  fflush(stdout);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status 
        = gsl_min_test_interval (a, b, 0.001, 0.0);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] "
              "%.7f %.7f\n",
              iter, a, b,
              m,  b - a);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free (s);

  return m;
}

void fbst_normal(gsl_rng *r,
		 gsl_vector *means,
		 gsl_vector *vars,
		 double *_ev,double *_err,
		 FBSTConfig *config) {

  assert(means->size == vars->size);

  size_t dim=means->size;

  if (config->mc_config==NULL) {
    config->mc_config=MC_DEFAULTS;
  }
  
  double xstar=normal_null_maximum(means,vars);
  gsl_vector *Xstar=gsl_vector_alloc(means->size);
  gsl_vector_set_all(Xstar,xstar);
  printf("x*=%lg\n",xstar);
  double cutoff=normal_prod(Xstar->data,means,vars);
  printf("cutoff=%lg\n",cutoff);

  gsl_vector_set_all(Xstar,xstar+0.1);
  double cutoff2=normal_prod(Xstar->data,means,vars);

  gsl_vector_set_all(Xstar,xstar-0.1);
  double cutoff3=normal_prod(Xstar->data,means,vars);

  printf("cutoff=%lg  +0.1=%lg  -0.1=%fl\n",cutoff,cutoff2,cutoff3);


  Params par;
  par.means=means;
  par.vars=vars;

  double a=xstar-10;
  double b=xstar+10;

  //double ixstar=iteractive_max(&par,0.0,1);
  //printf("icutoff=%lg cutoff=%lg\n",icutoff,cutoff);

  par.cutoff=cutoff;

  gsl_monte_function T = { *tangent_region_fun, dim, &par };
  double t,t_err;


  double *xl = (double*)malloc(dim*sizeof(double));
  double *xu = (double*)malloc(dim*sizeof(double));

  size_t i;

  for (i=0;i<dim;i++) {
    xl[i]=ELTd(means,i)-10*ELTd(vars,i);
    xu[i]=ELTd(means,i)+10*ELTd(vars,i);
  }

  //gsl_monte_function F = { *normal_fun, dim, &par };
  //double f,f_err;
  //mc_integrate1(&F,r,xl,xu,&f,&f_err,config->mc_config);

  //printf("I(f)=%lg+-%le\n",f,f_err);
  
  mc_integrate1(&T,r,xl,xu,&t,&t_err,config->mc_config);

  free(xl);
  free(xu);

  *_ev=1-t;
  *_err=t_err;

  printf("k*=%lg+-%le\n",t,t_err);
  printf("evidence=%lg+-%le\n",*_ev,*_err);
}

void fbst_loginorm(gsl_rng *r,
		   gsl_vector_uint *x,
		   gsl_vector_uint *sums,
		   double alpha,
		   double beta,
		   double *_ev,double *_err,
		   FBSTConfig *config) {

  size_t i,j;
  double mean;
  double var;
  double a,b;

  assert(x->size == sums->size);

  gsl_vector *means=gsl_vector_alloc(x->size);
  gsl_vector *vars=gsl_vector_alloc(x->size);

  
  for(i=0;i<x->size;i++) {
    a=ELT(x,i)+alpha;
    b=ELT(sums,i)+beta-a;
    mean=gsl_sf_psi(a)+gsl_sf_psi(b);
    var=gsl_sf_psi_1(a)+gsl_sf_psi_1(b);

    printf("x%i=N(%lg %lg)\n",i,mean,var);

    gsl_vector_set(means,i,mean);
    gsl_vector_set(vars,i,var);
  }

  fbst_normal(r,means,vars,_ev,_err,config);
  gsl_vector_free(means);
  gsl_vector_free(vars);  
}
