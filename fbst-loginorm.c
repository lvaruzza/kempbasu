#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_uint.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_blas.h>

#include <assert.h>
#include <stdio.h>

#include "kempbasu.h"
#include "fbst.h"

#define ELT gsl_vector_uint_get
#define ELTd gsl_vector_get

double normal_prod(double *x,gsl_vector *means,gsl_vector *vars) {
  size_t j;
  double lprod=0.0;
  double p;

  //printf("NP: x=%lf %lf\n",x[0],x[1]);

  for(j=0;j<means->size;j++) {
    p=gsl_ran_gaussian_pdf(x[j]-ELTd(means,j),ELTd(vars,j));

    lprod+=log(p);
  }

  return exp(lprod);
}

double gsl_vector_sum(gsl_vector *v) {
  size_t i;
  double sum=0.0;
  for(i=0;i<v->size;i++) {
    sum+=gsl_vector_get(v,i);
  }
  return sum;
}

void bhaskara(double a,double b,double c,double *x1,double *x2) {
  double delta=b*b-4.0*a*c;

  printf("delta=%lf\n",delta);

  if (delta > 0) {
    *x1=(-b-sqrt(delta))/(2.0*a);
    *x2=(-b+sqrt(delta))/(2.0*a);
  }  else
    if (delta == 0) {
      *x1=*x2=-b/(2.0*a);
    } else {
      printf("Bhaskara: Root not found\n");
      *x1=*x2=NAN;
    }
}

// Calulate p(j)=Prod_(i!=j) sigma_i
void mutual_prod(gsl_vector *x,gsl_vector *prods) {
  assert(x->size == prods->size);
  size_t i,j;
  size_t n=x->size;
  double prod;

  printf("mutual prod: ");
  for(i=0;i<n;i++) {
    prod=1.0;
    printf("(");
    for(j=0;j<n;j++) {
      printf("%lf ",gsl_vector_get(x,j));
      if (i!=j) prod*=gsl_vector_get(x,j);
    }
    gsl_vector_set(prods,i,prod);
    printf(")=%lg ",gsl_vector_get(prods,i));
  }
  printf("\n");
}

double normal_null_maximum(gsl_vector *means,gsl_vector *s) {
  assert(means->size == s->size);
  
  // Baskara coefs.
  double a,b,c;

  gsl_vector *s2=gsl_vector_alloc(means->size);
  gsl_blas_dcopy(s,s2);
  gsl_vector_mul(s2,s2);

  gsl_vector *B=gsl_vector_alloc(means->size);
  gsl_blas_dcopy(means,B);
  gsl_blas_dscal(-2.0,B);
  printf("B=%lg %lg\n",ELTd(B,0),ELTd(B,1));

  gsl_vector *C=gsl_vector_alloc(means->size);
  gsl_blas_dcopy(means,C);
  gsl_vector_mul(C,C);

  gsl_vector *prods=gsl_vector_alloc(means->size);
  mutual_prod(s2,prods);

  a=gsl_blas_dasum(prods);
  gsl_blas_ddot(B,prods,&b);
  gsl_blas_ddot(C,prods,&c);

  printf("null max: a=%lf b=%lf c=%lf\n",a,b,c);

  double delta=b*b-4*a*c;
  double x0=-b/(2.0*a);
  printf("null max: delta=%lg\n",delta);

  if (delta > 0) {
    if (delta < 1e-5) {
      return x0; 
    } else {
      double x1=(-b-sqrt(delta))/(2.0*a);
      double x2=(-b-sqrt(delta))/(2.0*a);
      
      printf("null max: x1=%lg x2=%lg\n",x1,x2);
      return x1;
    }
  } else {
    printf("WARNING: Null max not found!");


    return 0.0;
  }
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

  printf("f(a)=%lg f(b)=%lg f((a+b)/2)=%lg\n",max_fun(a,params),max_fun(b,params),max_fun((a+b)/2.0,params));

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
  double max=normal_prod(means->data,means,vars);
  printf("cutoff=%lg max=%lg\n",cutoff,max);

  gsl_vector_set_all(Xstar,xstar+ELTd(vars,0)/10.0);
  double cutoff2=normal_prod(Xstar->data,means,vars);

  gsl_vector_set_all(Xstar,xstar-ELTd(vars,0)/10.0);
  double cutoff3=normal_prod(Xstar->data,means,vars);

  printf("cutoff=%lg  +s/10=%lg  -s/10=%fl\n",cutoff,cutoff2,cutoff3);


  double *xl = (double*)malloc(dim*sizeof(double));
  double *xu = (double*)malloc(dim*sizeof(double));

  size_t i;
  double min_xl=1e30;
  double max_xu=-1e30;
  double xl1,xu1;

  for (i=0;i<dim;i++) {
    xl[i]=ELTd(means,i)-10*ELTd(vars,i);
    xu[i]=ELTd(means,i)+10*ELTd(vars,i);
    
    xl1=ELTd(means,i)-10*ELTd(vars,i);
    xu1=ELTd(means,i)+10*ELTd(vars,i);

    if (xl1 < min_xl) min_xl=xl1;
    if (xu1 > max_xu) max_xu=xu1;
  }
  printf("min_xl=%lg max_xu=%lg\n",min_xl,max_xu);

  Params par;
  par.means=means;
  par.vars=vars;

  /*double ixstar=iteractive_max(&par,min_xl,max_xu);
  double icutoff=normal_prod(Xstar->data,means,vars);
  printf("ix* = %lg\n",ixstar);

  gsl_vector_set_all(Xstar,ixstar);
  printf("icutoff=%lg cutoff=%lg\n",icutoff,cutoff);
  */

  par.cutoff=cutoff;

  gsl_monte_function T = { *tangent_region_fun, dim, &par };
  double t,t_err;


  /*gsl_monte_function F = { *normal_fun, dim, &par };
  double f,f_err;
  mc_integrate1(&F,r,xl,xu,&f,&f_err,config->mc_config);
  printf("I(f)=%lg+-%le\n",f,f_err);*/
  
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
