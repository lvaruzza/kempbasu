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
#include <gsl/gsl_matrix_uint.h>
#include <gsl/gsl_vector_uint.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <Judy.h>

#ifdef HAVE_PTHREAD
#  include <pthread.h>
#endif

#include "kempbasu.h"
#include "matrix.h"


// Hehe, how to put to programs in a sigle c file

#ifdef KEMP
static int calculate_evalue = 0;
static int calculate_pvalue = 1;
static char *suffix="-kemp";
#define BINARY "kemp"
#endif

#ifdef BASU
static int calculate_evalue = 1;
static int calculate_pvalue = 0;
static char *suffix="-basu";
#define BINARY "basu"
#endif 


/************
 *
 * CMD LINE
 *
 ************/

int verbose_flag=1;
double prior_alpha=1.0;
double prior_beta=19.0;

#ifdef HAVE_PTHREAD
pthread_mutex_t results_mutex;
size_t nprocs=2;
#endif

static void help() {
  printf("Use: " BINARY "[-h|--help] [-n|--nprocs N] input_matrix\n");
  exit(0);
}

#ifdef HAVE_PTHREAD
#  define PTHREADS_ARGS "n:"
#endif

static struct option long_options[] =
  {
    /* These options set a flag. */
    {"verbose", no_argument,       &verbose_flag, 1},
    {"brief",   no_argument,       &verbose_flag, 0},
#ifdef HAVE_PTHREAD
    {"nprocs",  required_argument,       0, 'n'},
#endif
    {"help",    no_argument,       0, 'h'},

    {"prior-alpha",    no_argument,       0, 'a'},
    {"prior-beta",    no_argument,       0, 'b'},
    {0, 0, 0, 0}
  };

void read_cmd_line(int *argc,char ***argv) {
  int c;

  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    c = getopt_long (*argc, *argv, PTHREADS_ARGS "ha:b:n:",
		     long_options, &option_index);
    
    /* Detect the end of the options. */
    if (c == -1)
      break;
    
    switch (c)
      {
      case 0:
	/* If this option set a flag, do nothing else now. */
	if (long_options[option_index].flag != 0)
	  break;
	printf ("option %s", long_options[option_index].name);
	if (optarg)
	  printf (" with arg %s", optarg);
	printf ("\n");
	break;
	
      case 'h':
	help();
	break;

      case 'a':
	prior_alpha=atof(optarg);
	break;

      case 'b':
	prior_beta=atof(optarg);
	break;

#ifdef HAVE_PTHREAD
      case 'n':
	nprocs=atoi(optarg);
	fprintf(stderr,"Running %i proccess\n",nprocs);
	break;
#endif
      default:
	abort ();
      }
  }    

  *argc -= optind-1;
  *argv += optind-1;
}



/***************************************************************************
 Global results and value, like in the good old time of fortran programs
  (it's not the most elegant solution, but is good for the pthread code) 
****************************************************************************/

typedef struct RESULT {
  double pvalue;
  double pvalue_alpha;
  double pvalue_alpha_score;
  double pvalue_beta_score;

  double evalue;
  double evalue_sigma;
  double evalue_alpha;
  double evalue_alpha_score;
  double evalue_beta_score;
} Result;


Result **results; 
Pvoid_t hash;
gsl_matrix_uint *matrix;


void calculate_index(gsl_rng * r,gsl_vector_uint *x,gsl_vector_uint *sums,size_t i) 
{  
  size_t j;

  Result *value;
  PWord_t ptr;

  
  fprintf(stderr,"====================\n");
  fprintf(stderr,"V%zu: ",i);
  for(j=0;j<x->size;j++)
    fprintf(stderr,"%u\t",gsl_vector_uint_get(x,j));
  fprintf(stderr,"\n");
  fprintf(stderr,"====================\n");
  
  fprintf(stderr,"Searching in cache... ");

  
#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&results_mutex);
#endif

  JHSI(ptr,hash,x->data,x->size*sizeof(unsigned int));
  
  if (*ptr) {

#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&results_mutex);
#endif

    fprintf(stderr,"Using calculed value. ptr = %lu (%p)\n",*ptr,ptr);
    results[i]=results[*ptr];
    value=results[i];

  } else {

#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&results_mutex);
#endif

    double ev=0,ev_err=0;
    double pv=0;
    double pvalue_alpha=1.0;
    double alpha_score=1.0;
    double beta_score=1.0;

    value=(Result*)malloc(sizeof(Result));
    
    if (calculate_pvalue) {
      pvalue(r,
	     x,sums,
	     &pv,
	     &pvalue_alpha,
	     &alpha_score,
	     &beta_score,
	     PVALUE_DEFAULTS);
      
      
      value->pvalue=pv;
      value->pvalue_alpha=pvalue_alpha;
      value->pvalue_alpha_score=alpha_score;
      value->pvalue_beta_score=beta_score;
    }
    
    if (calculate_evalue) {
      int tries=1;
      do {
	fprintf(stderr,"====> EVALUE: %i try\n",tries);
	fbst(r,
	     x,
	     sums,
	     prior_alpha,
	     prior_beta,
	     &ev,&ev_err,
	     FBST_DEFAULTS);
	
	tries++;
	// Try again on non-sense results
      } while (ev < -0.01 && tries <= 10);
	
      if (ev < -0.01) {
	// if still negative is a problematic one
	ev=1.0/0.0;
      }
	
      // Remove small negative numbers (integration errors)
      if (ev < 0.0) {
	ev=0.0;
      }
      
      value->evalue=ev;
      value->evalue_sigma=ev_err;
    }


#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&results_mutex);
    JHSI(ptr,hash,x->data,x->size*sizeof(unsigned int));
    results[i]=value;
#endif

    *ptr=i;

#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&results_mutex);
#endif

  }
}


void *worker(void *_params) {
  size_t k=(size_t)_params;
  size_t i;


  /* RNG setup */
  const gsl_rng_type *RT;
  gsl_rng *r;
  gsl_rng_env_setup();
     
  RT = gsl_rng_default;
  r = gsl_rng_alloc (RT);

  gsl_vector_uint_view sums=gsl_matrix_uint_row(matrix,0);

  for(i=1+k;i<matrix->size1;i+=nprocs) {
    gsl_vector_uint_view x=gsl_matrix_uint_row(matrix,i);

    calculate_index(r,&x.vector,&sums.vector,i);
  }

  gsl_rng_free (r);

  return NULL;  
}

void save_results(char *out_filename) {
  FILE *out;

  size_t i;

  //Open output file
  if (out_filename) {
    if ((out=fopen(out_filename,"w+"))==NULL) {
      fprintf(stderr,"ERROR: Can't save %s. Priting results to stdout",out_filename);
      out=stdout;
    }
  } else {
    out=stdout;
  }


  // Output Header
  if (calculate_pvalue) {
    fprintf(out,"pvalue\tscore\tcategory");
    if (calculate_evalue) {
      fprintf(out,"\t");
    }
  }
  
  if (calculate_evalue) {
    fprintf(out,"evalue\tev ie");
  }
  fprintf(out,"\n");
  
  // Write results
  for(i=1;i<matrix->size1;i++) {    
    if (calculate_pvalue) {
      char *category=NULL;
      double score=0.0;

      double pvalue=results[i]->pvalue;
      double alpha=results[i]->pvalue_alpha;

      if (pvalue < alpha) {
	 category="D";
	 score=10.0*results[i]->pvalue_alpha_score;
      }

      if (pvalue >= alpha) {
	category="U";
	score=10.0*results[i]->pvalue_beta_score;
      }

      fprintf(out,"%g\t%g\t%s",
	      pvalue,
	      score,
	      category);

      if (calculate_evalue) {
	fprintf(out,"\t");
      }
    }

    if (calculate_evalue) {
      fprintf(out,"%.5g\t%.5g",
	      results[i]->evalue,
	      results[i]->evalue_sigma);
    }
    fprintf(out,"\n");
  }

  if(out_filename) {
    fclose(out);
  }
}

int main (int argc,char **argv) {
  size_t i,j;
  FILE *in;
  char *in_filename=NULL;
  char *out_filename=NULL;


  /*************
   *
   * SETUP
   *
   **************/
  

  read_cmd_line(&argc,&argv);

  /* open input file */

  if (argc > 1) {
    in_filename=argv[1];
    if ((in=fopen(in_filename,"r"))==NULL) {
      fprintf(stderr,"ERROR: Can't open input file %s\n",in_filename);
      return -1;
    }
    char *buff=(char*)malloc(sizeof(char)*512);
    strncpy(buff,in_filename,256);
    strcat(buff,suffix);
    out_filename=buff;
  } else {
    in=stdin;
    out_filename=NULL;
  }
  


  /**************
   *
   * READ 
   *
   ***************/

  
  matrix_read(in,MAT_UINT,&matrix);
  //matrix_print(stderr,MAT_UINT,matrix);

  gsl_vector_uint_view sums=gsl_matrix_uint_row(matrix,0);

  for(i=0;i < sums.vector.size;i++) {
    if (gsl_vector_uint_get(&sums.vector,i)==0) {
      fprintf(stderr,"Invalid zero sum at column %zi\n",i);
      return -1;
    }
  }

  fprintf(stderr,"M: %zu %zu\n",matrix->size1,matrix->size2);

  results=(Result**)malloc(sizeof(void*)*matrix->size1);


  /******************
   * Check Dataset 
   *
   ******************/
  
  if (matrix_check(MAT_UINT,matrix)) {
    fprintf(stderr,"ERROR: Invalid input matrix\n");
    exit(-1);
  }

  /******************
   *
   * Calculations 
   *
   *******************/

#ifdef HAVE_PTHREAD
  pthread_t *threads=(pthread_t*)malloc(sizeof(pthread_t) * nprocs);
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  size_t k;
  int rc;

  for(k=0;k<nprocs;k++) {
    fprintf(stderr,"Starting thread %i\n",k);
    rc = pthread_create(&threads[k], &attr, worker, (void *)k);
    if (rc){
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  /* Wait for all threads to complete */
  for (i=0; i<nprocs; i++) {
    pthread_join(threads[i], NULL);
  }
  free(threads);
#else
  worker((void*)0);
#endif


  /******************
   *
   *  Save results
   *
   *******************/

  save_results(out_filename);
  free(out_filename);

  /*
   * Finish
   *
   */
  gsl_matrix_uint_free(matrix);

  return 0;
}
