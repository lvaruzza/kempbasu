#include <stdio.h>
#include <gsl/gsl_sf_psi.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_matrix_uint.h>
#include <gsl/gsl_vector_uint.h>
#include <gsl/gsl_vector_double.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>

#ifdef HAVE_PTHREAD
#  include <pthread.h>
#  define PTHREADS_ARGS "n:"
//pthread_mutex_t results_mutex;
size_t nprocs=2;
#endif

#include "kempbasu.h"
#include "matrix.h"

#define ELT gsl_matrix_uint_get

int verbose_flag=1;
gsl_matrix_uint *matrix;

static struct option long_options[] =
  {
    /* These options set a flag. */
    {"verbose", no_argument,       &verbose_flag, 1},
    {"brief",   no_argument,       &verbose_flag, 0},
#ifdef HAVE_PTHREAD
    {"nprocs",  required_argument,       0, 'n'},
#endif
    {"help",    no_argument,       0, 'h'},
    {0, 0, 0, 0}
  };

static void help() {
  printf("Use: loginorm [-h|--help] [-n|--nprocs N] input_matrix\n");
  exit(0);
}

void read_cmd_line(int *argc,char ***argv) {
  int c;

  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    c = getopt_long (*argc, *argv, PTHREADS_ARGS "hn:",
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

void loginorm_params(gsl_matrix_uint *m,char *in_filename) {
  size_t i,j,k;
  double mean,var,psi_i,psi0;
  gsl_vector *psi=gsl_vector_alloc(m->size2);
  gsl_vector *psi1=gsl_vector_alloc(m->size2);
  size_t nvars=m->size2-1;
  FILE *tbl;
  FILE *norm;

  char buff[512];

  strncpy(buff,in_filename,256);
  strcat(buff,"-gamma");

  if ((tbl=fopen(buff,"w+"))==NULL) {
    fprintf(stderr,"Can't open file %s\n",buff);
    exit(-1);
  }

  strncpy(buff,in_filename,256);
  strcat(buff,"-norm");

  if ((norm=fopen(buff,"w+"))==NULL) {
    fprintf(stderr,"Can't open file %s\n",buff);
    exit(-1);
  }

  fprintf(norm,"%i %i\n",m->size1,nvars*(1+nvars));


  for(i=0;i<m->size1;i++) {
    fprintf(tbl,"X\t");
    for(j=0;j<m->size2;j++) {
      fprintf(tbl,"%i\t",ELT(m,i,j));
    }
    fprintf(tbl,"\n");

    for(j=0;j<m->size2;j++) {
      psi_i = gsl_sf_psi(ELT(m,i,j) + 0.5);
      gsl_vector_set(psi,j,psi_i);

      psi_i = gsl_sf_psi_1(ELT(m,i,j) + 0.5);
      gsl_vector_set(psi1,j,psi_i);
    }

    // Gamma-table
    fprintf(tbl,"digamma\t");
    for(j=0;j<m->size2;j++) {
      fprintf(tbl,"%.16lf\t",gsl_vector_get(psi,j));
    }
    fprintf(tbl,"\n");

    fprintf(tbl,"trigamma\t");
    for(j=0;j<m->size2;j++) {
      fprintf(tbl,"%.16lf\t",gsl_vector_get(psi1,j));
    }
    fprintf(tbl,"\n\n");

    for(j=1;j<m->size2;j++) {
      mean=gsl_vector_get(psi,j)-psi0;
      fprintf(norm,"%.16lf ",mean);
    }
    for(j=1;j<m->size2;j++) {
      for(k=1;k<m->size2;k++) {
	if (j==k) {
	  var=gsl_vector_get(psi1,j)+gsl_vector_get(psi1,0);
	  fprintf(norm,"%.16lf ",var);
	} else {
	  fprintf(norm,"%.16lf ",gsl_vector_get(psi1,0));
	}
      }
    }
    fprintf(norm,"\n");
  }
}


int main(int argc,char **argv) {
  size_t i,j;
  FILE *in;
  char *in_filename=NULL;
  char *out_filename=NULL;

  read_cmd_line(&argc,&argv);

  if (argc > 1) {
    in_filename=argv[1];
    if ((in=fopen(in_filename,"r"))==NULL) {
      fprintf(stderr,"ERROR: Can't open input file %s\n",in_filename);
      return -1;
    }
    char *buff=(char*)malloc(sizeof(char)*512);
    strncpy(buff,in_filename,256);
    //strcat(buff,suffix);
    out_filename=buff;
  } else {
    in=stdin;
    out_filename=NULL;
  }

  matrix_read(in,MAT_UINT,&matrix);

  loginorm_params(matrix,in_filename);
}
