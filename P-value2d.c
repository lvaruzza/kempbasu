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
#include "params.h"
#include "kempbasu.h"

int main (int argc,char **argv)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();
     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  unsigned int a=atoi(argv[1]);
  unsigned int b=atoi(argv[2]);
  double Pvalue=0.0;

  gsl_vector_uint *x=gsl_vector_uint_alloc(2);

  gsl_vector_uint_set(x,0,a);
  gsl_vector_uint_set(x,1,b);

  gsl_vector_uint *sum=gsl_vector_uint_alloc(2);

  gsl_vector_uint_set(sum,0,Na);
  gsl_vector_uint_set(sum,1,Nb);

  P_value(r,x,sum,&Pvalue,P_VALUE_DEFAULTS);

  printf ("\n");
  gsl_rng_free (r);
  return 0;
}
