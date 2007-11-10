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
#include <math.h>
#include <stdlib.h>
#include "params.h"

#define NN (Na+Nb)
    
double logRV(unsigned int x,unsigned int y) {
  unsigned int n=x+y;

  double r=0.0;

  if (x!=0) {
    r+=x*(log(Na)-log(x));
  }

  if (y!=0) {
    r+=y*(log(Nb)-log(y));
  }

  return r+n*(log(n)-log(NN));
}


int main (int argc,char **argv)
{
  const gsl_rng_type * T;
  gsl_rng * r;
     
  int i, runs = 1000*1000;

  //unsigned int a=144;
  //unsigned int b=221;

  unsigned int a=atoi(argv[1]);
  unsigned int b=atoi(argv[2]);

  unsigned int N=a+b;
  unsigned int n[2];

  double p[]={Na*1.0/NN,Nb*1.0/NN};

  /* create a generator chosen by the 
     environment variable GSL_RNG_TYPE */
     
  gsl_rng_env_setup();
     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
     
  /* print n random variates chosen from 
     the poisson distribution with mean 
     parameter mu */

  double cutoff = logRV(a,b);
  printf("cutoff = %f\n",cutoff);

  double rv;
  unsigned int positives=0;

  for (i = 0; i < runs; i++)  {
    gsl_ran_multinomial (r, 2, N, p, n);
    rv=logRV(n[0],n[1]);
    if (rv < cutoff) 
      positives++;
    //printf("%i: %i %i\t%f\t%i\n",i,n[0],n[1],rv,);    
  }
  printf("%i %i %f",positives,runs,positives*1.0/runs);

  printf ("\n");
  gsl_rng_free (r);
  return 0;
}
