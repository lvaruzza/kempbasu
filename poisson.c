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
     
int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;
     
  int i, n = 10;
  double mu = 3.0;
     
  /* create a generator chosen by the 
     environment variable GSL_RNG_TYPE */
     
  gsl_rng_env_setup();
     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
     
  /* print n random variates chosen from 
     the poisson distribution with mean 
     parameter mu */
     
  for (i = 0; i < n; i++) 
    {
      unsigned int k = gsl_ran_poisson (r, mu);
      printf (" %u", k);
    }
     
  printf ("\n");
  gsl_rng_free (r);
  return 0;
}

