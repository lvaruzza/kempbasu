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
#include <stdlib.h>

#include <gsl/gsl_randist.h>
#include "kempbasu.h"

void usage() {
  printf("Use: pvalue_ab u v\n");
}


int main(int argc,char *argv[]) {
  unsigned int i,u,v;

  if (argc < 2) {
    usage();
    exit(-1);
  }

  u=atoi(argv[1]);
  v=atoi(argv[2]);

  printf("#N\talpha\n");
  for(i=1;i<=10000;i+=1) {
    double alpha;
    pvalue_alpha_beta(i,u,v,&alpha,NULL);
    printf("%i\t%lf\n",i,alpha);
  }

  return 0;
}
