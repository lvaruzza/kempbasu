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
#include <unistd.h>
#include <gsl/gsl_matrix_uint.h>
#include <gsl/gsl_errno.h>
#include "matrix.h"


int matrix_read_uint(FILE *in,size_t rows,size_t columns,void *_matrix) {
  gsl_matrix_uint *mat=gsl_matrix_uint_alloc(rows,columns);

  gsl_matrix_uint_fscanf(in,mat);

  gsl_matrix_uint **result=_matrix;
  *result=mat;

  return 0;
}


int matrix_read(FILE *in,MatrixType type,void *matrix) {
  size_t rows;
  size_t columns;

  fscanf(in,"%zi %zi",&rows,&columns);

  switch(type) {
  case MAT_UINT:
    matrix_read_uint(in,rows,columns,matrix);
    break;
  default:
    GSL_ERROR("Invalid matrix format",GSL_EINVAL);
    break;
  }

  return 0;
}

int matrix_check_uint(void *_matrix) {
  size_t i,j;

  gsl_matrix_uint *m=(gsl_matrix_uint *)_matrix;
  gsl_vector_uint *s=gsl_vector_uint_alloc(m->size2);

  gsl_vector_set_zero(s);

  for(i=1;i<m->size1;i++) {
    for(j=0;j<m->size2;j++) { 
      gsl_vector_uint_set(s,j,gsl_vector_uint_get(s,j)+gsl_matrix_uint_get(m,i,j));
    }
  }
  
  for(j=0;j<m->size2;j++) { 
    if (gsl_vector_uint_get(s,j) > gsl_matrix_uint_get(m,0,j)) {
      char buff[128];
      sprintf(buff,
	      "Sum of values in column %i is %i,"
	      "which is greater then declared total %i",
	      j,
	      gsl_vector_get(s,j),
	      gsl_matrix_uint_get(m,0,j));
	
      fprintf(stderr,buff);
      GSL_ERROR(buff,GSL_EINVAL);
    }
  }
  return 1;
}

int matrix_check(MatrixType type,void *matrix) {
  switch(type) {
  case MAT_UINT:
    return matrix_check_uint(matrix);
    break;
  default:
    GSL_ERROR("Invalid matrix format",GSL_EINVAL);
    break;
  }

  return 0;
}

int matrix_print_uint(FILE *out,void *_matrix) {
  size_t i,j;

  gsl_matrix_uint *m=_matrix;

  fprintf(out,"%zi %zi\n",m->size1,m->size2);

  for(i=0;i<m->size1;i++) {
    for(j=0;j<m->size2;j++) { 
      fprintf(out,"%u\t",gsl_matrix_uint_get(m,i,j));
    }
    fputc('\n',out);
  }

  return 0;
}

int matrix_print(FILE *out,MatrixType type,void *matrix) {
  switch(type) {
  case MAT_UINT:
    return matrix_print_uint(out,matrix);
    break;
  default:
    GSL_ERROR("Invalid matrix format",GSL_EINVAL);
    break;
  }
}

