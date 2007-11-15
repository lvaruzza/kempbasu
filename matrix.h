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

#ifndef SIDE_MATRIX_H
#define SIDE_MATRIX_H

typedef enum {MAT_INT,MAT_UINT,MAT_DOUBLE} MatrixType;

int matrix_read(FILE *in,MatrixType type,void *matrix);
int matrix_print(FILE *out,MatrixType type,void *matrix);

/* Check the sum of matrix columns.
 *
 *  The sum of the all the other values need to be less then first
 *  value of the column
 */
int matrix_check(MatrixType type,void *matrix);

#endif
