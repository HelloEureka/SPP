/*
 * CMat.h
 *
 *  Created on: 2018/3/12
 *      Author: doublestring
 */

#ifndef INCLUDE_COM_CMAT_H_
#define INCLUDE_COM_CMAT_H_

#include <stdio.h>

//// double precision part
static void
CMat_Matmul(const char *tr, int m, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta,
            double *c, int ldc);

static int CMat_Inverse(double *A, int lda, int n);

static void CMat_PrintMatrix(double *mat, int lda, int row, int col, const char *premsg);

double *mat(int n, int m);

static int *imat(int n, int m);

static void lubksb(const double *A, int n, const int *indx, double *b);

static int ludcmp(double *A, int n, int *indx, double *d);


void demo();

#endif /* INCLUDE_COM_CMAT_H_ */
