/*
 * CMat.h
 *
 *  Created on: 2018/3/12
 *      Author: doublestring
 */

#ifndef INCLUDE_COM_CMAT_H_
#define INCLUDE_COM_CMAT_H_

#include <stdio.h>
#include <math.h>

typedef struct {
    int row;
    int col;
    double *val;
} Mat;


//// double precision part
void CMat_Matmul(const char *tr, int m, int n, int k, double alpha, double *a, double *b, double beta, double *c);

void CMat_Multiply(double *left, double *right, double *res, int row, int mid, int col);

int CMat_Inverse(double *A, int lda, int n);

void CMat_Transpose(double *A, int row, int col);

void CMat_PrintMatrix(double *mat, int row, int col, const char *premsg);

double *mat(int n, int m);

static int *imat(int n, int m);

static void lubksb(const double *A, int n, const int *indx, double *b);

static int ludcmp(double *A, int n, int *indx, double *d);

void CMat_Copy(double *des, double *src, int col, int row);

void demo();

double Mat_Norm(double *A, int len);

void Mat_Add(double *A, double *B, int row, int col);

#endif /* INCLUDE_COM_CMAT_H_ */
