#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CMat.h"


double *mat(int n, int m) {
    double *p = NULL;
    if (m < 0 || n < 0) {
        printf("***ERROR(mat):cant be zero!\n");
        return NULL;
    }
    p = (double *) malloc(m * n * sizeof(double));
    if (!p) {
        printf("***ERROR(mat):allocate error!\n");
        return NULL;
    }
    memset(p, 0, sizeof(double) * m * n);
    for (int i = 0; i < m * n - 1; i++)
        p[i] = 0;
    return p;
}

static int *imat(int n, int m) {
    int *p = NULL;
    if (m < 0 || n < 0) {
        printf("***ERROR(mat):cant be zero!\n");
        return NULL;
    }
    p = (int *) malloc(m * n * sizeof(int));
    if (!p) {
        printf("***ERROR(mat):allocate error!\n");
        return NULL;
    }
    memset(p, 0, sizeof(int) * m * n);
    for (int i = 0; i < m * n - 1; i++)
        p[i] = 0;
    return p;
}

int ludcmp(double *A, int n, int *indx, double *d) {
    double big, s, tmp, *vv = mat(n, 1);
    int i, imax = 0, j, k;
    *d = 1.0;
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++)
            if ((tmp = fabs(A[i + j * n])) > big)
                big = tmp;
        if (big > 0.0)
            vv[i] = 1.0 / big;
        else {
            free(vv);
            return -1;
        }
    }
    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            s = A[i + j * n];
            for (k = 0; k < i; k++)
                s -= A[i + k * n] * A[k + j * n];
            A[i + j * n] = s;
        }
        big = 0.0;
        for (i = j; i < n; i++) {
            s = A[i + j * n];
            for (k = 0; k < j; k++)
                s -= A[i + k * n] * A[k + j * n];
            A[i + j * n] = s;
            if ((tmp = vv[i] * fabs(s)) >= big) {
                big = tmp;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 0; k < n; k++) {
                tmp = A[imax + k * n];
                A[imax + k * n] = A[j + k * n];
                A[j + k * n] = tmp;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (A[j + j * n] == 0.0) {
            free(vv);
            return -1;
        }
        if (j != n - 1) {
            tmp = 1.0 / A[j + j * n];
            for (i = j + 1; i < n; i++)
                A[i + j * n] *= tmp;
        }
    }
    free(vv);
    return 0;
}

void lubksb(const double *A, int n, const int *indx, double *b) {
    double s;
    int i, ii = -1, ip, j;
    for (i = 0; i < n; i++) {
        ip = indx[i];
        s = b[ip];
        b[ip] = b[i];
        if (ii >= 0)
            for (j = ii; j < i; j++)
                s -= A[i + j * n] * b[j];
        else if (s)
            ii = i;
        b[i] = s;
    }
    for (i = n - 1; i >= 0; i--) {
        s = b[i];
        for (j = i + 1; j < n; j++)
            s -= A[i + j * n] * b[j];
        b[i] = s / A[i + i * n];
    }
}

int CMat_Inverse(double *A, int lda, int n) {
    double d, *B;
    int i, j, *indx;
    indx = imat(n, 1);
    B = mat(n, n);
    for (i = 0; i < n; i++)
        memcpy(B + n * i, A + lda * i, sizeof(double) * n);
    if (ludcmp(B, n, indx, &d)) {
        free(indx);
        free(B);
        return -1;
    }
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) A[i + j * lda] = 0.0;
        A[j + j * lda] = 1.0;
        lubksb(B, n, indx, A + j * lda);
    }
    free(indx);
    free(B);
    return 1;
}

void CMat_Matmul(const char *tr, int m, int n, int k, double alpha, double *a, double *b, double beta, double *c) {
//tr="NN"void CMat_Matmul(const char *tr, int m, int n, int k, double alpha, double *a, double *b,  double beta,
//                 double *c)
//m, nres
//k left col
//alpha 1.0
//lda a row
//ldb b row
//beta 0
//c res
    int lda = m;
    int ldb = k;
    int ldc = m;
    double d;
    int i, j, x;
    int opr = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);
    for (i = 0;i < m;i++) {
        for (j = 0;j < n;j++) {
            d = 0;
            switch (opr) {
                case 1:
                    for (x = 0;x < k;x++)
                        d += a[lda * x+ i] * b[ldb * j+ x];
                    break;
                case 2:
                    for (x = 0;x < k;x++)
                        d += a[lda * x+ i] * b[ldb * x+ j];
                    break;
                case 3:
                    for (x = 0;x < k;x++)
                        d += a[lda * i+ x] * b[ldb * j+ x];
                    break;
                case 4:
                    for (x = 0;x < k;x++)
                        d += a[lda * i+ x] * b[ldb * x+ j];
                    break;
            }
            c[ldc * j+ i] =alpha * d+beta * c[ldc * j + i];
        }
    }
}



void CMat_Multiply(double *left, double *right, double *res, int row, int mid, int col)
{
    int i, j, k;
    for (i = 0; i < row; i++)
        for (j = 0; j < col; j++)
        {
            res[i+j*row] = 0.0;
            for (k = 0; k < mid; k++)
            {
                if (left[i+k*row] != 0 && right[k+j*mid] != 0)
                    res[i+j*row] += left[i+k*row] * right[k+j*mid];
            }
        }
}


void CMat_PrintMatrix(double *mat, int row, int col, const char *premsg) {
    if (row < 0 || col < 0) {
        printf("***print_matrix wrong input arguments for the row and col ***\n");
        return;
    }
    int i, j;
    if (strlen(premsg) != 0) {
        printf("%s\n", premsg);
    }
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%13.7lf ", mat[i + j * row]);
        }
        printf("\n");
    }
}

void CMat_Transpose(double *A, int row, int col) {
    double *ACopy = mat(row, col);
    CMat_Copy(ACopy, A, col, row);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            A[j + i * col] = ACopy[i + j * row];
        }
    }
}


void CMat_Copy(double *des, double *src, int col, int row) {
    int len = col * row;
    for (int i = 0; i < len; ++i) {
        des[i] = src[i];
    }
}


double Mat_Norm(double *A, int len) {
    double res = 0.0;
    for (int i = 0; i < len; i++)
        res += A[i] * A[i];
    return sqrt(res);
}

void Mat_Add(double *A, double *B, int row, int col) {
    int len = row * col;
    for (int i = 0; i < len; i++) {
        A[i] += B[i];
    }
}

void Mat_Minus(double *A, double *B, int row, int col) {
    int len = row * col;
    for (int i = 0; i < len; i++) {
        A[i] -= B[i];
    }
}


