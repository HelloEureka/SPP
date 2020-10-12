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
    return p;
}

int *imat(int n, int m) {
    int *p;
    if (n <= 0 || m <= 0)
        return NULL;
    if (!(p = (int *) malloc(sizeof(int) * n * m))) {
        printf("***ERROR(mat):allocate error!\n");
        return NULL;
    }
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
    return 0;
}

void CMat_Matmul(const char *tr, int m, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta,
                 double *c, int ldc) {
    double d;
    int i, j, x;
    int opr = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            d = 0;
            switch (opr) {
                case 1:
                    for (x = 0; x < k; x++)
                        d += a[lda * x + i] * b[ldb * j + x];
                    break;
                case 2:
                    for (x = 0; x < k; x++)
                        d += a[lda * x + i] * b[ldb * x + j];
                    break;
                case 3:
                    for (x = 0; x < k; x++)
                        d += a[lda * i + x] * b[ldb * j + x];
                    break;
                case 4:
                    for (x = 0; x < k; x++)
                        d += a[lda * i + x] * b[ldb * x + j];
                    break;
            }
            c[ldc * j + i] = alpha * d + beta * c[ldc * j + i];
        }
    }
}

void CMat_PrintMatrix(double *mat, int lda, int row, int col, const char *premsg) {
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
            printf("%13.7lf ", mat[i + j * lda]);
        }
        printf("\n");
    }
}


void demo() {
    // matrix is column major
    double norm[3 * 3] = {1, 2, 3, 6, 3, 5, 4, 8, 9};
    double A[3] = {4, 5, 6}, L[3];
    CMat_Matmul("NN", 3, 1, 3, 1.0, norm, 3, A, 3, 0.0, L, 3);
    printf("norm(3*3) * A(3*1) = L(3*1)\n")  ;
    CMat_PrintMatrix(norm, 3, 3, 3, "norm");
    CMat_PrintMatrix(A, 3, 3, 1, "A");
    CMat_PrintMatrix(L, 3, 3, 1, "L");

    int status = CMat_Inverse(norm, 3, 3);
    if (status == -1) {
        printf("Singular matrix of norm\n") ;
        exit(1);
    }
    printf("inverse of norm(3*3) is :\n");
    CMat_PrintMatrix(norm, 3, 3, 3, "norm");
}
