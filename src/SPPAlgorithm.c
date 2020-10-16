//
// Created by eureka on 2020/10/16.
//

#include "SPPAlgorithm.h"


int Adjust(double *B, double *P, double *l, int countSat, double *x, double tol) {
    double *W = mat(4, 1); //BTPl
    double *N = mat(4, 4);
    double *BT = mat(countSat, 4);
    double *BTP = mat(4, countSat);
    double *dx = mat(4, 1);
    int status;
    CMat_Copy(BT, B, countSat, 4);
    CMat_Transpose(BT, countSat, 4);
    CMat_Matmul("NN", 4, countSat, countSat, 1.0, BT, P, 0.0, BTP);
    CMat_Matmul("NN", 4, 4, countSat, 1.0, BTP, B, 0.0, N);
    CMat_Matmul("NN", 4, 1, countSat, 1.0, BTP, l, 0.0, W);
    CMat_Inverse(N, 4, 4);
    CMat_Matmul("NN", 4, 1, 4, 1.0, N, W, 0.0, dx);

    Mat_Add(x, dx, 4, 1);
    status = (Mat_Norm(dx, 4) < tol);
    return status;
}

void SPP(Coord *sitPos, Coord *satPos, double *obs, double *corr, int countSat, double *elev) {
    double *B = mat(countSat, 4);
    double *P = mat(countSat, countSat);
    double *l = mat(countSat, 1);
    double dZ, dY, dX, rho_0;
    int i;
    //initialize
    double x[] = {sitPos->X, sitPos->Y, sitPos->Z, 0.0};
    for (i = 0; i < countSat; i++) {
        dX = x[0] - satPos[i].X;
        dY = x[1] - satPos[i].Y;
        dZ = x[2] - satPos[i].Z;
        rho_0 = sqrt(dX * dX + dY * dY + dZ * dZ);

        B[i + 0 * countSat] = dX / rho_0;
        B[i + 1 * countSat] = dY / rho_0;
        B[i + 2 * countSat] = dZ / rho_0;
        B[i + 3 * countSat] = 1.0;
        l[i] = obs[i] - (rho_0 + corr[i]) - x[3];

        if (elev[i] < 30) //截止高度角, 低于此降权处理
            P[i + i * countSat] = 0.05 * (elev[i] - 10); //线性处理，简化运算
        else {
            P[i + i * countSat] = 1;
        }
    }

    //Adjust and Update B,l
    int isOver = Adjust(B, P, l, countSat, x, 1.0e-6);
    int iterCount = 0;
    int iterMax = 10;
    while (!isOver) {
        for (i = 0; i < countSat; i++) {
            dX = x[0] - satPos[i].X;
            dY = x[1] - satPos[i].Y;
            dZ = x[2] - satPos[i].Z;
            rho_0 = sqrt(dX * dX + dY * dY + dZ * dZ);

            B[i + 0 * countSat] = dX / rho_0;
            B[i + 1 * countSat] = dY / rho_0;
            B[i + 2 * countSat] = dZ / rho_0;
            l[i] = obs[i] - (rho_0 + corr[i]) - x[3];
        }
        isOver = Adjust(B, P, l, countSat, x, 1.0e-8);
        iterCount++;
        isOver=isOver | (iterCount>iterMax);
    }
    sitPos->X = x[0];
    sitPos->Y = x[1];
    sitPos->Z = x[2];
}