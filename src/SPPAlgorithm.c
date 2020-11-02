//
// Created by eureka on 2020/10/16.
//

#include "SPPAlgorithm.h"

/*
 * prepare data for adjusting
 */
void GetData(RNXOBS *OB, int *obsCountAva, Coord *satPos, Coord *groundRef, double t, double *obs_c1w, double *corr,
             double *elev) {
    for (int k = 0; k < MAXSAT; k++) {
        double c1w = OB->obs[k][3];
        double c2w = OB->obs[k][4];
        if (fabs(c1w) > 1e-3) {
            double toe = 99999999999999999999.0;
            double dt = fabs(t - toe);
            for (int i = 0; i < MAXEPH - 1; i++) {
                // 寻找最近参考时刻
                if (!strcmp(ephgps[0][i].cprn, OB->cprn[k])) {
                    if ((fabs(t - ephgps[0][i].toe) < dt)
                        && (ephgps[0][i + 1].cprn[2] == ephgps[0][i].cprn[2])) {  //or the last epoch
                        toe = ephgps[0][i].toe;
                        dt = fabs(t - toe);
                    } else //假定星历按时间顺序排列, 最佳时间差值为最小, 此前时间差值一直递减, 一旦不再出现递减, 则上次循环为最佳参考历元
                    {
                        satPos[*obsCountAva] = SatPos_CorrSpread(&ephgps[0][i - 1], groundRef, t);
                        double t_spread = Dist(&satPos[*obsCountAva], groundRef) / CLIGHT; //传播时间修正
                        Coord satVel = SatVel(&ephgps[0][i - 1], t - t_spread); //

                        double satClk = SatClk(&ephgps[0][i - 1], t - t_spread);
                        double earthRot_corr = SatPos_CorrEarthRot(&satPos[*obsCountAva], groundRef);
                        double relative_corr = SatPos_CorrRelative(&satPos[*obsCountAva], &satVel);
                        double azel[2];
                        elevazel(&satPos[*obsCountAva], groundRef, azel);

                        double pos[3];
                        ecef2pos(groundRef, pos);
                        double ionPara[] = {0.9313e-08, 0.1490e-07, -0.5960e-07, -0.1192e-06,
                                            0.8806e+05, 0.4915e+05, -0.1311e+06, -0.3277e+06};
                        int iyear, imonth, iday, ih, imin;
                        double sec;
                        mjd2date(OB->jd, OB->tsec, &iyear, &imonth, &iday, &ih, &imin, &sec);
                        double ep[] = {iyear, imonth, iday, ih, imin, sec};
                        gtime_t gtime = epoch2time(&ep);
                        double ion = ionmodel(gtime, ionPara, &pos, &azel);
                        double trp = tropmodel(gtime, &pos, &azel, 0.5);
                        double tgd = TGD(ephgps[0][i - 1].cprn);

                        obs_c1w[*obsCountAva] = c1w;
                        corr[*obsCountAva] = -satClk + trp + ion + earthRot_corr + relative_corr + tgd;
                        elev[*obsCountAva] = azel[1] * RAD2DEG;
                        (*obsCountAva)++;
                        break;
                    }
                }
            }
        }
    }
}


int Adjust(double *B, double *P, double *l, int countSat, double *x, double tol, double *sigma) {
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

    double *V = mat(countSat, 1);
    CMat_Matmul("NN", countSat, 1, 4, 1.0, B, dx, 0.0, V);
    Mat_Minus(V, l, countSat, 1);

//    CMat_PrintMatrix(V,countSat,1,"V");

    *sigma = Sigma(B, P, l, countSat, countSat, dx);
    status = (Mat_Norm(dx, 4) < tol);

    int satAva = countSat;
    if (status) {
        Mat_Minus(x, dx, 4, 1);
        for (int j = 0; j < 5; ++j) {
            satAva = countSat;
            Mat_Add(x, dx, 4, 1);
            for (int i = 0; i < countSat; ++i) {
                if ((*sigma > 2.0) & (fabs(V[i]) > 1 * (*sigma))) {
                    P[i + i * countSat] = 0;
                    satAva--;
                } else if (fabs(P[i + i * countSat] - 1) < 1e-8)
                    P[i + i * countSat] = 1;
            }

            CMat_Matmul("NN", 4, countSat, countSat, 1.0, BT, P, 0.0, BTP);
            CMat_Matmul("NN", 4, 4, countSat, 1.0, BTP, B, 0.0, N);
            CMat_Matmul("NN", 4, 1, countSat, 1.0, BTP, l, 0.0, W);
            CMat_Inverse(N, 4, 4);
            CMat_Matmul("NN", 4, 1, 4, 1.0, N, W, 0.0, dx);

            CMat_Matmul("NN", countSat, 1, 4, 1.0, B, dx, 0.0, V);
            Mat_Minus(V, l, countSat, 1);
            Mat_Add(x, dx, 4, 1);
            *sigma = Sigma(B, P, l, countSat, satAva, dx);
        }
        status = (Mat_Norm(dx, 4) < tol);
    }

    return status;
}

double Sigma(double *B, double *P, double *l, int countSat, int satAvalid, double *x) {
    //data quality
    //V=Bx-l
    //sigma^2=VT P V
    double *V = mat(countSat, 1);
    double *VT = mat(countSat, 1);
    double *VTP = mat(1, countSat);
    CMat_Matmul("NN", countSat, 1, 4, 1.0, B, x, 0.0, V);
    Mat_Minus(V, l, countSat, 1);
    CMat_Copy(VT, V, countSat, 1);
    CMat_Transpose(VT, countSat, 1);
    CMat_Matmul("NN", 1, countSat, countSat, 1.0, VT, P, 0.0, VTP);
    double sigma;
    for (int i = 0; i < countSat; ++i) {
        sigma += VTP[i] * V[i];
    }
    return sqrt(sigma / (satAvalid - 4));
}


void SPP(Coord *sitPos, Coord *satPos, double *obs, double *corr, int countSat, double *elev, double *sigma) {
    if (countSat < 5)
        printf("number of satllites is less than 4\n");
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

        if (elev[i] < 20) //截止高度角, 低于此降权处理
            P[i + i * countSat] = 0.05 * elev[i]; //线性处理，简化运算
        else {
            P[i + i * countSat] = 1.0;
        }

        P[i + i * countSat] = 1.0;
    }

    //Adjust and Update B,l
    int isOver = Adjust(B, P, l, countSat, x, 1.0e-6, sigma);
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
        isOver = Adjust(B, P, l, countSat, x, 1.0e-8, sigma);
        iterCount++;
        isOver = isOver | (iterCount > iterMax);
    }

    sitPos->X = x[0];
    sitPos->Y = x[1];
    sitPos->Z = x[2];

}