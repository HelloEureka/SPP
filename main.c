/*
 * main.cpp
 *
 *  Created on: 2020年9月20日
 *      Author: xps
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CMat.h"

#define GLOBAL_VAR

#include "rnx.h"
#include "global.h"
#include "SatPos.h"

void initialize() {
    int isat, isys;
    memset(&CKF, 0, sizeof(CKDCFG));
    strcpy(SYS, "GCRESJI");
    strcpy(OBSTYPE, "PWCIXSAQLDBYMZN ");
    // initialize the cprn
    for (isat = 0; isat < 32; isat++) {
        sprintf(CKF.cprn[isat], "G%02d", isat + 1);
    }
    CKF.nprn = 32;
    for (isys = 0; isys < MAXSYS; isys++) {
        switch (SYS[isys]) {
            case 'G':
                strcpy(CKF.freq[isys][0], "L1");
                strcpy(CKF.freq[isys][1], "L2");
                CKF.nfreq[isys] = 2;
                break;
            default:
                break;
        }
    }
    strcpy(CKF.cobs, "SF");
}

int main(int argc, char *args[]) {
    RNXHEAD HD = {0};
    RNXOBS OB = {0};
    BRDHEAD brdhd;
    initialize();

    CKF.mjd = (int) modified_julday(2020, 4, 22); /* the date should be set properly */
    CKF.sod = 30.0;
    const char *pfilobs = "../data/albh1130.20o";
    const char *pfilbrdm = "../data/brdm1130.20p";
    // const char *pfilbrdm = "../data/brdm2150.20p";
    FILE *fp = fopen(pfilobs, "r");
    if (!fp) {
        printf("can't open file to read!\n");
        exit(1);
    }
    read_rnxobs_head(&HD, fp);
    read_rnxobs(fp, CKF.mjd, CKF.sod, CKF.nprn, CKF.cprn, CKF.nfreq, CKF.freq, &HD, &OB);
    read_rnxnav('M', pfilbrdm, CKF.mjd, CKF.mjd + 1, &brdhd, neph, ephgps, ephgls);

    /*sat pos */
    Coord groundRef = {-2341333.112, -3539049.520, 4745791.248};
//    char gps_prn[32][4] = {"G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08",
//                           "G09", "G10", "G11", "G12", "G13", "G14", "G15", "G16",
//                           "G17", "G18", "G19", "G20", "G21", "G22", "G23", "G24",
//                           "G25", "G26", "G27", "G28", "G29", "G30", "G31", "G32"};
    int counter = 0, j, i; //加快循环, 使得循环一次 MAXEPH
//    for (j = 0; j < 32; j++) {
//        double t = 3 * 24 * 3600 + 30;
//        double toe = 99999999999999999999.0;
//        double dt = fabs(t - toe);
//        for (i = counter; i < MAXEPH; i++) {
//            // 寻找最近参考时刻
//            if (!strcmp(ephgps[0][i].cprn, gps_prn[j])) {
//                if (fabs(t - ephgps[0][i].toe) < dt) {
//                    toe = ephgps[0][i].toe;
//                    dt = fabs(t - toe);
//                } else //假定星历按时间顺序排列, 最佳时间差值为最小, 此前时间差值一直递减, 一旦不再出现递减, 则上次循环为最佳参考历元
//                {
//                    printf("%s  ", ephgps[0][i - 1].cprn);
////                    Coord satPos = SatPos_Cal(&ephgps[0][i - 1], t);
//                    Coord satPos = SatPos_CorrSpread(&ephgps[0][i - 1], &groundRef, 3 * 24 * 3600 + 30); //传播时间修正
//                    Coord satVel = SatVel(&ephgps[0][i - 1], &groundRef, 3 * 24 * 3600 + 30);
//                    counter = i;
//
//                    double earthRot_corr = SatPos_CorrEarthRot(&satPos, &groundRef);
//                    double relative_corr = SatPos_CorrRelative(&satPos, &satVel);
//                    printf("X = %14.3f  Y = %14.3f  Z = %14.3f  ", satPos.X, satPos.Y, satPos.Z);
//                    printf("VX = %10.3f  VY = %10.3f  VZ = %10.3f  ", satVel.X, satVel.Y, satVel.Z);
//                    printf("send time = %10.3f   earthRot = %8.3f  relative = %8.3f\n",
//                           fmod(t, 24 * 3600) - Dist(&satPos, &groundRef) / CLIGHT, earthRot_corr, relative_corr);
//                    break;
//                    // Coord satPos = SatPos_Cal_Corr(&ephgps[0][i], &groundRef, 3 * 24 * 3600 + 30); //传播时间修正
//
//                }
//            }
//        }
//    }




    //wite status.log
//    fp = fopen("../data/status.log", "w");
//    fprintf(fp, "PRN      X           Y            Z          TIME_SEND        OBS_L1       OBS_L2        S_X           S_Y          S_Z           S_VZ        S_VY       S_VZ        SATCLK      TRP    ION    EARTH_ROT   RELA    TGD     ELEV  \n");
//
//    double t = 3 * 24 * 3600 + 30;
//    for(int k = 0;k<MAXSAT;k++){
//        double c1w = OB.obs[k][3];
//        double c2w = OB.obs[k][4];
//        if(fabs(c1w)>1e-3)
//        {
//            fprintf(fp,"%s ", OB.cprn[k]);
//            fprintf(fp,"%12.3f  %12.3f  %12.3f ", groundRef.X,groundRef.Y,groundRef.Z);
//
//
//            double toe = 99999999999999999999.0;
//            double dt = fabs(t - toe);
//            for (int i = counter; i < MAXEPH; i++) {
//                // 寻找最近参考时刻
//                if (!strcmp(ephgps[0][i].cprn, OB.cprn[k])) {
//                    if (fabs(t - ephgps[0][i].toe) < dt) {
//                        toe = ephgps[0][i].toe;
//                        dt = fabs(t - toe);
//                    } else //假定星历按时间顺序排列, 最佳时间差值为最小, 此前时间差值一直递减, 一旦不再出现递减, 则上次循环为最佳参考历元
//                    {
//                        Coord satPos = SatPos_CorrSpread(&ephgps[0][i - 1], &groundRef, 3 * 24 * 3600 + 30); //传播时间修正
//                        Coord satVel = SatVel(&ephgps[0][i - 1], &groundRef, 3 * 24 * 3600 + 30);
//                        counter = i;
//
//                        double earthRot_corr = SatPos_CorrEarthRot(&satPos, &groundRef);
//                        double relative_corr = SatPos_CorrRelative(&satPos, &satVel);
//
//                        fprintf(fp,"%d %6.3f ", 58961, fmod(t, 24 * 3600) - Dist(&satPos, &groundRef) / CLIGHT);
//                        fprintf(fp, "%12.3f  %12.3f ", c1w, c2w);
//                        fprintf(fp, "%13.3f  %13.3f  %13.3f ", satPos.X, satPos.Y, satPos.Z);
//                        fprintf(fp,"%9.3f  %9.3f  %9.3f ", satVel.X, satVel.Y, satVel.Z);
//                        fprintf(fp,"   %6.2f     %6.2f  %6.2f  %6.2f      %6.2f  %6.2f %6.2f\n",
//                               0.0, 0.0, 0.0, earthRot_corr, relative_corr,0.0,0.0, 0.0);
//                        break;
//                    }
//                }
//            }
//        }
//    }
//    fclose(fp);



    demo();

    printf("process done!\n");
}

