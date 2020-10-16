/*
 * main.cpp
 *
 *  Created on: 2020年9月20日
 *      Author: xps
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GLOBAL_VAR

#include "rnx.h"
#include "global.h"

#include "SPPAlgorithm.h"

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

// //    wite status.log
    fp = fopen("../data/status.log", "w");
    fprintf(fp,
            "PRN      X           Y            Z         "
            " TIME_SEND        OBS_L1       OBS_L2        S_X         "
            "  S_Y          S_Z           S_VZ        S_VY       S_VZ     "
            "   SATCLK           TRP     ION   EARTH_ROT   RELA    TGD     ELEV  \n");

    int obsCountAva = 0;
    const int OBSCOUNTMAX = 40;  //
    Coord satPos[OBSCOUNTMAX];  //
    double corr[OBSCOUNTMAX];
    Coord sitPos = {0, 0, 0};
    double obs_c1w[OBSCOUNTMAX];
    double elev[OBSCOUNTMAX];

    double t = 3 * 24 * 3600 + 30;
    for (int k = 0; k < MAXSAT; k++) {
        double c1w = OB.obs[k][3];
        double c2w = OB.obs[k][4];
        if (fabs(c1w) > 1e-3) {
            double toe = 99999999999999999999.0;
            double dt = fabs(t - toe);
            for (int i = counter; i < MAXEPH; i++) {
                // 寻找最近参考时刻
                if (!strcmp(ephgps[0][i].cprn, OB.cprn[k])) {
                    if (fabs(t - ephgps[0][i].toe) < dt) {
                        toe = ephgps[0][i].toe;
                        dt = fabs(t - toe);
                    } else //假定星历按时间顺序排列, 最佳时间差值为最小, 此前时间差值一直递减, 一旦不再出现递减, 则上次循环为最佳参考历元
                    {
                        satPos[obsCountAva] = SatPos_CorrSpread(&ephgps[0][i - 1], &groundRef,
                                                                3 * 24 * 3600 + 30); //传播时间修正
                        double t_spread = Dist(&satPos[obsCountAva], &groundRef) / CLIGHT;
                        Coord satVel = SatVel(&ephgps[0][i - 1], t - t_spread);
                        counter = i;

                        double earthRot_corr = SatPos_CorrEarthRot(&satPos[obsCountAva], &groundRef);
                        double relative_corr = SatPos_CorrRelative(&satPos[obsCountAva], &satVel);
                        double azel[2];
                        elevazel(&satPos[obsCountAva], &groundRef, azel);

                        double satClk = SatClk(&ephgps[0][i - 1], t - t_spread);


                        double pos[3];
                        ecef2pos(&groundRef, pos);
                        double ionPara[] = {0.9313e-08, 0.1490e-07, -0.5960e-07, -0.1192e-06,
                                            0.8806e+05, 0.4915e+05, -0.1311e+06, -0.3277e+06};

                        int iyear, imonth, iday, ih, imin;
                        double sec;
                        mjd2date(OB.jd, OB.tsec, &iyear, &imonth, &iday, &ih, &imin, &sec);
                        double ep[]={iyear,imonth,iday,ih,imin,sec};
                        gtime_t gtime = epoch2time(&ep);
                        double ion = ionmodel(gtime, ionPara, &pos, &azel);
                        double trp = tropmodel(gtime, &pos, &azel, 0.5);
                        double tgd = 0;
                        printf("azim:%9.3lf elev %9.3lf\n", azel[0] * RAD2DEG, azel[1] * RAD2DEG);
                        fprintf(fp, "%s ", OB.cprn[k]);
                        fprintf(fp, "%12.3f  %12.3f  %12.3f ", groundRef.X, groundRef.Y, groundRef.Z);
                        fprintf(fp, "%d %6.3f ", 58961, fmod(t, 24 * 3600) - t_spread);
                        fprintf(fp, "%12.3f  %12.3f ", c1w, c2w);
                        fprintf(fp, "%13.3f  %13.3f  %13.3f ",
                                satPos[obsCountAva].X, satPos[obsCountAva].Y, satPos[obsCountAva].Z);
                        fprintf(fp, "%9.3f  %9.3f  %9.3f ", satVel.X, satVel.Y, satVel.Z);
                        fprintf(fp, "   %11.3f     %6.2f  %6.2f  %6.2f      %6.2f  %6.2f %5.2f\n",
                                satClk, trp, ion, earthRot_corr, relative_corr, tgd, azel[1] * RAD2DEG);

                        obs_c1w[obsCountAva] = c1w;
                        corr[obsCountAva] = -satClk + trp + ion + earthRot_corr + relative_corr + tgd;
                        elev[obsCountAva] = azel[1] * RAD2DEG;
                        obsCountAva++;

                        break;
                    }
                }
            }
        }
    }
    fclose(fp);
    SPP(&sitPos, &satPos, &obs_c1w, &corr, obsCountAva, &elev);

    printf("process done!\n");
}

