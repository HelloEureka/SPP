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


void GetData(RNXOBS *OB, int *obsCountAva, Coord *satPos, Coord *groundRef, double t, double *obs_c1w, double *corr,
             double *elev) {
    int counter = 0; //加快循环, 使得循环一次 MAXEPH
    for (int k = 0; k < MAXSAT; k++) {
        double c1w = OB->obs[k][3];
        double c2w = OB->obs[k][4];
        if (fabs(c1w) > 1e-3) {
            double toe = 99999999999999999999.0;
            double dt = fabs(t - toe);
            for (int i = counter; i < MAXEPH; i++) {
                // 寻找最近参考时刻
                if (!strcmp(ephgps[0][i].cprn, OB->cprn[k])) {
                    if (fabs(t - ephgps[0][i].toe) < dt) {
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
                        double tgd = 0;

                        obs_c1w[*obsCountAva] = c1w;
                        corr[*obsCountAva] = -satClk + trp + ion + earthRot_corr + relative_corr + tgd;
                        elev[*obsCountAva] = azel[1] * RAD2DEG;
                        (*obsCountAva)++;
                        counter = i;
                        break;
                    }
                }
            }
        }
    }
}


int main(int argc, char *args[]) {
    RNXHEAD HD = {0};
    RNXOBS OB = {0};
    BRDHEAD brdhd;
    initialize();

    CKF.mjd = (int) modified_julday(2020, 4, 22); /* the date should be set properly */


    const char *pfilobs = "../data/albh1130.20o";
    const char *pfilbrdm = "../data/brdm1130.20p";
    FILE *fp = fopen(pfilobs, "r");
    if (!fp) {
        printf("can't open file to read!\n");
        exit(1);
    }
    read_rnxobs_head(&HD, fp);
    read_rnxnav('M', pfilbrdm, CKF.mjd, CKF.mjd + 1, &brdhd, neph, ephgps, ephgls);

    const int OBSCOUNTMAX = 40;  //
    Coord sitPos = {0, 0, 0};
    Coord groundRef = {-2341333.112, -3539049.520, 4745791.248};
    Coord satPos[OBSCOUNTMAX];  //
    double corr[OBSCOUNTMAX];
    double obs_c1w[OBSCOUNTMAX];
    double elev[OBSCOUNTMAX];
    int week;
    double sow;

    FILE *fpPos = fopen("../data/ecef.txt", "w");
    for (CKF.sod = 0; CKF.sod <3600*23 ; CKF.sod += 1.0) {
        read_rnxobs(fp, CKF.mjd, CKF.sod, CKF.nprn, CKF.cprn, CKF.nfreq, CKF.freq, &HD, &OB);
        if(fabs(OB.tsec - CKF.sod)>0.1)
            continue;
        int obsCountAva = 0;

        mjd2wksow(CKF.mjd, CKF.sod, &week, &sow);
//        GetData(&OB, &obsCountAva, &satPos, &groundRef, sow, &obs_c1w, &corr, &elev);
        GetData(&OB, &obsCountAva, &satPos, &sitPos, sow, &obs_c1w, &corr, &elev);
        SPP(&sitPos, &satPos, &obs_c1w, &corr, obsCountAva, &elev);

        fprintf(fpPos, "%6.0f  %10.3f  %10.3f  %10.3f\n", CKF.sod, sitPos.X, sitPos.Y, sitPos.Z);
    }
    fclose(fpPos);
    fclose(fp);

    printf("process done!\n");
}

