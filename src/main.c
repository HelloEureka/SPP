#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GLOBAL_VAR

#include "rnx.h"
#include "global.h"

#include "SPPAlgorithm.h"

int main(int argc, char *args[]) {
    //initialize for renix
    RNXHEAD HD = {0};
    RNXOBS OB = {0};
    BRDHEAD brdhd;
    initialize();
    CKF.mjd = (int) modified_julday(2020, 9, 17); /* the date should be set properly */
//    const char *pfilobs = "../data/hkkt/hkkt1130.20o";
//    const char *pfilbrdm = "../data/hkkt/brdm1130.20p";

    //read renix
//    const char *pfilobs = "../data/albh/albh1130.20o";
//    const char *pfilbrdm = "../data/albh/brdm1130.20p";

    const char *pfilobs = "../data/decoder/obs.out";
    const char *pfilbrdm = "../data/decoder/nav.out";

    FILE *fp = fopen(pfilobs, "r");
    if (!fp) {
        printf("can't open file to read!\n");
        exit(1);
    }
    read_rnxobs_head(&HD, fp);
    read_rnxnav('M', pfilbrdm, CKF.mjd, CKF.mjd + 1, &brdhd, neph, ephgps, ephgls);
    SortEph();


    // SPP
    const int OBSCOUNTMAX = 40;  //
    //    Coord sitPos_ref = {-2405145.476, 5385196.812,2420034.840};
    Coord sitPos_ref = {-2341333.112,  -3539049.520,4745791.248};
    Coord sitPos = {0, 0, 0};
    Coord satPos[OBSCOUNTMAX];  //
    double corr[OBSCOUNTMAX];
    double obs_c1w[OBSCOUNTMAX];
    double elev[OBSCOUNTMAX];
    int week;
    double sow;
    double sigma;

    FILE *fpPos = fopen("../data/decoder/res_spp.txt", "w");
    for (CKF.sod = 0; CKF.sod < 3600 * 24; CKF.sod += 30.0) {
//    for (CKF.sod = 720; CKF.sod <= 720; CKF.sod += 30.0) {
        read_rnxobs(fp, CKF.mjd, CKF.sod, CKF.nprn, CKF.cprn, CKF.nfreq, CKF.freq, &HD, &OB);
        if (fabs(OB.tsec - CKF.sod) > 0.1)
            continue;
        int obsCountAva = 0;

        mjd2wksow(CKF.mjd, CKF.sod, &week, &sow);
        GetData(&OB, &obsCountAva, &satPos, &sitPos_ref, sow, &obs_c1w, &corr, &elev);

        SPP(&sitPos, &satPos, &obs_c1w, &corr, obsCountAva, &elev, &sigma);
        double errX = sitPos_ref.X - sitPos.X;
        double errY = sitPos_ref.Y - sitPos.Y;
        double errZ = sitPos_ref.Z - sitPos.Z;
        fprintf(fpPos, "%8.0f  %10.3f  %10.3f  %10.3f  %8.3f %8.3f %8.3f %6.3f  %d\n", CKF.sod, sitPos.X, sitPos.Y,
                sitPos.Z,
                errX, errY, errZ, sigma, obsCountAva);
    }
    fclose(fpPos);
    fclose(fp);

    printf("process done!\n");
}

