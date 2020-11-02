//
// Created by eureka on 2020/10/12.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rnx.h"
#include "global.h"

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


void read_rnxnav(char csys, const char *flnbrd, double mjd0, double mjd1, BRDHEAD *hd,
                 int *neph, GPS_BRDEPH ephm[MAXSYS][MAXEPH], GLONASS_BRDEPH ephg[MAXEPH])
{
    FILE *fp;
    int i, k, isys, iy, im, id, ih, imin, lastpos, len;
    char line[LEN_STRING], buf[LEN_STRING];
    char *ptr;
    double dt0, dt1, isec;
    double data1, data2, data3;
    int already;

    GPS_BRDEPH eph = {0};
    GLONASS_BRDEPH ephr = {0};

    if ((fp = fopen(flnbrd, "r")) == NULL)
    {
        printf("***ERROR(read_rnxnav):cant open brd file %s\n", flnbrd);
        exit(1);
    }
    memset(neph, 0, sizeof(int) * MAXSYS);
    strncpy(line, "\0", LEN_STRING);

    while (strstr(line, "END OF HEADER") == NULL)
    {
        strncpy(line, "\0", LEN_STRING);
        fgets(line, LEN_STRING, fp);
        if (strstr(line, "RINEX VERSION / TYPE") != NULL)
        {
            sscanf(line, "%lf", &hd->ver);
            if (hd->ver > 3.0)
            {
                if (csys == 'M' && line[40] != 'M')
                {
                    printf(
                            "###WARNING(read_rnxnav):there only broadcast for line[40]!\n");
                }
            }
        }
        else if (strstr(line, "PGM / RUN BY / DATE") != NULL)
            continue;
        else if (strstr(line, "COMMENT") != NULL)
            continue;
        else if (strstr(line, "ION ALPHA") != NULL)
        {
        }
        else if (strstr(line, "ION BETA") != NULL)
        {
        }
            //rinex 3.00 3.01 3.02
        else if (strstr(line, "IONOSPHERIC CORR") != NULL)
        {
            if (!strncmp(line, "GPS ", 4) || !strncmp(line, "GPSA", 4))
            {
                i = index_string(SYS, 'G');
                k = 0;
            }
            else if (!strncmp(line, "GPSB", 4))
            {
                i = index_string(SYS, 'G');
                k = 1;
            }
            else if (!strncmp(line, "GAL ", 4))
            {
                i = index_string(SYS, 'E');
                k = 0;
            }
            else if (!strncmp(line, "BDS ", 4) || !strncmp(line, "BDSA", 4))
            {
                i = index_string(SYS, 'C');
                k = 0;
            }
            else if (!strncmp(line, "BDSB", 4))
            {
                i = index_string(SYS, 'C');
                k = 1;
            }
            else if (!strncmp(line, "QZS ", 4) || !strncmp(line, "QZSA", 4))
            {
                i = index_string(SYS, 'J');
                k = 0;
            }
            else if (!strncmp(line, "QZSB", 4))
            {
                i = index_string(SYS, 'J');
                k = 1;
            }
            else
                printf("###WARNING(read_rnxnav):unknown ionospheric corr!\n");

            strncpy(hd->ionc[i][k], line, 4);
            ptr = line;
            while (*ptr != '\0')
            {
                if (*ptr == 'D')
                    *ptr = 'e';
                ptr++;
            }
            sscanf(line, "%*s%lf%lf%lf%lf", hd->ion[i][k], hd->ion[i][k] + 1,
                   hd->ion[i][k] + 2, hd->ion[i][k] + 3);
        }
        else if (strstr(line, "DELTA-UTC: A0,A1,T,W") != NULL)
        {
            i = index_string(SYS, csys);
            if (i == -1)
            {
                printf(
                        "$$$MESSAGE(read_rnxnav):DELTA-UTC: A0,A1,T,W is only valid for single system!\n");
                continue;
            }
            ptr = line;
            while (*ptr != '\0')
            {
                if (*ptr == 'D')
                    *ptr = 'e';
                ptr++;
            }
            sscanf(line, "%lf%lf%lf%lf", hd->tim[i][0], hd->tim[i][0] + 1,
                   hd->tim[i][0] + 2, hd->tim[i][0] + 3);
        }
        else if (strstr(line, "TIME SYSTEM CORR") != NULL)
        {
            if (!strncmp(line, "GPUT", 4))
            {
                i = index_string(SYS, 'G');
                k = 0;
            }
            else if (!strncmp(line, "GLUT", 4))
            {
                i = index_string(SYS, 'R');
                k = 0;
            }
            else if (!strncmp(line, "GAUT", 4))
            {
                i = index_string(SYS, 'E');
                k = 0;
            }
            else if (!strncmp(line, "BDUT", 4))
            {
                i = index_string(SYS, 'C');
                k = 0;
            }
            else if (!strncmp(line, "SBUT", 4))
            {
                i = index_string(SYS, 'S');
                k = 0;
            }
            else if (!strncmp(line, "QZUT", 4))
            {
                i = index_string(SYS, 'J');
                k = 0;
            }
            else if (!strncmp(line, "GAGP", 4))
            {
                i = index_string(SYS, 'G');
                k = 1;
            }
            else if (!strncmp(line, "GLGP", 4))
            {
                i = index_string(SYS, 'R');
                k = 1;
            }
            else if (!strncmp(line, "QZGP", 4))
            {
                i = index_string(SYS, 'J');
                k = 1;
            }
            else
            {
                printf(
                        "$$$MESSAGE(read_rnxnav):unknown unknown TIME SYSTEM CORR!\n");
                continue;
            }
            strncpy(hd->timc[i][k], line, 4);
            ptr = line;
            while (*ptr != '\0')
            {
                if (*ptr == 'D')
                    *ptr = 'e';
                ptr++;
            }
            sscanf(line, "%*5c%lf%lf%lf%lf", hd->tim[i][k], hd->tim[i][k] + 1,
                   hd->tim[i][k] + 2, hd->tim[i][k] + 3);
        }
        else if (strstr(line, "LEAP SECONDS") != NULL)
            sscanf(line, "%d", &hd->leap);
    }
    while (!feof(fp))
    {
        lastpos = ftell(fp);
        fgets(line, LEN_STRING, fp);
        fseek(fp, lastpos - ftell(fp), SEEK_CUR);
        memset(&eph, 0, sizeof(eph));
        if (hd->ver < 3.0)
        {
            strncpy(buf + 1, line, strlen(line));
            buf[0] = csys;
        }
        else
            buf[0] = line[0];

        switch (buf[0])
        {
            case 'G':
            case 'E':
            case 'C':
            case 'J':
                isys = index_string(SYS, buf[0]);
                len = 0;
                strncpy(buf, "", LEN_STRING);
                for (i = 0; i < 8; i++)
                {
                    fgets(line, LEN_STRING, fp);
                    if (i != 7)
                        filleph(line, hd->ver);
                    if (hd->ver < 3)
                    {
                        if (buf[0] != csys)
                            buf[0] = csys;
                        strncpy(buf + 1 + len, line, strlen(line));
                    }
                    else
                        strncpy(buf + len, line, strlen(line));
                    len += strlen(line);
                }
                ptr = buf;

                while (*ptr != '\0')
                {
                    if (*ptr == '\n' || *ptr == '\r')
                        *ptr = ' ';

                    if (*ptr == 'D')
                        *ptr = 'e';
                    ptr++;
                }
                if (hd->ver < 3.0)
                {
                    buf[0] = csys;
                    if (buf[1] == ' ')
                        buf[1] = '0';
                }
                if (buf[0] == 'C')
                {
                    sscanf(buf,
                           "%s%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                           eph.cprn, &iy, &im, &id, &ih, &imin, &isec, &eph.a0,
                           &eph.a1, &eph.a2, &eph.aode, &eph.crs, &eph.dn, &eph.m0,
                           &eph.cuc, &eph.e, &eph.cus, &eph.roota, &eph.toe,
                           &eph.cic, &eph.Omega0, &eph.cis, &eph.i0, &eph.crc,
                           &eph.omega, &eph.Omega_dot, &eph.i_dot, &eph.resvd0,
                           &eph.week, &eph.resvd1, &eph.accu, &eph.hlth, &eph.tgd,
                           &eph.aodc, &data1, &eph.iodc, &data3, &eph.signal_idx);
                    eph.signal_idx = (int)eph.signal_idx == 0 ? 1 : (int)eph.signal_idx;
                    if (eph.signal_idx > 6 && eph.signal_idx <= 13)
                    {
                        eph.delta_A = eph.roota;
                        eph.A_DOT = eph.resvd0;
                        eph.delta_n_dot = eph.resvd1;
                    }
                    if (eph.signal_idx == 1 || eph.signal_idx == 0)
                    { //B1I
                        eph.tgd_BDS[0] = eph.tgd;
                        eph.isc_BDS[0] = 0.0;
                    }
                    else if (eph.signal_idx == 2)
                    { //B2I
                        eph.tgd_BDS[1] = eph.aodc;
                        eph.isc_BDS[1] = 0.0;
                    }
                    else if (eph.signal_idx == 3)
                    { //B3I
                        eph.tgd_BDS[2] = 0.0;
                        eph.isc_BDS[2] = 0.0;
                    }
                    else if (eph.signal_idx == 4)
                    { //B1Q
                        eph.tgd_BDS[3] = eph.tgd;
                        eph.isc_BDS[3] = 0.0;
                    }
                    else if (eph.signal_idx == 5)
                    { //B1Q
                        eph.tgd_BDS[4] = eph.aodc;
                        eph.isc_BDS[4] = 0.0;
                    }
                    else if (eph.signal_idx == 6)
                    { //B3Q
                        eph.tgd_BDS[5] = 0.0;
                        eph.isc_BDS[5] = 0.0;
                    }
                    else if (eph.signal_idx == 7)
                    { //B1C
                        eph.tgd_BDS[6] = eph.tgd;
                        eph.isc_BDS[6] = data3;
                    }
                    else if (eph.signal_idx == 8)
                    { //B2a
                        eph.tgd_BDS[7] = eph.tgd;
                        eph.isc_BDS[7] = data3;
                    }
                    else if (eph.signal_idx == 9)
                    { //B2bI
                        eph.tgd_BDS[8] = eph.tgd;
                        eph.isc_BDS[8] = 0.0;
                    }
                    else if (eph.signal_idx == 10)
                    { //B2bQ
                        eph.tgd_BDS[9] = eph.tgd;
                        eph.isc_BDS[9] = 0.0;
                    }
                    else if (eph.signal_idx == 11)
                    { //B1A
                        eph.tgd_BDS[10] = eph.tgd;
                        eph.isc_BDS[10] = data3;
                    }
                    else if (eph.signal_idx == 12)
                    { //B3A
                        eph.tgd_BDS[11] = 0.0;
                        eph.isc_BDS[11] = data3;
                    }
                    else if (eph.signal_idx == 13)
                    { //B3AE
                        eph.tgd_BDS[12] = eph.aodc;
                        eph.isc_BDS[12] = data3;
                    }
                }
                else
                {
                    sscanf(buf,
                           "%s%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                           eph.cprn, &iy, &im, &id, &ih, &imin, &isec, &eph.a0,
                           &eph.a1, &eph.a2, &eph.aode, &eph.crs, &eph.dn, &eph.m0,
                           &eph.cuc, &eph.e, &eph.cus, &eph.roota, &eph.toe, &eph.cic,
                           &eph.Omega0, &eph.cis, &eph.i0, &eph.crc, &eph.omega,
                           &eph.Omega_dot, &eph.i_dot, &eph.resvd0, &eph.week,
                           &eph.resvd1, &eph.accu, &eph.hlth, &eph.tgd, &eph.aodc);
                }
                if (eph.hlth > 0)
                    continue;
                yr2year(&iy);
                if (iy < 2000)
                    continue;

                // time of clock
                eph.mjd = modified_julday(iy, im, id); // TIME IN BDS TIME AND SHOULD NOT CHANGE IT INTO GPST BECAUSE THERE ARE OTHER TIME TAG IN THE EPHEMERIS
                eph.sod = ih * 3600.0 + imin * 60.0 + isec;



                // adapter to bds
                if (eph.cprn[0] == 'C')
                {
                    //the week in broadcast file generated by WHU is GPS week,
                    //but that in IGS meraged file is BDS week
                    mjd2wksow(eph.mjd, eph.sod, &k, &dt1);
                    if (k != (int)eph.week)
                        eph.week = 1356 + eph.week;
                    eph.tgd1 = eph.aodc;
                }
                //check time
                dt0 = 0.0;
                dt1 = 0.0;
                if (mjd0 != 0)
                {
                    dt0 = eph.mjd + eph.sod / 86400.0 - mjd0;
                }

                if (mjd1 != 0)
                {
                    dt1 = eph.mjd + eph.sod / 86400.0 - mjd1;
                }
                if (dt0 < -1.0 / 24.0 || dt1 > 1.0 / 24.0)
                    continue;

                already = false;
                for (i = 0; i < neph[isys]; i++)
                {
                    if (strstr(ephm[isys][i].cprn, eph.cprn) && ephm[isys][i].mjd == eph.mjd && ephm[isys][i].sod == eph.sod)
                    {
                        already = true;
                        break;
                    }
                }
                if (!already)
                {
                    neph[isys] = neph[isys] + 1;
                    if (neph[isys] > MAXEPH)
                    {
                        printf(
                                "ERROR(read_rnxnav):exceed the maxium ephemeris number!\n");
                        exit(1);
                    }
                    memcpy(&ephm[isys][neph[isys] - 1], &eph, sizeof(eph));
                }
                break;
            case 'R':
                isys = index_string(SYS, buf[0]);
                len = 0;
                strncpy(buf, "", LEN_STRING);
                for (i = 0; i < 4; i++)
                {
                    fgets(line, LEN_STRING, fp);
                    if (hd->ver < 3)
                    {
                        if (buf[0] != csys)
                            buf[0] = csys;
                        strncpy(buf + 1 + len, line, strlen(line));
                    }
                    else
                        strncpy(buf + len, line, strlen(line));
                    len += strlen(line);
                }
                ptr = buf;

                while (*ptr != '\0')
                {
                    if (*ptr == '\n' || *ptr == '\r')
                    {
                        *ptr = ' ';
                    }
                    if (*ptr == 'D')
                    {
                        *ptr = 'e';
                    }
                    ptr++;
                }

                if (hd->ver < 3.0)
                {
                    buf[0] = csys;
                    if (buf[1] == ' ')
                        buf[1] = '0';
                }

                sscanf(buf,
                       "%s%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                       ephr.cprn, &iy, &im, &id, &ih, &imin, &isec, &ephr.tau,
                       &ephr.gamma, &ephr.tk, &ephr.pos[0], &ephr.vel[0],
                       &ephr.acc[0], &ephr.health, &ephr.pos[1], &ephr.vel[1],
                       &ephr.acc[1], &ephr.frenum, &ephr.pos[2], &ephr.vel[2],
                       &ephr.acc[2], &ephr.age);

                if (ephr.health > 0.0)
                    continue;

                yr2year(&iy);
                if (iy < 2000)
                    continue;

                ephr.mjd = modified_julday(iy, im, id);
                ephr.sod = ih * 3600.0 + imin * 60.0 + isec;

                //check time
                dt0 = 0.0;
                dt1 = 0.0;
                if (mjd0 != 0)
                {
                    dt0 = ephr.mjd + ephr.sod / 86400.0 - mjd0;
                }

                if (mjd1 != 0)
                {
                    dt1 = ephr.mjd + ephr.sod / 86400.0 - mjd1;
                }

                if (dt0 < -1.0 / 24.0 || dt1 > 1.0 / 24.0)
                    continue;

                already = false;
                for (i = 0; i < neph[isys]; i++)
                {
                    if (strstr(ephg[i].cprn, ephr.cprn) && ephg[i].mjd == ephr.mjd && ephg[i].sod == ephr.sod)
                    {
                        already = true;
                        break;
                    }
                }
                if (!already)
                {
                    neph[isys] = neph[isys] + 1;
                    if (neph[isys] > MAXEPH)
                    {
                        printf(
                                "***ERROR(read_rnxnav):exceed the maxium ephemeris number!\n");
                        exit(1);
                    }
                    memcpy(&ephg[neph[isys] - 1], &ephr, sizeof(ephr));
                }
                break;
            default:
                fgets(line, LEN_STRING, fp);
                break;
        }
    }
    fclose(fp);
}

void read_rnxobs(FILE *fp, int jd0, double sod0, int nprn0,
                 char (*cprn0)[LEN_PRN], int *nfreq,
                 char freq[MAXSYS][MAXFREQ][LEN_FREQ], RNXHEAD *HD, RNXOBS *OB)
{
    char line[1024] = {'\0'};
    char varword[LEN_STRING], code[LEN_OBSTYPE] = {'\0'};
    int nprn = 0, iflag = 0, i, j, k, isat, isys;
    int iy, im, id, ih, imi, lastpos, nline, ntemp, count, iobs;
    double obs[MAXOBSTYP];
    double sec, ds;
    char cprn[MAXSAT][LEN_PRN];
    for (isat = 0; isat < MAXSAT; isat++)
        strncpy(cprn[isat], "\0", LEN_PRN);
    int lfind = false;
    for (i = 0; i < CKF.nprn; i++)
    {
        for (j = 0; j < 2 * MAXFREQ; j++)
            OB->obs[i][j] = 0.0;
    }
    while (!lfind)
    {
        lastpos = ftell(fp);
        fgets(line, LEN_STRING, fp);
        if (feof(fp))
        {
            return;
        }
        if (len_trim(line) == 0)
            continue;
        if (HD->ver < 3.0)
        {
            substringEx(varword, line, 29, 3);
            nprn = atoi(varword);
            substringEx(varword, line, 26, 3);
            iflag = atoi(varword);
        }
        else
        {
            substringEx(varword, line, 32, 3);
            nprn = atoi(varword);
            substringEx(varword, line, 29, 3);
            iflag = atoi(varword);
        }
        if (nprn > MAXSAT)
        {
            printf("satellite number > MAXSAT!\n");
            exit(1);
        }
        if (iflag > 1)
        {
            if (HD->ver < 3.0)
            {
                for (i = 0; i < nprn; i++)
                {
                    fgets(line, LEN_STRING, fp);
                    if (strstr(line, "ANTENNA: DELTA H/E/N") != NULL)
                        sscanf(line, "%lf%lf%lf", &HD->h, &HD->e, &HD->n);
                }
                continue;
            }
        }
        // initialize obs
        for (i = 0; i < MAXSAT; i++)
            for (j = 0; j < 2 * MAXFREQ; j++)
                OB->obs[i][j] = 0.0;
        if (HD->ver >= 3.0)
            sscanf(line, "%*s%d%d%d%d%d%lf", &iy, &im, &id, &ih, &imi, &sec);
        else
            sscanf(line, "%d%d%d%d%d%lf", &iy, &im, &id, &ih, &imi, &sec);
        yr2year(&iy);
        OB->jd = modified_julday(iy, im, id);
        OB->tsec = nint(ih * 3600.0 + imi * 60.0 + sec);
        if (HD->ver < 3.0)
        {
            if (nprn % 12 != 0)
                nline = (int)(nprn / 12.0) + 1;
            else
                nline = nprn / 12;
            if (nprn == 0)
                nline = 1;
        }

        if (jd0 != 0)
        {
            ds = timdif(jd0, sod0, OB->jd, OB->tsec);
            // the requested time  is behind the current time in the obsfile
            if (ds < -MAXWND)
            {
                fseek(fp, lastpos - ftell(fp), SEEK_CUR);
                /*printf(
                 "###Warning(read_rnx_obs): lost time %8.2lf!\n",
                 sod0, OB->tsec);*/
                return;
            }
            else if (ds > MAXWND)
            {
                if (HD->ver < 3.0)
                {
                    for (i = 0; i < nline - 1; i++)
                        fgets(line, LEN_STRING, fp);
                    if (HD->nobstype[0] % 5 != 0)
                        ntemp = nprn * ((int)(HD->nobstype[0] / 5) + 1);
                    else
                        ntemp = nprn * (HD->nobstype[0] / 5);

                    for (i = 0; i < ntemp; i++)
                        fgets(line, LEN_STRING, fp);
                }
                else
                {
                    for (i = 0; i < nprn; i++)
                        fgets(line, LEN_STRING, fp);
                }
            }
            else
            {
                lfind = true;
                CKF.mjd = OB->jd;
                CKF.sod = ih * 3600.0 + imi * 60.0 + sec;
            }
        }
        else
        {
            break;
        }
    }

    if (HD->t0[0] == 0)
    {
        HD->t0[0] = iy;
        HD->t0[1] = im;
        HD->t0[2] = id;
        HD->t0[3] = ih;
        HD->t0[4] = imi;
        HD->t0[5] = (int)sec;
    }
    HD->t1[0] = iy;
    HD->t1[1] = im;
    HD->t1[2] = id;
    HD->t1[3] = ih;
    HD->t1[4] = imi;
    HD->t1[5] = (int)sec;

    if (HD->ver < 3.0)
    {
        fseek(fp, lastpos - ftell(fp), SEEK_CUR);
        for (i = 0; i < nline; i++)
        {
            fgets(line, LEN_STRING, fp);
            // get the correct cprn list
            for (j = 0; j < (MIN((nprn - 12 * i), 12)); j++)
            {
                substringEx(varword, line, 32 + j * 3, 3);
                if (varword[0] == ' ')
                    varword[0] = 'G';
                if (varword[1] == ' ')
                    varword[1] = '0';
                strcpy(cprn[j + i * 12], varword);
            }
        }
        if (HD->nobstype[0] % 5 != 0)
            nline = (int)(HD->nobstype[0] / 5.0) + 1;
        else
            nline = HD->nobstype[0] / 5;
        // for some observations maybe empty,so decided to split string
        for (i = 0; i < nprn; i++)
        {
            memset(obs, 0, sizeof(double) * MAXOBSTYP);
            for (j = 0; j < nline; j++)
            {
                memset(line, 0, sizeof(char) * LEN_STRING);
                fgets(line, LEN_STRING, fp);
                fillobs(line, MIN((HD->nobstype[0] - j * 5), 5), 16, HD->ver);
                for (k = 0; k < (MIN((HD->nobstype[0] - j * 5), 5)); k++)
                {
                    //	memset(varword,0,sizeof(char)*LEN_STRING);
                    substringEx(varword, line, 16 * k, 14);
                    obs[k + j * 5] = atof(varword);
                }
            }
            isat = -1;
            if (nprn0 > 0)
                isat = pointer_string(nprn0, LEN_PRN, (char *)cprn0, cprn[i]);
            else
                isat = i;
            if (isat != -1)
            {
                isys = index_string(SYS, cprn[i][0]);
                for (j = 0; j < nfreq[isys]; j++)
                {
                    OB->obs[isat][MAXFREQ + j] = 0.0;
                    OB->obs[isat][j] = 0.0;
                    code[0] = 'P';
                    code[1] = freq[isys][j][1];
                    code[2] = '\0';
                    k = pointer_string(HD->nobstype[isys], LEN_OBSTYPE,
                                       (char *)HD->obstype[isys], code);
                    if (k != -1)
                    {
                        OB->obs[isat][MAXFREQ + j] = obs[k];
                        OB->fob[isat][MAXFREQ + j][0] = 'C';
                        OB->fob[isat][MAXFREQ + j][1] = freq[isys][j][1];
                        OB->fob[isat][MAXFREQ + j][2] = 'P';
                    }
                    if (OB->obs[isat][MAXFREQ + j] == 0)
                    {
                        code[0] = 'C';
                        code[1] = freq[isys][j][1];
                        code[2] = '\0';
                        k = pointer_string(HD->nobstype[isys], LEN_OBSTYPE,
                                           (char *)HD->obstype[isys], code);
                        if (k != -1)
                        {
                            OB->obs[isat][MAXFREQ + j] = obs[k];
                            OB->fob[isat][MAXFREQ + j][0] = 'C';
                            OB->fob[isat][MAXFREQ + j][1] = freq[isys][j][1];
                            OB->fob[isat][MAXFREQ + j][2] = 'C';
                        }
                    }
                    code[0] = 'L';
                    code[1] = freq[isys][j][1];
                    code[2] = '\0';
                    k = pointer_string(HD->nobstype[isys], LEN_OBSTYPE,
                                       (char *)HD->obstype[isys], code);
                    if (k != -1)
                    {
                        OB->obs[isat][j] = obs[k];
                        OB->fob[isat][j][0] = 'L';
                        OB->fob[isat][j][1] = freq[isys][j][1];
                        OB->fob[isat][j][2] = 'P';
                    }
                }
                count = 0;

                if (strstr(CKF.cobs, "SF") == NULL)
                {
                    for (k = 0; k < 2 * MAXFREQ; k++)
                        if (OB->obs[isat][k] != 0)
                            count++;
                    // remember to initialize obs at the end of the circulation
                    ds = ABS(OB->obs[isat][MAXFREQ] - OB->obs[isat][MAXFREQ + 1]);
                    // please be careful the therhold 50.d0, it is effeced by ionosphere
                    if (count < 4 || ds > 50)
                    {
                        for (k = 0; k < 2 * MAXFREQ; k++)
                        {
                            OB->obs[isat][k] = 0.0;
                            memset(OB->fob[isat][k], 0, sizeof(char) * LEN_OBSTYPE);
                        }
                    }
                }
                else
                {
                    if (OB->obs[isat][0] == 0 || OB->obs[isat][MAXFREQ] == 0)
                    {
                        OB->obs[isat][0] = 0.0;
                        memset(OB->fob[isat][0], 0, sizeof(char) * LEN_OBSTYPE);

                        OB->obs[isat][MAXFREQ] = 0.0;
                        memset(OB->fob[isat][MAXFREQ], 0, sizeof(char) * LEN_OBSTYPE);
                    }
                }
            }
        }
    }

    if (HD->ver >= 3.0)
    {
        isat = -1;
        for (i = 0; i < nprn; i++)
        {
            memset(line, 0, sizeof(char) * LEN_STRING);
            fgets(line, LEN_STRING, fp);
            isys = index_string(SYS, line[0]);
            if (-1 == isys)
                continue;
            if (nprn0 > 0)
            {
                if (line[1] == ' ')
                    line[1] = '0';
                // clearstring(varword);
                memset(varword, 0, sizeof(varword));
                substringEx(varword, line, 0, 3);
                isat = pointer_string(nprn0, LEN_PRN, (char *)cprn0, varword);
                // if isat equals -1 then the observations is totally disordered
                // the entire program can't be runned corrected
                if (isat == -1)
                    continue;
                strcpy(cprn[isat], varword);
            }
            else
            {
                // the final is the index (len-1)
                isat += 1;
                memset(varword, 0, sizeof(varword));
                substringEx(varword, line, 0, 3);
                strcpy(cprn[isat], varword);
            }
            substringEx(varword, line, 3, strlen(line) - 3);
            fillobs(varword, HD->nobstype[isys], 16, HD->ver);
            strcpy(line + 3, varword);
            memset(obs, 0, sizeof(double) * MAXOBSTYP);
            for (j = 0; j < HD->nobstype[isys]; j++)
            {
                substringEx(varword, line, 3 + j * 16, 14);
                obs[j] = atof(varword);
            }
            if (isat != -1)
            {
                for (j = 0; j < nfreq[isys]; j++)
                {
                    OB->obs[isat][MAXFREQ + j] = 0.0;
                    OB->obs[isat][j] = 0.0;
                    // remember to initialize the HD->usetype
                    if (!OB->lstored[isys][j])
                    {
                        // because of the \0
                        for (iobs = 0; iobs < strlen(OBSTYPE); iobs++)
                        {
                            code[0] = 'C';
                            code[1] = freq[isys][j][1];
                            code[2] = OBSTYPE[iobs];
                            code[3] = '\0';
                            k = pointer_string(HD->nobstype[isys], LEN_OBSTYPE,
                                               (char *)HD->obstype[isys], code);
                            if (k != -1)
                            {
                                if (fabs(obs[k]) > 1.0)
                                {
                                    OB->obs[isat][MAXFREQ + j] = obs[k];
                                    strcpy(OB->fob[isat][MAXFREQ + j],
                                           code);
                                }
                            }
                            code[0] = 'L';
                            code[1] = freq[isys][j][1];
                            code[2] = OBSTYPE[iobs];
                            code[3] = '\0';
                            k = pointer_string(HD->nobstype[isys], LEN_OBSTYPE,
                                               (char *)HD->obstype[isys], code);
                            if (k != -1)
                            {
                                if (fabs(obs[k]) > 1.0)
                                {
                                    OB->obs[isat][j] = obs[k];
                                    strcpy(OB->fob[isat][j], code);
                                }
                            }
                            if (OB->obs[isat][j] != 0 && OB->obs[isat][MAXFREQ + j] != 0)
                                break;
                        }
                        if (iobs < strlen(OBSTYPE))
                        {
                            HD->usetype[isys][j] = OBSTYPE[iobs];
                            OB->lstored[isys][j] = true;
                        }
                    }
                    else
                    {

                        code[0] = 'C';
                        code[1] = freq[isys][j][1];
                        code[2] = HD->usetype[isys][j];
                        code[3] = '\0';
                        k = pointer_string(HD->nobstype[isys], LEN_OBSTYPE,
                                           (char *)HD->obstype[isys], code);
                        if (k != -1)
                        {
                            if (fabs(obs[k]) > 1.0)
                            {
                                OB->obs[isat][MAXFREQ + j] = obs[k];
                                strcpy(OB->fob[isat][MAXFREQ + j], code);
                            }
                        }
                        code[0] = 'L';
                        code[1] = freq[isys][j][1];
                        code[2] = HD->usetype[isys][j];
                        code[3] = '\0';
                        k = pointer_string(HD->nobstype[isys], LEN_OBSTYPE,
                                           (char *)HD->obstype[isys], code);
                        if (k != -1)
                        {
                            if (fabs(obs[k]) > 1.0)
                            {
                                OB->obs[isat][j] = obs[k];
                                strcpy(OB->fob[isat][j], code);
                            }
                        }
                    }
                }
                count = 0;
                if (strstr(CKF.cobs, "SF") == NULL)
                {
                    for (k = 0; k < 2 * MAXFREQ; k++)
                        if (OB->obs[isat][k] != 0)
                            count++;
                    // remember to initialize obs at the end of the circulation
                    ds = ABS(OB->obs[isat][MAXFREQ] - OB->obs[isat][MAXFREQ + 1]);
                    // please be careful the therhold 50.d0, it is effeced by ionosphere
                    if (count < 4 || ds > 50)
                    {
                        for (k = 0; k < 2 * MAXFREQ; k++)
                        {
                            OB->obs[isat][k] = 0.0;
                            memset(OB->fob[isat][k], 0, sizeof(char) * LEN_OBSTYPE);
                        }
                    }
                }
                else
                {
                    if (OB->obs[isat][0] == 0 || OB->obs[isat][MAXFREQ + 0] == 0)
                    {
                        OB->obs[isat][0] = 0.0;
                        memset(OB->fob[isat][0], 0, sizeof(char) * LEN_OBSTYPE);

                        OB->obs[isat][MAXFREQ] = 0.0;
                        memset(OB->fob[isat][MAXFREQ], 0, sizeof(char) * LEN_OBSTYPE);
                    }
                }
            }
        }
    }
    if (nprn0 == 0)
    {
        OB->nprn = nprn;
        for (i = 0; i < nprn; i++)
            strcpy(OB->cprn[i], cprn[i]);
    }
    else
    {
        OB->nprn = nprn0;
        for (i = 0; i < nprn0; i++)
            strcpy(OB->cprn[i], cprn0[i]);
    }
}

void read_rnxobs_head(RNXHEAD *HD, FILE *fp)
{
    int i, j, id, lastpos;
    char line[LEN_STRING] = {'\0'};
    char keyword[LEN_STRING];
    char varword[LEN_STRING];
    HD->nsys = 0;
    memset(HD->nobstype, 0, sizeof(HD->nobstype));
    memset(HD->obstype, 0, sizeof(HD->obstype));

    while (true)
    {
        lastpos = ftell(fp);
        fgets(line, LEN_STRING, fp);
        substringEx(keyword, line, 60, strlen(line) - 60);
        if (strstr(keyword, "END OF HEADER") != NULL) //strstr()未找到子字符串返回NULL
            return;
        if ((strstr(keyword, "END OF HEADER") != NULL && HD->ver >= 2.0) || (len_trim(keyword) == 0 && HD->ver < 2.0))
            return;

        if (strstr(keyword, "RINEX VERSION") != NULL)
        {
            //版本号
            substringEx(varword, line, 0, 9);
            HD->ver = atof(varword);

            //导航电文类型
            substringEx(varword, line, 40, 1);
            HD->sys[0] = varword[0];
            HD->sys[1] = '\0';

            if (HD->sys[0] == ' ')
                HD->sys[0] = 'G'; //缺省为GPS
            if (HD->ver < 1.0 || HD->ver > 4.0)
                printf("read rinex version error!\n");
        }
        else if (strstr(keyword, "COMMENT") != NULL)
            continue;
        else if (strstr(keyword, "MARKER NAME") != NULL)
        {
            substringEx(varword, line, 0, 4);
            strncpy(HD->mark, varword, 4);
        }
        else if (strstr(keyword, "REC #") != NULL)
        {

            substringEx(varword, line, 0, 20);
            strncpy(HD->recnum, varword, 20);
            substringEx(varword, line, 20, 20);
            strncpy(HD->rectype, varword, 20);
        }
        else if (strstr(keyword, "ANT #") != NULL)
        {
            substringEx(varword, line, 0, 20);
            strncpy(HD->antnum, varword, 20);
            substringEx(varword, line, 20, 20);
            strncpy(HD->anttype, varword, 20);
        }
        else if (strstr(keyword, "APPROX POSITION") != NULL)
        {
            sscanf(line, "%lf%lf%lf", &HD->x, &HD->y, &HD->z);
        }
        else if (strstr(keyword, "ANTENNA: DELTA") != NULL)
        {
            sscanf(line, "%lf%lf%lf", &HD->h, &HD->e, &HD->n);
        }
        else if (strstr(keyword, "WAVELENGTH FACT") != NULL)
        {
        }
        else if (strstr(keyword, "SYS / # / OBS TYPES") != NULL)
        {
            id = index_string(SYS, line[0]);
            if (-1 == id)
                continue;
            substringEx(varword, line, 3, 3);
            HD->nobstype[id] = atoi(varword);
            int nline = 0;
            if (HD->nobstype[id] % 13 != 0)
                nline = (int)(HD->nobstype[id] / 13) + 1;
            else
                nline = (int)(HD->nobstype[id] / 13);
            //maybe can use fseek
            fseek(fp, lastpos - ftell(fp), SEEK_CUR);
            for (j = 0; j < nline; j++)
            {
                memset(line, 0, sizeof(line));
                fgets(line, LEN_STRING, fp);

                for (i = 0; i < (MIN((HD->nobstype[id] - 13 * j), 13)); i++)
                {
                    substringEx(varword, line, 7 + i * 4, 3);
                    strcpy(HD->obstype[id][i + 13 * j], varword);
                }
            }
            for (i = HD->nobstype[id]; i < MAXOBSTYP; i++)
                HD->obstype[id][i][0] = '\0';
        }
        else if (strstr(keyword, "TYPES OF OBSERV") != NULL)
        {

            substringEx(varword, line, 0, 6);
            HD->nobstype[0] = atoi(varword);

            if (HD->nobstype[0] == 0)
                continue;
            int nline = 0;
            if (HD->nobstype[0] % 9 != 0)
                nline = (int)(HD->nobstype[0] / 9) + 1;
            else
                nline = (int)(HD->nobstype[0] / 9);

            fseek(fp, lastpos - ftell(fp), SEEK_CUR);
            for (j = 0; j < nline; j++)
            {
                memset(line, 0, sizeof(line));
                fgets(line, LEN_STRING, fp);
                for (i = 0; i < (MIN((HD->nobstype[0] - 9 * j), 9)); i++)
                {
                    substringEx(varword, line, 6 + 6 * i, 6);
                    strcpy(HD->obstype[0][i + 9 * j], left_justify_string(varword));
                }
            }
            for (i = HD->nobstype[0]; i < MAXOBSTYP; i++)
                HD->obstype[0][i][0] = '\0';

            for (i = 1; i < MAXSYS; i++)
            {
                HD->nobstype[i] = HD->nobstype[0];
                for (j = 0; j < HD->nobstype[0]; j++)
                    strcpy(HD->obstype[i][j], HD->obstype[0][j]);
            }

            //for whu's tracking stations,the C1 is C2,C2 is C7
            id = index_string(SYS, 'C');
            if (-1 != id)
            {
                for (i = 0; i < HD->nobstype[id]; i++)
                {
                    if ('2' == HD->obstype[id][i][1])
                        HD->obstype[id][i][1] = '7';
                    if ('1' == HD->obstype[id][i][1])
                        HD->obstype[id][i][1] = '2';
                }
            }
        }
        else if (strstr(keyword, "INTERVAL") != NULL)
        {
            substringEx(varword, line, 0, 16);
            HD->intv = atof(varword);
        }
        else if (strstr(keyword, "TIME OF FIRST OBS") != NULL)
        {
            double sec;
            char tmetag[256] = {0};
            sscanf(line, "%d%d%d%d%d%lf%s", HD->t0, HD->t0 + 1, HD->t0 + 2, HD->t0 + 3,
                   HD->t0 + 4, &sec, tmetag);
            HD->t0[5] = (int)sec;
        }
        else if (strstr(keyword, "TIME OF LAST OBS") != NULL)
        {
            for (i = 0; i < 5; i++)
            {
                substringEx(varword, line, 0 + i * 6, 6);
                HD->t1[i] = atoi(varword);
            }
            substringEx(varword, line, 30, 13);
            HD->t1[5] = (int)atof(varword);
        }
        memset(line, 0, sizeof(line));
    }
    if (HD->ver >= 3.0)
    {
        for (i = 0; i < MAXSYS; i++)
            if (HD->nobstype[i] != 0)
                HD->nsys++;
    }
}