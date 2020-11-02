/*
 * all kinds of corrections for observation */
#ifndef SPP_CORR_H
#define SPP_CORR_H

#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rnx.h"
#include "global.h"
#include <math.h>
#include "SatPos.h"


double SatPos_CorrEarthRot(Coord *satPos, Coord *groundRef);

double SatPos_CorrRelative(Coord *satPos, Coord *satPos_v);

extern void ecef2pos(Coord *r, double *pos);

extern void rot_enu2xyz(double lat, double lon, double *rotmat);

void matmul(const char *tr, int m, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta,
            double *c, int ldc);

int azel_demo();

double SatClk(GPS_BRDEPH *gps_brd, double t0);

double TGD(char cprn[4]);

double ionmodel(gtime_t t, const double *ion, const double *pos,
                const double *azel);

double tropmodel(gtime_t time, const double *pos, const double *azel,
                 double humi);
#endif //SPP_CORR_H
