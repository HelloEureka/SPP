#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rnx.h"
#include "global.h"
#include <math.h>

typedef struct {
    double X;
    double Y;
    double Z;
} Coord;

Coord SatPos_Cal(GPS_BRDEPH *gps_brd, double t);

double Calculate_E(double E_initial, double M, double e, double tol, int iter);

Coord SatPos_CorrSpread(GPS_BRDEPH *gps_brd, Coord *groundRef, double t);

Coord SatVel(GPS_BRDEPH *gps_brd,  double t);

double SatPos_CorrEarthRot(Coord *satPos, Coord *groundRef);

double SatPos_CorrRelative(Coord *satPos, Coord *satPos_v);

double Dist(Coord *p1, Coord *p2);





extern double dot(const double *a, const double *b, int n);

extern void ecef2pos(Coord *r, double *pos);

extern void rot_enu2xyz(double lat, double lon, double *rotmat);

void matmul(const char *tr, int m, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta,
            double *c, int ldc);

void elevazel(Coord *sat, Coord *rec, double *azel);

int azel_demo();

double SatClk(GPS_BRDEPH *gps_brd, double t0);

double TGD(char cprn[4]);

double ionmodel(gtime_t t, const double *ion, const double *pos,
                const double *azel);

double tropmodel(gtime_t time, const double *pos, const double *azel,
                 double humi);