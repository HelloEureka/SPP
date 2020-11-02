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

void elevazel(Coord *sat, Coord *rec, double *azel);

double Dist(Coord *p1, Coord *p2);

