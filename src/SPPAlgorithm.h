//
// Created by eureka on 2020/10/16.
//

#ifndef SPP_SPPALGORITHM_H
#define SPP_SPPALGORITHM_H

#include "CMat.h"
#include "corr.h"

void GetData(RNXOBS *OB, int *obsCountAva, Coord *satPos, Coord *groundRef, double t, double *obs_c1w, double *corr,
             double *elev);

int Adjust(double *B, double *P, double *l, int countSat, double *x, double tol, double *sigma);

void SPP(Coord *sitPos, Coord *satPos, double *obs, double *corr, int countSat, double *elev, double *sigma);

double Sigma(double *B, double *P, double *l, int countSat, int satAvalid, double *x);

#endif //SPP_SPPALGORITHM_H
