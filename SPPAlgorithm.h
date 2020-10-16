//
// Created by eureka on 2020/10/16.
//

#ifndef SPP_SPPALGORITHM_H
#define SPP_SPPALGORITHM_H
#include "CMat.h"
#include "SatPos.h"


int Adjust(double *B, double *P, double *l, int countSat, double *x, double tol);

void SPP(Coord *sitPos, Coord *satPos, double *obs, double *corr, int countSat, double *elev);

#endif //SPP_SPPALGORITHM_H
