/*
 * global.h
 *
 *  Created on: 2020年9月20日
 *      Author: xps
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_
#include "rnx.h"
#ifdef GLOBAL_VAR
#define GLOBAL_EXTERN
#else
#define GLOBAL_EXTERN extern
#endif

GLOBAL_EXTERN char SYS[16];
GLOBAL_EXTERN char OBSTYPE[17];
GLOBAL_EXTERN CKDCFG CKF; //TODO what?

GLOBAL_EXTERN int neph[MAXSYS]; //导航星历总数
GLOBAL_EXTERN GLONASS_BRDEPH ephgls[MAXEPH];
GLOBAL_EXTERN GPS_BRDEPH ephgps[MAXSYS][MAXEPH];

static  double CLIGHT = 299792458.0;
static  double GM = 3.986005e14;
static  double OMEGA_E = 7.2921151467e-5; //地球自转角速度

#endif /* GLOBAL_H_ */
