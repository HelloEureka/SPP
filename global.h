/*
 * global.h
 *
 *  Created on: 2020年9月20日
 *      Author: xps
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_
#include "rnx.h"
#include <time.h>
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

typedef struct {        /* time struct */
    time_t time;        /* time (s) expressed by standard time_t */
    double sec;         /* fraction of second under 1 s */
} gtime_t;


#define PI 3.1415926535897932384626433832795028841971693993
#define DEG2RAD  (PI/180.0)
#define RAD2DEG  (180.0/PI)
#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */


static  double CLIGHT = 299792458.0;
static  double GM = 3.986005e14;
static  double OMEGA_E = 7.2921151467e-5; //地球自转角速度
static const double gpst0[]={1980,1,6,0,0,0}; //GPS
static const double gst0[]={1999,8,22,0,0,0};   //GALILEO
static const double bdt0[]={2006,1,1,0,0,0};    //BEIDOU


gtime_t epoch2time(const double *ep);
#endif /* GLOBAL_H_ */
