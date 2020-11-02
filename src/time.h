/* time.h
 * different from built-in <time.h>
 * Created by eureka on 2020/11/1.
 * */
#ifndef SPP_TIME_H
#define SPP_TIME_H

#include <time.h>

typedef struct {        /* time struct */
    time_t time;        /* time (s) expressed by standard time_t */
    double sec;         /* fraction of second under 1 s */
} gtime_t;

static const double gpst0[]={1980,1,6,0,0,0}; //GPS
static const double gst0[]={1999,8,22,0,0,0};   //GALILEO
static const double bdt0[]={2006,1,1,0,0,0};    //BEIDOU


gtime_t epoch2time(const double *ep);
gtime_t gpst2time(int week, double sec);
double time2gpst(gtime_t t, int *week);







#endif //SPP_TIME_H
