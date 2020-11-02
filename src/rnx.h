/*
 * rnx.h
 * read renix file
 */

#ifndef RNX_H_
#define RNX_H_

#include <stdio.h>

#define MIN(a, b)   ((a>=b) ? b : a)
#define MAX(a, b)   ((a>=b) ? a : b)
#define ABS(a)     ((a>=0) ? a : -(a))
#define SIGN(a, b)  ((b>=0) ? a:-(a))
#define nint(a)    ((int)(a+SIGN(1,a)*0.5))

#define true        1
#define false        0
#define MAXSIT        150
#define MAXSYS        7
#define LEN_PRN        4
#define MAXFREQ        3
#define LEN_STRING    1024
#define LEN_OBSTYPE 4
#define MAXSAT        150
#define MAXOBSTYP    50
#define LEN_ANTENNA        21
#define MAXWND 0.01
#define LEN_FREQ 4
#define MAXEPH     (MAXSAT * 24)
#define LEN_EPHION        4

//////////////////////////
typedef struct {
    double ver;
    char sys[4];
    char mark[4];
    char rectype[LEN_ANTENNA];
    char anttype[LEN_ANTENNA];
    char recnum[LEN_ANTENNA];
    char antnum[LEN_ANTENNA];
    double x, y, z, h, e, n;
    int fact1, fact2, nobstype[MAXSYS];
    char obstype[MAXSYS][MAXOBSTYP][LEN_OBSTYPE];
    double intv;
    int nprn;
    char cprn[MAXSAT][LEN_PRN];
    int t0[6], t1[6];
    int nsys;
    char tsys[8];
    char usetype[MAXSYS][MAXFREQ];
} RNXHEAD;

typedef struct {
    int nprn;
    char cprn[MAXSAT][LEN_PRN];
    char obstyp[MAXOBSTYP][LEN_OBSTYPE];
    int jd;
    double tsec;
    double obs[MAXSAT][2 * MAXFREQ];
    char fob[MAXSAT][2 * MAXFREQ][LEN_OBSTYPE];
    int lstored[MAXSYS][MAXFREQ]; //eturn bool
} RNXOBS;

typedef struct {
    char cobs[20];
    // the GNSS system L1 OR L2 nfreq is the count
    //int nsys;
    //char system[MAXSYS];
    int nfreq[MAXSYS];
    char freq[MAXSYS][MAXFREQ][LEN_FREQ];

    // time of epoch and interval
    int mjd;
    double sod;

    // stations and satellites
    int nprn;
    char cprn[MAXSAT][LEN_PRN]; //the cprn in session.obj
} CKDCFG;

// GPS Ephemeris
typedef struct {
    char cprn[4];
    int mjd;
    double sod;
    // a[0]: SV clock offset
    // a[1]: SV clock drift
    // a[2]: SV clock drift rate
    double a0, a1, a2;
    // aode: age of ephemeris upload
    // crs, crc: Ortital radius correction
    // dn: Mean motion difference
    // m0: Mean anomaly at reference epoch
    // e: Eccentricity
    // cuc, cus: Latitude argument correction
    // roota: Square root of semi-major axis
    double aode, crs, dn;
    double m0, cuc, e;
    double cus, roota;
    // toe, week: Ephemerides reference epoch in seconds with the week
    // cis, cic: Inclination correction
    // omega0: Longtitude of ascending node at the begining of the week
    // i0: Inclination at reference epoch
    // omega: Argument of perigee
    // omegadot: Rate of node's right ascension
    double toe, cic, Omega0;
    double cis, i0, crc;
    double omega, Omega_dot;
    // idot: Rate of inclination angle
    // sesvd0:
    // resvd1:
    // accu: SV accuracy
    // hlth: SV health
    // tgd: Time group delay
    // aodc: Age of clock parameter upload
    double i_dot, resvd0, week, resvd1;
    double accu, hlth, tgd, tgd1, aodc;
    // added for b2b
    int signal_idx;
    unsigned int iodc;
    double tgd_BDS[13], isc_BDS[13];
    double delta_A, A_DOT, delta_n_dot;
} GPS_BRDEPH;

// GLONASS Ephemeris
typedef struct {
    char cprn[4];
    int mjd;
    double sod;
    // tau: SV clock bias
    // gama: SV relative frequency bias
    // tk: message frame time (tk+nd*86400) in seconds of the UTC week
    // pos: coordinate at ephemerides reference epoch in PZ-90
    // vel: velocity at ephemerides reference epoch in PZ-90
    // acc: acceleration at ephemerides reference epoch in PZ-90
    // health: SV health
    // frenum: frequency number
    // age: age of operation information
    double tau;
    double gamma;
    double tk;
    double pos[3];
    double vel[3];
    double acc[3];
    double health;
    double frenum;
    double age;
    double aode;
} GLONASS_BRDEPH;

// GNSS NAVIGATION MESSAGE FILE     HEADER SECTION
typedef struct {
    double ver;
    // ionospheric correction parameters
    char ionc[MAXSYS][2][LEN_EPHION];
    double ion[MAXSYS][2][4];
    // corrections to transform the system to UTC or other time systems

    char timc[MAXSYS][2][LEN_EPHION];
    double tim[MAXSYS][2][4];
    // number of leap second since 6-Jan-980
    int leap;
} BRDHEAD;

void initialize();

extern void read_rnxobs(FILE *fp, int jd0, double sod0, int nprn0,
                        char (*cprn0)[LEN_PRN], int *nfreq,
                        char freq[MAXSYS][MAXFREQ][LEN_FREQ], RNXHEAD *HD, RNXOBS *OB);

void read_rnxnav(char csys, const char *flnbrd, double mjd0, double mjd1, BRDHEAD *hd,
                 int *neph, GPS_BRDEPH ephm[MAXSYS][MAXEPH], GLONASS_BRDEPH ephg[MAXEPH]);

extern void read_rnxobs_head(RNXHEAD *HD, FILE *fp);

extern char *trim(char *pStr);

extern int len_trim(char *pStr);

extern int pointer_string(int row, int col, char *string_array, char *string);

extern int index_string(char *src, char key);

extern char *substringEx(char *dest, char *src, int start, int length);

extern void fillobs(char *line, int nobs, int itemlen, double ver);

extern void filleph(char *line, double ver);

extern char *left_justify_string(char *string);

void mjd2wksow(int mjd, double sod, int *week, double *sow);

int modified_julday(int iyear, int imonth, int iday);

void timinc(int jd, double sec, double delt, int *jd1, double *sec1);

double timdif(int jd2, double sod2, int jd1, double sod1);

void yeardoy2monthday(int iyear, int idoy, int *imonth, int *iday);

void wksow2mjd(int week, double sow, int *mjd, double *sod);

void mjd2date(int jd, double sod, int *iyear, int *imonth, int *iday, int *ih, int *imin, double *sec);

void mjd2doy(int jd, int *iyear, int *idoy);

void yr2year(int *yr);

void Swap(double *a, double *b);

void SortEph();

#endif /* RNX_H_ */
