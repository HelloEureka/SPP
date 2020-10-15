#include "SatPos.h"

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n) {
    double c = 0.0;

    while (--n >= 0)
        c += a[n] * b[n];

    return c;
}

/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void ecef2pos(Coord *r, double *pos) {
    double e2 = FE_WGS84 * (2.0 - FE_WGS84);
    double r2 = r->X * r->X + r->Y * r->Y ; //+ r->Z * r->Z;
    double z, zk, v = RE_WGS84, sinp;
    for (z = r->Z, zk = 0.0; fabs(z - zk) >= 1E-4;) {
        zk = z;
        sinp = z / sqrt(r2 + z * z);
        v = RE_WGS84 / sqrt(1.0 - e2 * sinp * sinp);
        z = r->Z + v * e2 * sinp;
    }
    pos[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (r->Z > 0.0 ? PI / 2.0 : -PI / 2.0);
    pos[1] = r2 > 1E-12 ? atan2(r->Y, r->X) : 0.0;
    pos[2] = sqrt(r2 + z * z) - v;
}

extern void rot_enu2xyz(double lat, double lon, double *rotmat) {
    double coslat, sinlat, coslon, sinlon;
    coslat = cos(lat - PI / 2);
    sinlat = sin(lat - PI / 2);
    coslon = cos(-PI / 2 - lon);
    sinlon = sin(-PI / 2 - lon);

    rotmat[0] = coslon;
    rotmat[1] = -sinlon;
    rotmat[2] = 0;

    rotmat[3] = sinlon * coslat;
    rotmat[4] = coslon * coslat;
    rotmat[5] = -sinlat;

    rotmat[6] = sinlon * sinlat;
    rotmat[7] = coslon * sinlat;
    rotmat[8] = coslat;
}

void matmul(const char *tr, int m, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta,
            double *c, int ldc) {
    double d;
    int i, j, x;
    int opr = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            d = 0;
            switch (opr) {
                case 1:
                    for (x = 0; x < k; x++)
                        d += a[lda * x + i] * b[ldb * j + x];
                    break;
                case 2:
                    for (x = 0; x < k; x++)
                        d += a[lda * x + i] * b[ldb * x + j];
                    break;
                case 3:
                    for (x = 0; x < k; x++)
                        d += a[lda * i + x] * b[ldb * j + x];
                    break;
                case 4:
                    for (x = 0; x < k; x++)
                        d += a[lda * i + x] * b[ldb * x + j];
                    break;
            }
            c[ldc * j + i] = alpha * d + beta * c[ldc * j + i];
        }
    }
}

void elevazel(Coord *sat, Coord *rec, double *azel) {
    double r[3], rot_l2f[9], dump[3], pos[3];
    ecef2pos(rec, pos);
    rot_enu2xyz(pos[0], pos[1], rot_l2f);
    r[0] = sat->X - rec->X;
    r[1] = sat->Y - rec->Y;
    r[2] = sat->Z - rec->Z;
    matmul("NN", 1, 3, 3, 1.0, r, 1, rot_l2f, 3, 0.0, dump, 1);
    azel[0] = atan2(dump[0], dump[1]);
    azel[1] = atan(dump[2] / sqrt(dump[0] * dump[0] + dump[1] * dump[1]));
}

int azel_demo() {
    double azel[2];
    Coord rr= {-2341333.112, -3539049.520, 4745791.248};
    Coord rs = {-3425190.629, -17360562.979, 20466189.117};
    elevazel(&rs, &rr, azel);

    printf("azim:%9.3lf elev %9.3lf\n", azel[0] * RAD2DEG, azel[1] * RAD2DEG);
}