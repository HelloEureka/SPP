/*
 * Satellite position and velocity
 * in ECEF
 * */

#include "SatPos.h"

/*
 * Calculate E by iterating with Newton's Method
 * */
double Calculate_E(double E_initial, double M, double e, double tol, int iter) {
    double E = E_initial;
    int i;
    double err;
    for (i = 0; i < iter; i++) {
        err = (E - e * sin(E) - M) / (1 - e * cos(E));
        if (fabs(err) > tol)
            E -= err;
        else
            break;
    }
    return E;
}


/*卫星轨道参数改正及计算*/
Coord SatPos_Cal(GPS_BRDEPH *gps_brd, double t) {
    //1. 计算卫星运动的平均角速度n
    double n = sqrt(GM) / pow(gps_brd->roota, 3) + gps_brd->dn;

    //2. 计算观测瞬间卫星平近点角M
    double dt = t - gps_brd->toe;
    double M = gps_brd->m0 + n * dt;

    //3. 计算偏近点角E
    double E = Calculate_E(M, M, gps_brd->e, 1e-20, 10); //初始值为M

    //4. 计算真近点角f
    double f = atan2((sqrt(1 - gps_brd->e * gps_brd->e) * sin(E)), (cos(E) - gps_brd->e));

    //5. 升交角距u_prime
    double u_prime = gps_brd->omega + f;

    //6. 计算摄动改正项
    double cos_2u_prime = cos(2 * u_prime);
    double sin_2u_prime = sin(2 * u_prime);
    double delta_u = gps_brd->cuc * cos_2u_prime + gps_brd->cus * sin_2u_prime;
    double delta_r = gps_brd->crc * cos_2u_prime + gps_brd->crs * sin_2u_prime;
    double delta_i = gps_brd->cic * cos_2u_prime + gps_brd->cis * sin_2u_prime;

    //7. 对u_prime, r_prime, i_0进行摄动改正
    double u = u_prime + delta_u;
    double r = gps_brd->roota * gps_brd->roota * (1 - gps_brd->e * cos(E)) + delta_r;
    double i = gps_brd->i0 + delta_i + gps_brd->i_dot * dt;

    //8. 计算观测瞬间升交点经度L
    double L = gps_brd->Omega0 + gps_brd->Omega_dot * dt - OMEGA_E * t;

    /* 卫星轨道面坐标计算 */
    double x = r * cos(u);
    double y = r * sin(u);

    /* 轨道面坐标 -> WGS84 */
    Coord satPos = {0, 0, 0};
    satPos.X = x * cos(L) - y * sin(L) * cos(i);
    satPos.Y = x * sin(L) + y * cos(L) * cos(i);
    satPos.Z = y * sin(i);

    return satPos;
}

double Dist(Coord *p1, Coord *p2) {
    double dX = p1->X - p2->X;
    double dY = p1->Y - p2->Y;
    double dZ = p1->Z - p2->Z;
    double dist = sqrt(dX * dX + dY * dY + dZ * dZ);
    return dist;
}


/*
 * t: time on ground
 * considering of signal spreading delay(about 0.07s)
 * */
Coord SatPos_CorrSpread(GPS_BRDEPH *gps_brd, Coord *groundRef, double t) {
    double t_spread, t_spread_prime, s2g_dist;
    t_spread_prime = 0;
    Coord satPos = SatPos_Cal(gps_brd, t - 0.07); //signal spreads for about 0.07s
    for (int i = 0; i < 5; i++) {
        s2g_dist = Dist(&satPos, groundRef);
        t_spread = s2g_dist / CLIGHT;
        if (fabs(t_spread_prime - t_spread) < 1e-10) {
            break;
        }
        t_spread_prime = t_spread;
        satPos = SatPos_Cal(gps_brd, t - t_spread);
    }
    return satPos;
}


Coord SatVel(GPS_BRDEPH *gps_brd, double t) {
    double dt = 1e-8;
    Coord posLast = SatPos_Cal(gps_brd, t - dt);
    Coord posNext = SatPos_Cal(gps_brd, t + dt);
    Coord vel = {posNext.X - posLast.X, posNext.Y - posLast.Y, posNext.Z - posLast.Z};
    vel.X /= (2 * dt);
    vel.Y /= (2 * dt);
    vel.Z /= (2 * dt);
    return vel;
}


