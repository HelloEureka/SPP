#include "SatPos.h"

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
    // printf("------------------calculate E, Newton's Method--------------------\n");
    // printf("iter = %d \n", i);
    // printf("计算偏近点角 E = %f\n", E);
    // if (i == iter)
    //     printf("    no convergence guarantee\n");
    // else
    //     printf("    convergence guarantee\n");
    // printf("-------------------------------------------------------------------\n\n");
    return E;
}

Coord SatPos_Cal(GPS_BRDEPH *gps_brd, double t) {
    /*卫星轨道参数改正及计算*/

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
    // double cL = cos(L);
    // double sL = sin(L);
    // double ci = cos(i);
    // double si = sin(i);
    // satPos.X = x * cL - y * sL * ci;
    // satPos.Y = x * sL + y * cL * ci;
    // satPos.Z = y * si;

    satPos.X = x * cos(L) - y * sin(L) * cos(i);
    satPos.Y = x * sin(L) + y * cos(L) * cos(i);
    satPos.Z = y * sin(i);

    //中间过程
    // printf("计算卫星运动的平均角速度 n = %.15f \n", n);
    // printf("计算观测瞬间卫星平近点角 M = %.15f\n", M);
    // printf("E = %.15f\n", E);
    // printf("计算真近点角 f = %f\n", f);
    // printf("升交角距 u_prime = %f\n", u_prime);
    // printf("\n--------------计算摄动改正项--------------------\n", M);
    // printf("delta_u = %f\n", delta_u);
    // printf("delta_r = %f\n", delta_r);
    // printf("delta_i = %f\n", delta_i);
    // printf("------------------------------------------------\n", M);
    // printf("--------------对u_prime, r_prime, i_0进行摄动改正---------------\n", M);
    // printf("u = %f\n", u);
    // printf("r = %f\n", r);
    // printf("i = %f\n", i);
    // printf("计算观测瞬间升交点经度 L = %f\n", L);
    // printf("\n----------------卫星轨道面坐标------------------\n");
    // printf("x = %f\n", x);
    // printf("y = %f\n", y);

    // printf("\n------------------------WGS84坐标-------------------------\n");

    return satPos;
}

double Dist(Coord *p1, Coord *p2) {
    double dX = p1->X - p2->X;
    double dY = p1->Y - p2->Y;
    double dZ = p1->Z - p2->Z;
    double dist = sqrt(dX * dX + dY * dY + dZ * dZ);
    return dist;
}

Coord SatPos_CorrSpread(GPS_BRDEPH *gps_brd, Coord *groundRef, double t) {
    double t_spread, t_spread_prime, s2g_dist;
    t_spread_prime = 0;
    Coord satPos = SatPos_Cal(gps_brd, t - 0.07); //signal spreads for about 0.07s
    for (int i = 0; i < 10000; i++) {
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

double SatPos_CorrEarthRot(Coord *satPos, Coord *groundRef) {
    double delta_rho = -
            OMEGA_E / CLIGHT * (satPos->Y * (groundRef->X - satPos->X) - satPos->X * (groundRef->Y - satPos->Y));
    return delta_rho;
    //以下直接改卫星位置
//    double s2g_dist = Dist(&satPos, groundRef);
    //double tau = s2g_dist / clight;
    //double alpha = omega_e * tau;
    //double sin_alpha = sin(alpha);
    //double cos_alpha = cos(alpha);
    // double X = satPos->X;
    // satPos->X = X * cos_alpha + satPos->Y * sin_alpha;
    // satPos->Y = -sin_alpha * X + cos_alpha * satPos->Y;
}

double SatPos_CorrRelative(Coord *satPos, Coord *satPos_v) {
    double delta_rho = 2.0 / CLIGHT * (satPos->X * satPos_v->X +
                                        satPos->Y * satPos_v->Y +
                                        satPos->Z * satPos_v->Z);
    return delta_rho;
}


Coord SatVel(GPS_BRDEPH *gps_brd, Coord *groundRef, double t) {
    double dt = 1e-8;
    Coord posLast = SatPos_CorrSpread(gps_brd, groundRef, t - dt);
    Coord posNext = SatPos_CorrSpread(gps_brd, groundRef, t + dt);
    Coord vel = {posNext.X - posLast.X, posNext.Y - posLast.Y, posNext.Z - posLast.Z};
    vel.X /= (2 * dt);
    vel.Y /= (2 * dt);
    vel.Z /= (2 * dt);
    return vel;
}