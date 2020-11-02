//
// Created by eureka on 2020/11/1.
//
#include "corr.h"

double SatPos_CorrEarthRot(Coord *satPos, Coord *groundRef) {
    double delta_rho = -
                               OMEGA_E / CLIGHT *
                       (satPos->Y * (groundRef->X - satPos->X) - satPos->X * (groundRef->Y - satPos->Y));
    return delta_rho;
}

double SatPos_CorrRelative(Coord *satPos, Coord *satPos_v) {
    double delta_rho = 2.0 / CLIGHT * (satPos->X * satPos_v->X +
                                       satPos->Y * satPos_v->Y +
                                       satPos->Z * satPos_v->Z);
    return delta_rho;
}

double SatClk(GPS_BRDEPH *gps_brd, double t0) {
    double dt = t0 - gps_brd->toe;
    return (gps_brd->a0 + gps_brd->a1 * dt + gps_brd->a2 * dt * dt) * CLIGHT;
}

double TGD(char cprn[4]){
    FILE *fp = fopen("../data/tgd_sav", "r");
    char prn[4];
    char buf[12];
    double tgd=0;
    while (fgets(buf, 12, fp) != NULL)
    {
        if (buf[1]==cprn[1] && buf[2]==cprn[2]){
            sscanf(buf,"%s %lf\n", &prn, &tgd);
            break;
        }
    }
    fclose(fp);
    return tgd;

}



/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/

/*euclid norm*/
double norm(const double *a,int n)
{
    double res=0;
    for (int i = 0; i < n; i++) {
        res+=a[i]*a[i];
    }
    return sqrt(res);
}


double ionmodel(gtime_t t, const double *ion, const double *pos, const double *azel)
{
    const double ion_default[]={ /* 2004/1/1 */
            0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
            0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
    };
    double tt,f,psi,phi,lam,amp,per,x;
    int week;

    if (pos[2]<-1E3||azel[1]<=0) return 0.0;
    if (norm(ion,8)<=0.0) ion=ion_default;

    /* earth centered angle (semi-circle) */
    psi=0.0137/(azel[1]/PI+0.11)-0.022;

    /* subionospheric latitude/longitude (semi-circle) */
    phi=pos[0]/PI+psi*cos(azel[0]);
    if      (phi> 0.416) phi= 0.416;
    else if (phi<-0.416) phi=-0.416;
    lam=pos[1]/PI+psi*sin(azel[0])/cos(phi*PI);

    /* geomagnetic latitude (semi-circle) */
    phi+=0.064*cos((lam-1.617)*PI);

    /* local time (s) */
    double a = time2gpst(t, &week);
    tt=43200.0*lam + time2gpst(t, &week);
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

    /* slant factor */
    f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);

    /* ionospheric delay */
    amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
    per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
    amp=amp<    0.0?    0.0:amp;
    per=per<72000.0?72000.0:per;
    x=2.0*PI*(tt-50400.0)/per;

    return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
}


/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by standard atmosphere and saastamoinen model
* args   : gtime_t time     I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double humi      I   relative humidity
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
double tropmodel(gtime_t time, const double *pos, const double *azel,
                 double humi)
{
    const double temp0=15.0; /* temparature at sea level */
    double hgt,pres,temp,e,z,trph,trpw;

    if (pos[2]<-100.0||1E4<pos[2]||azel[1]<=0) return 0.0;

    /* standard atmosphere */
    hgt=pos[2]<0.0?0.0:pos[2];

    pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
    temp=temp0-6.5E-3*hgt+273.16;
    e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));

    /* saastamoninen model */
    z=PI/2.0-azel[1];
    trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3)/cos(z);
    trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
    return trph+trpw;
}




