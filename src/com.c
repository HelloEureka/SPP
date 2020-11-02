/*
 * com.c
 * common functions for getting data from renix
 */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rnx.h"
#include "global.h"

void Swap(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}


char *trim(char *pStr) {
    int len = strlen(pStr);
    char *pIndex = pStr + len - 1;
    while (*pIndex == '\n' || *pIndex == '\r' || *pIndex == '\0' || *pIndex == ' ') {
        *pIndex = '\0';
        pIndex--;
    }
    return pStr;
}

//计算有效字符数, 去除空格换行回车
int len_trim(char *pStr) {
    int length = strlen(pStr);
    int count = length;
    int i;
    for (i = length - 1; i >= 0; i--) {
        if (pStr[i] == '\0' || pStr[i] == '\n' || pStr[i] == '\r' || pStr[i] == ' ')
            //				空格				换行				回车			??? 和空格一样
            count--;
        else
            break;
    }
    return count;
}


int pointer_string(int row, int col, char *string_array, char *string) {
    int i;
    char *pStr = (char *) string_array;
    trim(string);

    for (i = 0; i < row; i++) {
        if (strcmp(pStr + i * col, string) == 0)
            break;
    }
    if (i == row)
        i = -1;
    return i;
}

//from left to right
int index_string(char *src, char key) {
    int len = strlen(src);
    int i;
    for (i = 0; i < len; i++) {
        if (src[i] == key)
            break;
    }
    if (i == len)
        return -1;
    else
        return i;
}

// start: started with index of zero
char *substringEx(char *dest, char *src, int start, int length) {
    int i, j = 0;
    int len = strlen(src);
    if (start < 0 || start >= len || start + length > len) {
        dest[0] = '\0';
        return NULL;
    }

    for (i = start; i < start + length; i++) {
        dest[j] = src[i];
        j++;
    }

    dest[j] = '\0';
    return dest;
}

void filleph(char *line, double ver) {
    int i, len, cpre = 4, ncount;
    char tmp[128];
    char xline[1024];
    if (ver < 3.0)
        cpre = 3;
    len = len_trim(line);
    for (i = len; i < cpre + 19 * 4; i++) {
        line[i] = ' ';
    }
    line[cpre + 19 * 4] = '\0';
    for (i = 0; i < 4; i++) {
        ncount = cpre + 19 * (i + 1);
        if (strlen(line) < ncount) {
            line[cpre + 19 * i + 1] = '0';
        } else {
            substringEx(tmp, line + cpre + 19 * i, 0, 19);
            if (len_trim(tmp) == 0) {
                line[cpre + 19 * i + 1] = '0';
            }
        }
    }
}

void fillobs(char *line, int nobs, int itemlen, double ver) {
    int OFFSET = 0, len, i;
    char tmp[256];
    if (ver > 3.0)
        OFFSET = 3;
    len = len_trim(line);
    for (i = len; i < itemlen * nobs + OFFSET; i++)
        line[i] = ' ';
    line[itemlen * nobs + OFFSET] = '\0';
    for (i = 0; i < nobs; i++) {
        memset(tmp, 0, sizeof(char) * 256);
        substringEx(tmp, line + OFFSET + i * itemlen, 0, itemlen - 2); // last is signal strength
        if (len_trim(tmp) == 0) {
            line[itemlen * i + OFFSET] = '0';
        }
    }
}

void mjd2doy(int jd, int *iyear, int *idoy) {
    *iyear = (jd + 678940) / 365;
    *idoy = jd - modified_julday(*iyear, 1, 1);
    while (*idoy < 0) {
        (*iyear)--;
        *idoy = jd - modified_julday(*iyear, 1, 1) + 1;
    }
}

char *left_justify_string(char *string) {
    int p = 0;
    while (*(string + p) == ' ' || *(string + p) == '\n' || *(string + p) == '\r')
        p++;
    return string + p;
}

void yr2year(int *yr) {
    if (*yr > 1900)
        return;
    if (*yr <= 20)
        *yr += 2000;
    else
        *yr += 1900;
}

void mjd2date(int jd, double sod, int *iyear, int *imonth, int *iday, int *ih,
              int *imin, double *sec) {
    int doy = 0;
    mjd2doy(jd, iyear, &doy);
    yeardoy2monthday(*iyear, doy, imonth, iday);

    *ih = (int) sod / 3600.0;
    *imin = (int) ((sod - (*ih) * 3600.0) / 60.0);
    *sec = sod - (*ih) * 3600.0 - (*imin) * 60.0;
}

void wksow2mjd(int week, double sow, int *mjd, double *sod) {
    if (mjd != NULL)
        *mjd = (int) (sow / 86400.0) + week * 7 + 44244;
    if (sod != NULL)
        *sod = fmod(sow, 86400.0);
}

void yeardoy2monthday(int iyear, int idoy, int *imonth, int *iday) {
    int days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int id, i;
    if ((iyear % 4 == 0 && iyear % 100 != 0) || iyear % 400 == 0)
        days_in_month[1] = 29;
    id = idoy;
    for (i = 0; i < 12; i++) {
        id = id - days_in_month[i];
        if (id > 0)
            continue;
        *iday = id + days_in_month[i];
        *imonth = i + 1;
        break;
    }
}

double timdif(int jd2, double sod2, int jd1, double sod1) {
    return 86400.0 * (jd2 - jd1) + sod2 - sod1;
}

int modified_julday(int iyear, int imonth, int iday) {
    int iyr, result;
    int doy_of_month[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304,
                            334};
    if (iyear < 0 || imonth < 0 || iday < 0 || imonth > 12 || iday > 366 || (imonth != 0 && iday > 31)) {
        printf("***ERROR(modified_julday): incorrect arguments!%d %d %d\n",
               iyear, imonth, iday);
        exit(1);
    }
    iyr = iyear;
    if (imonth <= 2)
        iyr -= 1;
    result = 365 * iyear - 678941 + iyr / 4 - iyr / 100 + iyr / 400 + iday;

    if (imonth != 0)
        result = result + doy_of_month[imonth - 1];
    return result;
}

void mjd2wksow(int mjd, double sod, int *week, double *sow) {
    *week = (int) ((mjd + sod / 86400.0 - 44244.0) / 7.0);
    *sow = (mjd - 44244.0 - *week * 7) * 86400.0 + sod;
}


/*
 * Sort ephemeris by time
 */
void SortEph() {
    for (int i = 0; i < MAXEPH - 1; i++) {
        for (int j = 0; j < MAXEPH - 1 - i; ++j) {
            if ((ephgps[0][j].cprn[0] != 0)
                && (ephgps[0][j].cprn[1] == ephgps[0][j + 1].cprn[1])
                && (ephgps[0][j].cprn[2] == ephgps[0][j + 1].cprn[2])
                && (ephgps[0][j].toe > ephgps[0][j + 1].toe)
                && (ephgps[0][j].mjd == ephgps[0][j + 1].mjd)) {
                Swap(&(ephgps[0][j].toe), &(ephgps[0][j + 1].toe));
                Swap(&(ephgps[0][j].sod), &(ephgps[0][j + 1].sod));
                Swap(&(ephgps[0][j].a0), &(ephgps[0][j + 1].a0));
                Swap(&(ephgps[0][j].a1), &(ephgps[0][j + 1].a1));
                Swap(&(ephgps[0][j].a2), &(ephgps[0][j + 1].a2));
                Swap(&(ephgps[0][j].aode), &(ephgps[0][j + 1].aode));
                Swap(&(ephgps[0][j].crs), &(ephgps[0][j + 1].crs));
                Swap(&(ephgps[0][j].dn), &(ephgps[0][j + 1].dn));
                Swap(&(ephgps[0][j].m0), &(ephgps[0][j + 1].m0));
                Swap(&(ephgps[0][j].e), &(ephgps[0][j + 1].e));
                Swap(&(ephgps[0][j].cus), &(ephgps[0][j + 1].cus));
                Swap(&(ephgps[0][j].roota), &(ephgps[0][j + 1].roota));
                Swap(&(ephgps[0][j].cic), &(ephgps[0][j + 1].cic));
                Swap(&(ephgps[0][j].Omega0), &(ephgps[0][j + 1].Omega0));
                Swap(&(ephgps[0][j].cic), &(ephgps[0][j + 1].cic));
                Swap(&(ephgps[0][j].i0), &(ephgps[0][j + 1].i0));
                Swap(&(ephgps[0][j].cic), &(ephgps[0][j + 1].cic));
                Swap(&(ephgps[0][j].omega), &(ephgps[0][j + 1].omega));
                Swap(&(ephgps[0][j].Omega_dot), &(ephgps[0][j + 1].Omega_dot));
                Swap(&(ephgps[0][j].i_dot), &(ephgps[0][j + 1].i_dot));
                Swap(&(ephgps[0][j].resvd0), &(ephgps[0][j + 1].resvd0));
                Swap(&(ephgps[0][j].week), &(ephgps[0][j + 1].week));
                Swap(&(ephgps[0][j].resvd1), &(ephgps[0][j + 1].resvd1));
                Swap(&(ephgps[0][j].accu), &(ephgps[0][j + 1].accu));
                Swap(&(ephgps[0][j].hlth), &(ephgps[0][j + 1].hlth));
                Swap(&(ephgps[0][j].tgd), &(ephgps[0][j + 1].tgd));
                Swap(&(ephgps[0][j].aodc), &(ephgps[0][j + 1].aodc));
                Swap(&(ephgps[0][j].delta_A), &(ephgps[0][j + 1].delta_A));
                Swap(&(ephgps[0][j].A_DOT), &(ephgps[0][j + 1].A_DOT));
                Swap(&(ephgps[0][j].delta_n_dot), &(ephgps[0][j + 1].delta_n_dot));
            }
        }
    }
}