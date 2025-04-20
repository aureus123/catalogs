/*
 * TRANSFORM - Transforma coordenadas de J2000 a Bxxxx (FK4)
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "misc.h"

/* Para uso de la libreria WCS: */
#define WCS_J2000 1 /* J2000(FK5) right ascension and declination */
#define WCS_B1950 2 /* B1950(FK4) right ascension and declination */
extern "C" void wcsconp(int sys1, int sys2, double eq1, double eq2, double ep1, double ep2,
             double *dtheta, double *dphi, double *ptheta, double *pphi);

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("TRANSFORM - Transform coordinates from J2000 to Besselian year.\n");
    printf("Made in 2025 by Daniel Severin.\n");

    if (argc < 8) {
        printf("\nUsage: transform xxxx RAh RAm 100*RAs +Decld Declm 10*Decls\n");
        printf("\nExample: transform 1875 12 35 4789 +24 15 83\n");
        printf("  finds the B1875 coordinates of 12h 35m 47s89 -24° 15' 8''3\n");
        return -1;
    }

    /* producimos la coordenada */
    int targetYear = atoi(argv[1]);
    if (targetYear < 1700 || targetYear > 2100) {
        printf("targetYear must be between 1700 and 2100\n");
        return -1;
    }
    int RAh = atoi(argv[2]);
    if (RAh < 0 || RAh > 23) {
        printf("RAh must be between 0 and 23\n");
        return -1;
    }
    int RAm = atoi(argv[3]);
    if (RAm < 0 || RAm > 59) {
        printf("RAm must be between 0 and 59\n");
        return -1;
    }
    int RAs = atoi(argv[4]);
    if (RAs < 0 || RAs > 5999) {
        printf("RAs must be between 0 and 5999\n");
        return -1;
    }
    bool negativeDecl = argv[5][0] == '-' ? true : false;
    int Decld = atoi(&argv[5][1]);
    if (Decld < 0 || Decld > 89) {
        printf("Decld must be between 0 and 89\n");
        return -1;
    }
    int Declm = atoi(argv[6]);
    if (Declm < 0 || Declm > 59) {
        printf("Declm must be between 0 and 59\n");
        return -1;
    }
    int Decls = atoi(argv[7]);
    if (Decls < 0 || Decls > 599) {
        printf("Decls must be between 0 and 599\n");
        return -1;
    }
    
    /* transformamos de J2000 a Bxxxx */
    double RA = (double) RAh + ((double) RAm)/60.0 + (((double) RAs)/100.0)/3600.0;
    double Decl = (double) Decld + ((double) Declm)/60.0 + (((double) Decls)/10.0)/3600.0;
    double pmRA = 0.0;
    double pmDecl = 0.0;
    RA *= 15.0; /* <-- RA en grados */
    if (negativeDecl) Decl = -Decl;
    wcsconp(WCS_J2000, WCS_B1950, 0.0, (float) targetYear, 2000.001278, (float) targetYear, &RA, &Decl, &pmRA, &pmDecl);
    RA /= 15.0;  /* RA <-- RA en horas */
    double dRAh = floor(RA);
    RAh = (int) dRAh;
    double dRAm = floor((RA-dRAh) * 60.0);
    RAm = (int) dRAm;
    double dRAs = floor((((RA-dRAh) * 60.0) - dRAm) * 6000.0);
    RAs = (int) dRAs;
    if (Decl < 0) {
        Decl = -Decl;
        negativeDecl = true;
    } else {
        negativeDecl = false;
    }
    double dDecld = floor(Decl);
    Decld = (int) dDecld;
    double dDeclm = floor((Decl-dDecld) * 60.0);
    Declm = (int) dDeclm;
    double dDecls = floor((((Decl-dDecld) * 60.0) - dDeclm) * 600.0);
    Decls = (int) dDecls;

    printf("Coordinates in B%d: %02dh %02dm %02ds%02d %c%02d° %02d' %02d''%01d\n",
           targetYear, RAh, RAm, RAs / 100, RAs % 100, negativeDecl ? '-' : '+', Decld, Declm, Decls / 10, Decls % 10);

    return 0;
}
