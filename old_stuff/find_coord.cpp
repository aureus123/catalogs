
/*
 * FIND_COORD - Busca por coordenadas en catálogos antiguos
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "read_old.h"
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
    printf("FIND_COORD - Find near stars around given coordinate in %s.\n",
        TRANSFORM ? "J2000.0" : "B1875.0 FK4");
    printf("Made in 2024 by Daniel Severin.\n");

    if (argc < 7) {
        printf("\nUsage: find_coord RAh RAm 100*RAs -Decld Declm 10*Decls\n");
        printf("\nExample: find_coord 12 35 4789 24 15 83\n");
        printf("  find by coordinates 12h 35m 47s89 -24° 15' 8''3 (%d)\n",
            TRANSFORM ? 2000 : 1875);
        return -1;
    }

    /* producimos la coordenada */
    int RAh = atoi(argv[1]);
    if (RAh < 0 || RAh > 23) {
        printf("RAh must be between 0 and 23\n");
        return -1;
    }
    int RAm = atoi(argv[2]);
    if (RAm < 0 || RAm > 59) {
        printf("RAm must be between 0 and 59\n");
        return -1;
    }
    int RAs = atoi(argv[3]);
    if (RAs < 0 || RAs > 5999) {
        printf("RAs must be between 0 and 5999\n");
        return -1;
    }
    int Decld = atoi(argv[4]);
    if (Decld < 0 || Decld > 89) {
        printf("Decld must be between 0 and 89\n");
        return -1;
    }
    int Declm = atoi(argv[5]);
    if (Declm < 0 || Declm > 59) {
        printf("Declm must be between 0 and 59\n");
        return -1;
    }
    int Decls = atoi(argv[6]);
    if (Decls < 0 || Decls > 599) {
        printf("Decls must be between 0 and 599\n");
        return -1;
    }

    if (TRANSFORM) {
        /* transformamos de J2000 y a B1875 */
        double RA = (double) RAh + ((double) RAm)/60.0 + (((double) RAs)/100.0)/3600.0;
        double Decl = (double) Decld + ((double) Declm)/60.0 + (((double) Decls)/10.0)/3600.0;
        double pmRA = 0.0;
        double pmDecl = 0.0;
        RA *= 15.0; Decl = -Decl;
        wcsconp(WCS_J2000, WCS_B1950, 0.0, 1875.0, 2000.001278, 1875.0, &RA, &Decl, &pmRA, &pmDecl);
        RA /= 15.0; Decl = -Decl;
        double dRAh = floor(RA);
        RAh = (int) dRAh;
        double dRAm = floor((RA-dRAh) * 60.0);
        RAm = (int) dRAm;
        double dRAs = floor((((RA-dRAh) * 60.0) - dRAm) * 6000.0);
        RAs = (int) dRAs;
        double dDecld = floor(Decl);
        Decld = (int) dDecld;
        double dDeclm = floor((Decl-dDecld) * 60.0);
        Declm = (int) dDeclm;
        double dDecls = floor((((Decl-dDecld) * 60.0) - dDeclm) * 600.0);
        Decls = (int) dDecls;
    }

    /* leemos catalogo CD */
    readDM("cat/cd.txt");
    struct DMstar_struct *CDstar = getDMStruct();

    /* leemos catalogo GC */
	readGC(false, RAh, RAm, RAs, Decld, Declm, Decls);
	struct GCstar_struct *GCstar = getGCStruct();
    int cdIndex = GCstar[0].cdIndex;
    printf("*) Nearest CD %d°%d (mag=%.1f) separated in %.1f arcsec.\n",
        CDstar[cdIndex].declRef,
        CDstar[cdIndex].numRef,
        CDstar[cdIndex].vmag,
        GCstar[0].dist);
	writeRegisterGC(0);
    return 0;
}
