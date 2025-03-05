
/*
 * GEN_TYCHO2 - Genera una identificación cruzada entre Tycho-2 y catálogos antiguos
 * El flag isCD() determina si es para el hemiferio norte o sur (en J2000)
 * Para el hemisferio norte en J2000, primero se utiliza PPM y luego BD (pero solo,
 * hasta la declinación +19, esto quiere decir que de la declinación +20 en adelante
 * solo se reportan las estrellas de BD cruzadas con PPM).
 * Para el hemisferio sur en J2000, primero se utiliza PPM y luego CD (obsérvese que
 * no están BD/SD por lo que sólo se mencionan estos catálogos a través de la
 * identificación cruzada con PPM). Por otra parte, se utilizan las identificaciones
 * del Catálogo General Argentino (GC), Oeltzen-Argelander (OA) y las de Gilliss (G).
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "read_ppm.h"
#include "trig.h"
#include "misc.h"

/* Para uso de la libreria WCS: */
#define WCS_J2000 1 /* J2000(FK5) right ascension and declination */
#define WCS_B1950 2 /* B1950(FK4) right ascension and declination */
extern "C" void wcsconp(int sys1, int sys2, double eq1, double eq2, double ep1, double ep2,
             double *dtheta, double *dphi, double *ptheta, double *pphi);

#define STRING_SIZE 14

static char dmStringDM[MAXDMSTAR][STRING_SIZE];

/*
 * readCrossFile - lee archivos de identificaciones cruzadas en formato CSV
 */
void readCrossFile(const char *ppm_file, struct PPMstar_struct *PPMstar, int PPMstars, const char *cd_file, struct DMstar_struct *DMstar, int DMstars) {
    char buffer[1024];
    char targetRef[STRING_SIZE];
    int ppmRef, cdDeclRef, cdNumRef;
    float minDistance;

    /* read cross file between PPM and target */
    printf("Reading cross file %s... ", ppm_file);
    FILE *stream = fopen(ppm_file, "rt");
    if (stream == NULL) {
        snprintf(buffer, 1024, "Cannot read %s", ppm_file);
        perror(buffer);
        exit(1);
    }
    bool first_line = true;
    while (fgets(buffer, 1023, stream) != NULL) {
        if (first_line) {
            // omit first line (header)
            first_line = false;
            continue;
        }
        sscanf(buffer, "%13[^,],PPM %d,%f\n", targetRef, &ppmRef, &minDistance);
        // printf("Cross: %s, PPM %d, dist = %.1f arcsec.\n", targetRef, ppmRef, minDistance);
        if (minDistance < __FLT_EPSILON__ || minDistance > 15.0) {
            // omit identifications with zero distance (bug) or too far away
            continue;
        }
        for (int i = 0; i < PPMstars; i++) {
            if (PPMstar[i].ppmRef == ppmRef) {
                strncpy(PPMstar[i].dmString, targetRef, STRING_SIZE);
                break;
            }
        }
    }
    fclose(stream);
    printf("done!\n");

    /* read cross file between CD and target (optative) */
    if (DMstars == 0) return;
    printf("Reading cross file %s... ", cd_file);
    stream = fopen(cd_file, "rt");
    if (stream == NULL) {
        snprintf(buffer, 1024, "Cannot read %s", cd_file);
        perror(buffer);
        exit(1);
    }
    first_line = true;
    while (fgets(buffer, 1023, stream) != NULL) {
        if (first_line) {
            // omit first line (header)
            first_line = false;
            continue;
        }
        sscanf(buffer, "%13[^,],CD %d°%d,%f\n", targetRef, &cdDeclRef, &cdNumRef, &minDistance);
        // printf("Cross: %s, CD %d°%d, dist = %.1f arcsec.\n", targetRef, cdDeclRef, cdNumRef, minDistance);
        if (minDistance < __FLT_EPSILON__ || minDistance > 30.0) {
            // omit identifications with zero distance (bug) or too far away
            continue;
        }
        for (int i = 0; i < DMstars; i++) {
            if (DMstar[i].declRef == cdDeclRef && DMstar[i].numRef == cdNumRef) {
                strncpy(dmStringDM[i], targetRef, STRING_SIZE);
                break;
            }
        }
    }
    fclose(stream);
    printf("done!\n");
}

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv) {
    FILE *stream;
    char buffer[1024];
    char cell[256];

    /* leemos catalogo BD/CD */
    readDM(isCD() ? "cat/cd_curated.txt" : "cat/bd_curated.txt");
    struct DMstar_struct *DMstar = getDMStruct();
    int DMstars = getDMStars();
    bool is_north = !isCD();

    /* generamos identificadores de DM */
    for (int i = 0; i < DMstars; i++) {
        snprintf(dmStringDM[i], STRING_SIZE, "%cD %c%d°%d",
                isCD() ? 'C' : 'B',
                DMstar[i].signRef ? '-' : '+',
                abs(DMstar[i].declRef),
                DMstar[i].numRef);
    }

    /* leemos catalogo PPM (2000 en B1950) */
    readPPM(false, true, !is_north, is_north, 2000.0);
    struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    if (isCD()) {
        /* leemos identificaciones cruzadas y renombramos designaciones */
        readCrossFile("results/cross_gilliss_ppm.csv", PPMstar, PPMstars, "results/cross_gilliss_cd.csv", DMstar, DMstars);
        readCrossFile("results/cross_usno_ppm.csv", PPMstar, PPMstars, "results/cross_usno_cd.csv", DMstar, DMstars);
        readCrossFile("results/cross_oa_ppm.csv", PPMstar, PPMstars, "results/cross_oa_cd.csv", DMstar, DMstars);
        readCrossFile("results/cross_gc_ppm.csv", PPMstar, PPMstars, "results/cross_gc_cd.csv", DMstar, DMstars);
    } else {
        /* leemos identificaciones cruzadas y renombramos designaciones, solo para PPM-USNO */
        readCrossFile("results/cross_usno_ppm.csv", PPMstar, PPMstars, "", nullptr, 0);
    }

    stream = fopen("cat/tyc2.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read tyc2.txt");
        exit(1);
    }

    snprintf(buffer, 1024, "tycho2/cross_tyc2_%s.csv", is_north ? "north" : "south");
    FILE *crossStream = openCrossFile(buffer);

    int TYCstarsPPM = 0;
    int TYCstarsDM = 0;
    int TYCunidentified = 0;
    int entry = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
        // if (TYCstarsPPM > 100) break;
		entry++;
		if (entry % 10000 == 0) {
			printf("Progress (%.2f%%): PPM = %d, DM = %d, unidentified = %d\n",
                (100.0 * (float) entry) / 2539913.0,
                TYCstarsPPM,
                TYCstarsDM,
                TYCunidentified);
		}

        /* lee declinación y descarta tempranamente (con un márgen de 2 grados) */
        readField(buffer, cell, 166, 12);
        double Decl = atof(cell);
        if (is_north) {
            if (Decl < -2.0) continue;
        } else {
            if (Decl > +2.0) continue;
        }

        /* lee numeración */
        readField(buffer, cell, 1, 4);
        int tyc1Ref = atoi(cell);
        readField(buffer, cell, 6, 5);
        int tyc2Ref = atoi(cell);
        readField(buffer, cell, 12, 1);
        int tyc3Ref = atoi(cell);

        char tycString[STRING_SIZE];
        snprintf(tycString, STRING_SIZE, "TYC %d-%d-%d", tyc1Ref, tyc2Ref, tyc3Ref);

        /* lee RA y Decl (epoch) */
        readField(buffer, cell, 153, 12);
        double RA = atof(cell);
        readField(buffer, cell, 179, 4);
        double epRA = atof(cell);
        readField(buffer, cell, 184, 4);
        double epDecl = atof(cell);
        double epoch = 1990.0 + (epRA + epDecl) / 2.0;

        /* lee mov. propios, si los tiene */
        double pmRA = 0.0;
        double pmDecl = 0.0;
        readField(buffer, cell, 14, 1);
        if (cell[0] != 'X') {
            readField(buffer, cell, 42, 7);
            pmRA = atof(cell);
            readField(buffer, cell, 50, 7);
            pmDecl = atof(cell);
            
            pmRA /= 1000 * 3600 * dcos(epDecl); /* conversion de mas/yr a grados/yr (juliano) */
            pmDecl /= 1000 * 3600; /* conversion de mas/yr a grados/yr (juliano) */
        }

        // printf("TYC %d-%d-%d: RA = %f, Decl = %f, pmRA = %f, pmDecl = %f, epoch = %f\n", tyc1Ref, tyc2Ref, tyc3Ref, RA, Decl, pmRA, pmDecl, epoch);

        /* convertir a J2000 y verificar si se descarta por hemisferio
           TODO: convertir epoch a de Julian a Besselian, y usar 2000.001278 */
        double RAtarget = RA;
        double Decltarget = Decl;
        double pmRAtarget = pmRA;
        double pmDecltarget = pmDecl;
        wcsconp(WCS_J2000, WCS_J2000, 0.0, 0.0, epoch, 2000.0, &RAtarget, &Decltarget, &pmRAtarget, &pmDecltarget);
        if (is_north) {
            if (Decltarget < 0.0) continue;
        } else {
            if (Decltarget > 0.0) continue;
        }

        /* convertir a 2000 (B1950) */
        RAtarget = RA;
        Decltarget = Decl;
        pmRAtarget = pmRA;
        pmDecltarget = pmDecl;
        wcsconp(WCS_J2000, WCS_B1950, 0.0, 2000.0, epoch, 2000.0, &RAtarget, &Decltarget, &pmRAtarget, &pmDecltarget);
        // printf("    RA = %f, Decl = %f, pmRA = %f, pmDecl = %f, epoch = %f\n", RAtarget, Decltarget, pmRAtarget, pmDecltarget, 2000.0);

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RAtarget, Decltarget, &x, &y, &z);

        /* halla PPM más cercana, dentro de 15 arcsec */
        int ppmIndex = -1;
        double minDistance = HUGE_NUMBER;
        findPPMByCoordinates(x, y, z, &ppmIndex, &minDistance);
        if (ppmIndex != -1 && minDistance < 15.0) {
            if (PPMstar[ppmIndex].dmString[0] != 0) {
                /* se almacena la identificación cruzada con la DM dada por PPM */
                writeCrossEntry(crossStream, tycString, PPMstar[ppmIndex].dmString, minDistance);
                TYCstarsPPM++;
                continue;
            }
        }

        /* convertir a 1875 o 1855 */
        RAtarget = RA;
        Decltarget = Decl;
        pmRAtarget = pmRA;
        pmDecltarget = pmDecl;
        double targetYear = isCD() ? 1875.0 : 1855.0;
        wcsconp(WCS_J2000, WCS_B1950, 0.0, targetYear, epoch, targetYear, &RAtarget, &Decltarget, &pmRAtarget, &pmDecltarget);
        
        /* calcula coordenadas rectangulares */
        sph2rec(RAtarget, Decltarget, &x, &y, &z);

        /* halla la DM más cercana, dentro de 60 arcsec */
        int dmIndex = -1;
        minDistance = HUGE_NUMBER;
        findDMByCoordinates(x, y, z, &dmIndex, &minDistance);
        if (dmIndex != -1 && minDistance < 60.0) {
            /* se almacena la identificación cruzada con la DM */
            writeCrossEntry(crossStream, tycString, dmStringDM[dmIndex], minDistance);
            TYCstarsDM++;
            continue;
        }
        TYCunidentified++;
    }
    printf("Stars read and identified of Tycho-2 from PPM: %d\n", TYCstarsPPM);
    printf("Stars read and identified of Tycho-2 from DM: %d\n", TYCstarsDM);
    printf("Stars read and unidentified of Tycho-2: %d\n", TYCunidentified);
    fclose(crossStream);
    fclose(stream);
    return 0;
}
