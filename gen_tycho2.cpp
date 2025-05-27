
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
 * 
 * También hay una versión alternativa para el hemisferio Sur (sólo declinaciones -22
 * a -31) donde se generan tres archivos CSV, el primero con las estrellas de CD simples,
 * el segundo con las CD dobles y el tercero con las de color. 
 * 
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "read_cpd.h"
#include "read_ppm.h"
#include "trig.h"
#include "misc.h"

/* Para uso de la libreria WCS: */
#define WCS_J2000 1 /* J2000(FK5) right ascension and declination */
#define WCS_B1950 2 /* B1950(FK4) right ascension and declination */
extern "C" void wcsconp(int sys1, int sys2, double eq1, double eq2, double ep1, double ep2,
             double *dtheta, double *dphi, double *ptheta, double *pphi);

#define STRING_SIZE 14
#define MAXUNSTAR 1000
#define THRESHOLD_PPM 15.0
#define THRESHOLD_CPD 30.0
#define THRESHOLD_CD 45.0
#define DIST_PPM_TYC 15.0
#define DIST_CPD_TYC 30.0 
#define DIST_DM_TYC 60.0
#define DIST_UN_TYC 15.0

// Nombres alternativos para BD/CD y CPD, también atributos de CD
char dmStringDM[MAXDMSTAR][STRING_SIZE];
bool dmIsColor[MAXDMSTAR];
bool dmIsDouble[MAXDMSTAR];
char dmStringCPD[MAXCPDSTAR][STRING_SIZE];

// Datos para estrellas no identificadas
char unidentifiedName[MAXUNSTAR][STRING_SIZE];
double unidentifiedX[MAXUNSTAR];
double unidentifiedY[MAXUNSTAR];
double unidentifiedZ[MAXUNSTAR];
bool alsoUnidentifiedFromTYC[MAXUNSTAR];
int countUnidentified = 0; 

/*
 * readCrossFile - lee archivos de identificaciones cruzadas en formato CSV
 */
void readCrossFile(
        const char *ppm_file, struct PPMstar_struct *PPMstar, int PPMstars,
        const char *cd_file, struct DMstar_struct *DMstar, int DMstars,
        const char *cpd_file, struct CPDstar_struct *CPDstar, int CPDstars) {
    char buffer[1024];
    char targetRef[STRING_SIZE];
    int ppmRef, cdDeclRef, cdNumRef, cpdDeclRef, cpdNumRef;
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
        if (minDistance < __FLT_EPSILON__ || minDistance > THRESHOLD_PPM) {
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

    /* read cross file between CPD and target (optative) */
    if (CPDstars != 0) {
        printf("Reading cross file %s... ", cpd_file);
        stream = fopen(cpd_file, "rt");
        if (stream == NULL) {
            snprintf(buffer, 1024, "Cannot read %s", cpd_file);
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
            sscanf(buffer, "%13[^,],CPD %d°%d,%f\n", targetRef, &cpdDeclRef, &cpdNumRef, &minDistance);
            // printf("Cross: %s, CPD %d°%d, dist = %.1f arcsec.\n", targetRef, cpdDeclRef, cpdNumRef, minDistance);
            if (minDistance < __FLT_EPSILON__ || minDistance > THRESHOLD_CPD) {
                // omit identifications with zero distance (bug) or too far away
                continue;
            }
            for (int i = 0; i < DMstars; i++) {
                if (CPDstar[i].declRef == cpdDeclRef && CPDstar[i].numRef == cpdNumRef) {
                    strncpy(dmStringCPD[i], targetRef, STRING_SIZE);
                    break;
                }
            }
        }
        fclose(stream);
        printf("done!\n");
    }

    /* read cross file between CD and target (optative) */
    if (DMstars != 0) {
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
            if (minDistance < __FLT_EPSILON__ || minDistance > THRESHOLD_CD) {
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
}

/*
 * readUnidentifiedFile - lee archivos de estrellas no identificadas
 */
void readUnidentifiedFile(const char *file) {
    char buffer[1024];
    char targetRef[STRING_SIZE];
    double x, y, z;

    printf("Reading file %s... ", file);
    FILE *stream = fopen(file, "rt");
    if (stream == NULL) {
        snprintf(buffer, 1024, "Cannot read %s", file);
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
        sscanf(buffer, "%13[^,],%lf,%lf,%lf\n", targetRef, &x, &y, &z);
        // printf("Unidentified: %s -> (%.12f, %.12f, %.12f)\n", targetRef, x, y, z);

        if (countUnidentified >= MAXUNSTAR) {
            printf("Error: too many unidentified stars.\n");
            exit(1);
        }
        strncpy(unidentifiedName[countUnidentified], targetRef, STRING_SIZE);
        unidentifiedX[countUnidentified] = x;
        unidentifiedY[countUnidentified] = y;
        unidentifiedZ[countUnidentified] = z;
        alsoUnidentifiedFromTYC[countUnidentified] = true;
        countUnidentified++;
    }
    fclose(stream);
    printf("done!\n");
}

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv) {
    char buffer[1024];
    char cell[256];

    /* leemos catalogo BD/CD */
#ifdef ALTERNATIVE
    readDM("cat/cd_vol1_curated.txt");
#else    
    readDM(isCD() ? "cat/cd_curated.txt" : "cat/bd_curated.txt");
#endif
    struct DMstar_struct *DMstar = getDMStruct();
    int DMstars = getDMStars();
    bool is_north = !isCD();

    /* generamos identificadores de DM */
    for (int i = 0; i < DMstars; i++) {
#ifdef ALTERNATIVE
        const char *dmCat = "";
#else       
        const char *dmCat = isCD() ? "CD " : "BD ";
#endif
        snprintf(dmStringDM[i], STRING_SIZE, "%s%c%d°%d",
                dmCat,
                DMstar[i].signRef ? '-' : '+',
                abs(DMstar[i].declRef),
                DMstar[i].numRef);
        dmIsColor[i] = false;
        dmIsDouble[i] = false;
    }

    struct CPDstar_struct *CPDstar = nullptr;
    int CPDstars = 0;
    struct PPMstar_struct *PPMstar = nullptr;
    int PPMstars = 0;
#ifdef ALTERNATIVE
    /* leemos atributos de color y dobles */
    for (int decl = 22; decl <= 24; decl++) {
        int cdRef;
        // Color
        snprintf(buffer, 1024, "cd/color_%d.txt", decl);
        FILE *stream = fopen(buffer, "rt");
        if (stream == NULL) {
            snprintf(buffer, 1024, "Cannot read %s", buffer);
            perror(buffer);
            exit(1);
        }
        while (fgets(buffer, 1023, stream) != NULL) {
            sscanf(buffer, "%d\n", &cdRef);
            if (cdRef == 0) continue;
            int index = getDMindex(true, decl, cdRef);
            dmIsColor[index] = true;
        }
        fclose(stream);
        // Dobles
        snprintf(buffer, 1024, "cd/dpl_%d.txt", decl);
        stream = fopen(buffer, "rt");
        if (stream == NULL) {
            snprintf(buffer, 1024, "Cannot read %s", buffer);
            perror(buffer);
            exit(1);
        }
        while (fgets(buffer, 1023, stream) != NULL) {
            sscanf(buffer, "%d\n", &cdRef);
            if (cdRef == 0) continue;
            int index = getDMindex(true, decl, cdRef);
            dmIsDouble[index] = true;
        }
        fclose(stream);
    }
#else
    /* leemos catalogo CPD (en caso de ubicarnos en el Sur) */
    if (!is_north) {
        readCPD(false, false);
        CPDstar = getCPDStruct();
        CPDstars = getCPDStars();

        /* generamos identificadores de CPD */
        for (int i = 0; i < CPDstars; i++) {
            snprintf(dmStringCPD[i], STRING_SIZE, "CPD %d°%d",
                CPDstar[i].declRef,
                CPDstar[i].numRef);
        }
    }

    /* leemos catalogo PPM (2000 en B1950) */
    readPPM(false, true, false, false, 2000.0);
    sortPPM();
    PPMstar = getPPMStruct();
    PPMstars = getPPMStars();

    if (isCD()) {
        /* leemos identificaciones cruzadas y renombramos designaciones.
         * Order de prioridad: UA (solo PPM), GC, ZC, OA, U, G, G2 */
        readCrossFile(
            "results/cross/cross_gc2_ppm.csv", PPMstar, PPMstars,
            "results/cross/cross_gc2_cd.csv", DMstar, DMstars,
            "results/cross/cross_gc2_cpd.csv", CPDstar, CPDstars);
        readCrossFile(
            "results/cross/cross_gilliss_ppm.csv", PPMstar, PPMstars,
            "results/cross/cross_gilliss_cd.csv", DMstar, DMstars,
            "results/cross/cross_gilliss_cpd.csv", CPDstar, CPDstars);
        readCrossFile(
            "results/cross/cross_usno_ppm.csv", PPMstar, PPMstars,
            "results/cross/cross_usno_cd.csv", DMstar, DMstars,
            "results/cross/cross_usno_cpd.csv", CPDstar, CPDstars);
        readCrossFile(
            "results/cross/cross_oa_ppm.csv", PPMstar, PPMstars,
            "results/cross/cross_oa_cd.csv", DMstar, DMstars,
            "results/cross/cross_oa_cpd.csv", CPDstar, CPDstars);
        readCrossFile(
            "results/cross/cross_zc_ppm.csv", PPMstar, PPMstars,
            "results/cross/cross_zc_cd.csv", DMstar, DMstars,
            "results/cross/cross_zc_cpd.csv", CPDstar, CPDstars);
        readCrossFile(
            "results/cross/cross_gc_ppm.csv", PPMstar, PPMstars,
            "results/cross/cross_gc_cd.csv", DMstar, DMstars,
            "results/cross/cross_gc_cpd.csv", CPDstar, CPDstars);
        readCrossFile(
            "results/cross/cross_ua_ppm.csv", PPMstar, PPMstars,
            "", nullptr, 0,
            "", nullptr, 0);
        /* también leemos las no identificadas */
        readUnidentifiedFile("results/cross/gc_unidentified.csv");
        readUnidentifiedFile("results/cross/zc_unidentified.csv");
        readUnidentifiedFile("results/cross/oa_unidentified.csv");
        readUnidentifiedFile("results/cross/usno_unidentified.csv");
        readUnidentifiedFile("results/cross/gilliss_unidentified.csv");
        readUnidentifiedFile("results/cross/gc2_unidentified.csv");
    } else {
        /* leemos identificaciones cruzadas y renombramos designaciones,
         * Solo UA y USNO contra PPM. */
        readCrossFile(
            "results/cross/cross_usno_ppm.csv", PPMstar, PPMstars,
            "", nullptr, 0,
            "", nullptr, 0);
        readCrossFile(
            "results/cross/cross_ua_ppm.csv", PPMstar, PPMstars,
            "", nullptr, 0,
            "", nullptr, 0);
        readUnidentifiedFile("results/cross/usno_unidentified.csv");
    }
#endif

    FILE *stream = fopen("cat/tyc2.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read tyc2.txt");
        exit(1);
    }

    FILE *stream2 = fopen("cat/tyc2_suppl.txt", "rt");
    if (stream2 == NULL) {
        perror("Cannot read tyc2.txt");
        exit(1);
    }

#ifdef ALTERNATIVE
    FILE *crossStream = openCrossFile("tycho2/cross_tyc2_south_plain.csv");
    FILE *crossStreamColor = openCrossFile("tycho2/cross_tyc2_south_color.csv");
    FILE *crossStreamDpl = openCrossFile("tycho2/cross_tyc2_south_dpl.csv");
#else
    snprintf(buffer, 1024, "tycho2/cross_tyc2_%s.csv", is_north ? "north" : "south");
    FILE *crossStream = openCrossFile(buffer);
#endif

    int TYCstarsPPM = 0;
    int TYCstarsDM = 0;
    int TYCstarsCPD = 0;
    int TYCstarsOther = 0;
    int TYCunidentified = 0;
    int entry = 0;

    printf("Starting with TYC supplementary catalog...\n");
    bool readSupplement = true;
    for (;;) {
        if (readSupplement) {
            /* lee del catálogo suplemento hasta consumirlo */
            if (fgets(buffer, 1023, stream2) == NULL) {
                printf("Now reading main TYC catalog...\n");
                readSupplement = false;
                continue;
            }
        } else {
            /* lee del catálogo principal */
            if (fgets(buffer, 1023, stream) == NULL) break;
        }
        // if (TYCstarsPPM > 100) break;
		entry++;
		if (entry % 10000 == 0) {
			printf("Progress (%.2f%%): PPM = %d, DM = %d, CPD = %d, other = %d, unidentified = %d\n",
                (100.0 * (float) entry) / 2539913.0,
                TYCstarsPPM,
                TYCstarsDM,
                TYCstarsCPD,
                TYCstarsOther,
                TYCunidentified);
		}

        /* lee declinación y descarta tempranamente */
        readField(buffer, cell, readSupplement ? 29 : 166, 12);
        double Decl = atof(cell);
        if (is_north) {
            if (Decl < 0) continue;
        } else {
            if (Decl > 0) continue;
        }

        /* lee numeración */
        readField(buffer, cell, 1, 4);
        int tyc1Ref = atoi(cell);
        readField(buffer, cell, 6, 5);
        int tyc2Ref = atoi(cell);
        readField(buffer, cell, 12, 1);
        int tyc3Ref = atoi(cell);

        char tycString[20];
        snprintf(tycString, 20, "TYC %d-%d-%d", tyc1Ref, tyc2Ref, tyc3Ref);

        /* lee RA y Decl (epoch) */
        readField(buffer, cell, readSupplement ? 16 : 153, 12);
        double RA = atof(cell);
        double epRA, epDecl, epoch;
        if (readSupplement) {
            epoch = 1991.25;
            epRA = epoch;
            epDecl = epoch;
        } else {
            readField(buffer, cell, 179, 4);
            epRA = atof(cell);
            readField(buffer, cell, 184, 4);
            epDecl = atof(cell);
            /* la época es un promedio de las de RA y Decl */
            epoch = 1990.0 + (epRA + epDecl) / 2.0;
        }
        /* lee mov. propios, si los tiene */
        double pmRA = 0.0;
        double pmDecl = 0.0;
        readField(buffer, cell, 14, 1);
        if (cell[0] != readSupplement ? 'T' : 'X') {
            readField(buffer, cell, 42, 7);
            pmRA = atof(cell);
            readField(buffer, cell, 50, 7);
            pmDecl = atof(cell);
            
            pmRA /= 1000 * 3600 * dcos(epDecl); /* conversion de mas/yr a grados/yr (juliano) */
            pmDecl /= 1000 * 3600; /* conversion de mas/yr a grados/yr (juliano) */
        }

        // printf("TYC %d-%d-%d: RA = %f, Decl = %f, pmRA = %f, pmDecl = %f, epoch = %f\n", tyc1Ref, tyc2Ref, tyc3Ref, RA, Decl, pmRA, pmDecl, epoch);

        /* convertir a 2000.0 (B1950) ya que las PPM fueron también convertidas ahí */
        double RAtarget = RA;
        double Decltarget = Decl;
        double pmRAtarget = pmRA;
        double pmDecltarget = pmDecl;
        wcsconp(WCS_J2000, WCS_B1950, 0.0, 2000.0, epoch + 0.001278, 2000.0, &RAtarget, &Decltarget, &pmRAtarget, &pmDecltarget);
        // printf("    RA = %f, Decl = %f, pmRA = %f, pmDecl = %f, epoch = %f\n", RAtarget, Decltarget, pmRAtarget, pmDecltarget, 2000.0);

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RAtarget, Decltarget, &x, &y, &z);

        /* halla PPM más cercana, dentro de 15 arcsec */
        int ppmIndex = -1;
        double minDistance = HUGE_NUMBER;
        findPPMByCoordinates(x, y, z, Decltarget, &ppmIndex, &minDistance);
        if (ppmIndex != -1 && minDistance < DIST_PPM_TYC) {
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
        if (is_north || (Decltarget <= -22.0
#ifdef ALTERNATIVE
            && Decltarget > -32.0
#endif
        )) {
            int dmIndex = -1;
            minDistance = HUGE_NUMBER;
            findDMByCoordinates(x, y, z, Decltarget, &dmIndex, &minDistance);
            if (dmIndex != -1 && minDistance < DIST_DM_TYC) {
                /* se almacena la identificación cruzada con la DM */
                FILE *chosenStream = crossStream;
#ifdef ALTERNATIVE
                if (dmIsDouble[dmIndex]) {
                    chosenStream = crossStreamDpl;
                } else if (dmIsColor[dmIndex]) {
                    chosenStream = crossStreamColor;
                }
#endif
                writeCrossEntry(chosenStream, tycString, dmStringDM[dmIndex], minDistance);
                TYCstarsDM++;
                continue;
            }
        }

        /* halla la CPD más cercana, dentro de 30 arcsec */
        if (Decltarget <= -18.0) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            findCPDByCoordinates(x, y, z, Decltarget, &cpdIndex, &minDistance);
            if (cpdIndex != -1 && minDistance < DIST_CPD_TYC) {
                /* se almacena la identificación cruzada con la CPD */
                writeCrossEntry(crossStream, tycString, dmStringCPD[cpdIndex], minDistance);
                TYCstarsCPD++;
                continue;
            }
        }

        /* Escanea las estrellas que no fueron identificadas con PPM/CD/CPD.
         * Para poder priorizar los catálogos, si se encuentra una más cerca,
         * debe sobrepasar un threshold de 1 arco de segundo más para elegirse. */
        int index = -1;
        minDistance = HUGE_NUMBER;
        for (int i = 0; i < countUnidentified; i++) {
            double dist = 3600.0 * calcAngularDistance(x, y, z, unidentifiedX[i], unidentifiedY[i], unidentifiedZ[i]);
            if (dist > DIST_UN_TYC) continue;
            alsoUnidentifiedFromTYC[i] = false;
            if (minDistance - 1.0 > dist) {
                index = i;
                minDistance = dist;
            }
        }
        if (index != -1) {
            /* se almacena la identificación cruzada con la estrella */
            writeCrossEntry(crossStream, tycString, unidentifiedName[index], minDistance);
            TYCstarsOther++;
            continue;
        }

        TYCunidentified++;
    }
    printf("\nStars read and identified of Tycho-2 from PPM: %d\n", TYCstarsPPM);
    printf("Stars read and identified of Tycho-2 from DM: %d\n", TYCstarsDM);
    printf("Stars read and identified of Tycho-2 from CPD: %d\n", TYCstarsCPD);
    printf("Stars read and identified of Tycho-2 from other catalogues: %d\n", TYCstarsOther);
    printf("Stars read and remain unidentified of Tycho-2: %d\n", TYCunidentified);
    fclose(crossStream);
#ifdef ALTERNATIVE
    fclose(crossStreamColor);
    fclose(crossStreamDpl);
#endif
    fclose(stream2);
    fclose(stream);

    printf("\nStars from other catalogues yet not identified:\n");
    for (int i = 0; i < countUnidentified; i++) {
        if (!alsoUnidentifiedFromTYC[i]) continue;

        double RA, Decl;
        rec2sph(unidentifiedX[i], unidentifiedY[i], unidentifiedZ[i], &RA, &Decl);
        printf("  %s in %s hemisphere (decl = %.2f°)\n", unidentifiedName[i],
            Decl >= 0 ? "northern" : "southern", Decl);
    }
    return 0;
}
