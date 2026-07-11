
/*
 * CROSS_SOUTH - Compara registros de varios catalogos
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "read_ppm.h"
#include "read_cpd.h"
#include "read_dm.h"
#include "read_gc.h"
#include "trig.h"
#include "misc.h"
#include "find_gsc.h"
#include "cross_utils.h"

#define MAXOASTAR 19000
#define MAXLALSTAR 47400
#define MAXLACSTAR 10000
#define MAXSTSTAR 12500
#define MAXBRISTAR 8000
#define MAXTAYLORSTAR 11100
#define MAXUSNOSTAR 11000
#define MAXGILLISSSTAR 17000
#define MAXZCSTAR 15000
#define MAXCLSTAR 100
#define EPOCH_GC2 1900.0
#define EPOCH_OA 1850.0
#define EPOCH_LAL 1800.0
#define EPOCH_ST 1880.0
#define EPOCH_BRI 1825.0
#define EPOCH_TAYLOR 1835.0
#define EPOCH_USNO 1860.0
#define EPOCH_UA 1875.0
#define EPOCH_GILLISS 1850.0
#define EPOCH_CL 1872.0
#define MAX_DIST_OA_PPM 45.0
#define MAX_DIST_LAL_PPM 90.0
#define MAX_DIST_BRI_PPM 90.0
#define MAX_DIST_TAY_PPM 40.0
#define MAX_DIST_OA_CPD 45.0
#define MAX_DIST_ST_CPD 30.0
#define MAX_DIST_USNO_PPM 30.0
#define MAX_DIST_USNO_CPD 30.0
#define MAX_DIST_TH_CPD 30.0
#define MAX_DIST_UA_PPM 45.0
#define MAX_DIST_ZC_ZC 15.0
#define CURATED true // true if curated CD catalog should be used

// Here, we save 1875.0 coordinates of OA stars in rectangular form
double oaX[MAXOASTAR], oaY[MAXOASTAR], oaZ[MAXOASTAR], oaMag[MAXOASTAR];
int oaRef[MAXOASTAR];
int countOA = 0;
StarList oaList = {&countOA, oaRef, oaX, oaY, oaZ};

// Here, we save 1875.0 coordinates of Lalande stars in rectangular form
double lalX[MAXLALSTAR], lalY[MAXLALSTAR], lalZ[MAXLALSTAR], lalMag[MAXLALSTAR];
int lalRef[MAXLALSTAR];
int countLal = 0;
StarList lalList = {&countLal, lalRef, lalX, lalY, lalZ};

// Here, we save 1875.0 coordinates of Lacaille stars in rectangular form
double lacX[MAXLACSTAR], lacY[MAXLACSTAR], lacZ[MAXLACSTAR], lacMag[MAXLACSTAR];
int lacRef[MAXLACSTAR];
int countLac = 0;
StarList lacList = {&countLac, lacRef, lacX, lacY, lacZ};

// Here, we save 1875.0 coordinates of Taylor stars in rectangular form
double tayX[MAXTAYLORSTAR], tayY[MAXTAYLORSTAR], tayZ[MAXTAYLORSTAR], tayMag[MAXTAYLORSTAR];
int tayRef[MAXTAYLORSTAR];
int countTaylor = 0;
StarList tayList = {&countTaylor, tayRef, tayX, tayY, tayZ};

// Here, we save 1875.0 coordinates of Stone stars in rectangular form
double stX[MAXSTSTAR], stY[MAXSTSTAR], stZ[MAXSTSTAR];
int stRef[MAXSTSTAR];
int countSt = 0;
StarList stList = {&countSt, stRef, stX, stY, stZ};

// Maps a Taylor designation (number) to the index (into stX/stY/stZ/stRef)
// of the Stone star that cross-references it.
int stTayRef[MAXTAYLORSTAR];

// Here, we save 1875.0 coordinates of Circumpolar List (CL) stars in rectangular form
double clX[MAXCLSTAR], clY[MAXCLSTAR], clZ[MAXCLSTAR];
int clRef[MAXCLSTAR];
int countCL = 0;
StarList clList = {&countCL, clRef, clX, clY, clZ};

// Here, we save 1875.0 coordinates of Brisbane stars in rectangular form
double briX[MAXBRISTAR], briY[MAXBRISTAR], briZ[MAXBRISTAR], briMag[MAXBRISTAR];
int briRef[MAXBRISTAR];
int countBri = 0;
StarList briList = {&countBri, briRef, briX, briY, briZ};

// Here, we save 1875.0 coordinates of USNO stars in rectangular form
double usnoX[MAXUSNOSTAR], usnoY[MAXUSNOSTAR], usnoZ[MAXUSNOSTAR], usnoMag[MAXUSNOSTAR];
int usnoRef[MAXUSNOSTAR];
int countUsno = 0;
StarList usnoList = {&countUsno, usnoRef, usnoX, usnoY, usnoZ};

// Here, we save 1875.0 coordinates of Gilliss stars in rectangular form
double gilX[MAXGILLISSSTAR], gilY[MAXGILLISSSTAR], gilZ[MAXGILLISSSTAR], gilMag[MAXGILLISSSTAR];
int gilRef[MAXGILLISSSTAR];
int countGil = 0;

// Here, we save 1875.0 coordinates of Gould's ZC stars in rectangular form
double zcX[MAXZCSTAR], zcY[MAXZCSTAR], zcZ[MAXZCSTAR];
int zcHour[MAXZCSTAR], zcNum[MAXZCSTAR];
int countZC = 0;
int countPPMZC = 0, countGSCZC = 0, countCDZC = 0, countCPDZC = 0;
FILE *crossPPMZCStream;
FILE *crossCDZCStream;
FILE *crossCPDZCStream;
FILE *crossSAOZCStream;
FILE *crossHDZCStream;
FILE *unidentifiedZCStream;

void saveZC(int RAh, int RAs, int Decls,
        int numRef, double x, double y, double z,
        bool ppmFound, int ppmIndex, double minPPMDistance,
        bool gscFound, const char *gscId, double nearestGSCDistance,
        bool cdFound, int cdIndex, double minCDDistance,
        bool cpdFound, int cpdIndex, double minCPDDistance) {
    char zcName[20], catName[20];

    /* check if the star was previously registered */
    for (int i = 0; i < countZC; i++) {
        if (zcNum[i] != numRef || zcHour[i] != RAh) continue;
        double dist = 3600.0 * calcAngularDistance(x, y, z, zcX[i], zcY[i], zcZ[i]);
        if (dist > MAX_DIST_ZC_ZC) {
            printf("**) Warning: ZC %dh %d is FAR from previous registration (dist = %.1f arcsec).\n",
                RAh,
                numRef,
                dist);
        }
        return;
    }

    /* save cross identifications */
    snprintf(zcName, 20, "ZC %dh %d", RAh, numRef);
    if (ppmFound) {
        countPPMZC++;
        struct PPMstar_struct *PPMstar = getPPMStruct();
        writePPMCrossEntry(crossPPMZCStream, crossSAOZCStream, crossHDZCStream, zcName, &PPMstar[ppmIndex], 0.0, minPPMDistance);
    }
    if (gscFound) {
        countGSCZC++;
        snprintf(catName, 20, "GSC %s", gscId);
        writeCrossEntry(crossPPMZCStream, zcName, catName, 0.0, nearestGSCDistance);
    }
    if (cdFound) {
        countCDZC++;
        struct DMstar_struct *CDstar = getDMStruct();
        snprintf(catName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
        writeCrossEntry(crossCDZCStream, zcName, catName, 0.0, minCDDistance);
    }
    if (cpdFound) {
        countCPDZC++;
        struct CPDstar_struct *CPDstar = getCPDStruct();
        snprintf(catName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
        writeCrossEntry(crossCPDZCStream, zcName, catName, 0.0, minCPDDistance);
    }

    if (!ppmFound && !cdFound && !cpdFound) {
        if (gscFound) {
            writeUnidentified(unidentifiedZCStream, zcName, x, y, z);
        } else {
            printf("**) %s is also ALONE.\n", zcName);
            logCauses(zcName, true,
                false, false, RAs, -50.0, Decls, -1, 0.0);
        }
    }

    /* add the new star */
    int i = storeStar(&countZC, MAXZCSTAR, "ZC", zcNum, zcX, zcY, zcZ, NULL, numRef, x, y, z, 0.0);
    zcHour[i] = RAh;
}

/*
 * readGC2 - lee y cruza segundo catálogo argentino
 * (ya se deben haber leidos los catalogs CD y CPD)
 */
void readGC2() {
    char buffer[1024], cell[256], catName[20];

    printf("\n***************************************\n");
    printf("Perform comparison between GC2 and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_gc2_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_gc2_cpd.csv");
    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("gc2", &crossPPMStream, &crossSAOStream, &crossHDStream);

    CrossStats stats;
    int countGC = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(EPOCH_GC2, false);

    /* leemos catalogo */
    FILE *stream = fopen("cat/gc2.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read gc2.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		/* descarta otras obs que no sean la primera */
		readField(buffer, cell, 12, 2);
		if (atoi(cell) > 0) continue;

		/* lee numeración */
		readField(buffer, cell, 8, 4);
		int gcRef = atoi(cell);

        snprintf(catName, 20, "G2 %d", gcRef);

        /* descarta declinacion positiva */
		readField(buffer, cell, 44, 1);
        if (cell[0] != '-') {
            printf("Note: %s discarded due to positive declination.\n", catName);
            continue;
        }

		/* lee ascension recta B1900.0 */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 21, &RAh, &RAm, &RAs);

		/* lee declinacion B1900.0 (incorpora signo negativo; en nuestro caso, siempre) */
        int Decld, Declm, Decls;
        double Decl = -readDeclField(buffer, 45, 2, 47, 49, 3, 10.0, &Decld, &Declm, &Decls);

		/* lee magnitud */
		readField(buffer, cell, 14, 2);
		float vmag = atof(cell);
		readFieldSanitized(buffer, cell, 16, 1);
        vmag += atof(cell)/10.0;

		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z, minDistance;
		int ppmIndex;
		sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_PPM, catName, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;

        /* si no encuentra una PPM cercana, prueba con GSC */
        bool gscFound = tryGSC(ppmFound, RA, Decl, EPOCH_GC2, crossPPMStream, catName, vmag, &stats);

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_GC2, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        /* busca la CPD y la CD mas cercanas y genera los cruzamientos */
        bool cpdFound = false;
        if (Decl1875 <= -18.0)
            cpdFound = crossWithCPD(x, y, z, Decl1875, vmag, MAX_DIST_CPD, crossCPDStream, catName, &stats, NULL, NULL);

        bool cdFound = false;
        if (Decl1875 <= -22.0)
            cdFound = crossWithCD(x, y, z, Decl1875, vmag, catName, crossCDStream, catName, &stats, NULL, NULL);

        if (!ppmFound && !cdFound && !cpdFound && !gscFound) {
            warnAlone(&stats.errors, catName, NULL, NULL,
                catName, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }
        countGC++;
    }
    fclose(stream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available G2 stars = %d\n", countGC);
    printf("Stars from G2 identified with PPM = %d, GSC-PPM = %d, CD = %d and CPD = %d\n", stats.countDist, stats.countGSC, stats.countCD, stats.countCPD);
    printRSMEDist(&stats);
    printRSMEMag(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * readWeiss - lee y cruza catalogo de Weiss
 * (ya se deben haber leidos los catalogs CD y CPD)
 */
void readWeiss() {
    char buffer[1024], cell[256], catName[20], warnName[20], aloneName[32];
    char catLine[64];

    printf("\n***************************************\n");
    printf("Perform comparison between Weiss and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_oa_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_oa_cpd.csv");
    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("oa", &crossPPMStream, &crossSAOStream, &crossHDStream);
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/oa_unidentified.csv");
    FILE *catalogStream = openCatalogFile("likelihood/cat1875/oa.csv");

    CrossStats stats;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(EPOCH_OA, false);

    /* leemos catalogo Weiss (pero solo nos interesa su identificación OA) */
    FILE *stream = fopen("cat/weiss.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read weiss.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		/* descarta otras obs que no sean la primera */
		readField(buffer, cell, 6, 1);
		if (cell[0] != '1') continue;

		/* lee numeración catálogo Weiss y Oeltzen-Argelander */
		readField(buffer, cell, 49, 5);
		int oeltzenRef = atoi(cell);
        if (oeltzenRef == 0) continue;
		readField(buffer, cell, 1, 5);
		int weissRef = atoi(cell);

        snprintf(catName, 20, "OA %d", oeltzenRef);
        snprintf(warnName, 20, "W %d", weissRef);

		/* lee ascension recta B1850.0 */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 13, &RAh, &RAm, &RAs);

		/* lee declinacion B1850.0 (incorpora signo negativo; en nuestro caso, siempre) */
        int Decld, Declm, Decls;
        double Decl = -readDeclField(buffer, 22, 2, 24, 26, 3, 10.0, &Decld, &Declm, &Decls);

        formatCatLine(catLine, RAh, RAm, RAs, '-', Decld, Declm, Decls);

		/* lee magnitud */
		readField(buffer, cell, 8, 2);
		int magInt = atoi(cell);
		float vmag = (float) magInt;
		readField(buffer, cell, 10, 2);
		if (cell[0] != ' ') {
		    int magFrac = atoi(cell);
		    if (magFrac != magInt + 1) {
		        printf("Error: Weiss %d has unexpected value %d in magnitude columns 10-11 (expected %d or blank).\n",
		            weissRef, magFrac, magInt + 1);
		        exit(1);
		    }
		    vmag += 0.5;
		}

		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z, minDistance;
		int ppmIndex;
		sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_OA_PPM, warnName, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;

        /* si no encuentra una PPM cercana, prueba con GSC */
        bool gscFound = tryGSC(ppmFound, RA, Decl, EPOCH_OA, crossPPMStream, catName, vmag, &stats);

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_OA, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        /* busca la CPD y la CD mas cercanas y genera los cruzamientos */
        bool cpdFound = false;
        if (Decl1875 <= -18.0)
            cpdFound = crossWithCPD(x, y, z, Decl1875, vmag, MAX_DIST_OA_CPD, crossCPDStream, catName, &stats, NULL, NULL);

        bool cdFound = false;
        if (Decl1875 <= -22.0)
            cdFound = crossWithCD(x, y, z, Decl1875, vmag, warnName, crossCDStream, catName, &stats, NULL, NULL);

        if (!ppmFound && !cdFound && !cpdFound) {
            if (gscFound) {
                writeUnidentified(unidentifiedStream, catName, x, y, z);
            } else {
                snprintf(aloneName, 32, "W %d (OA %d)", weissRef, oeltzenRef);
                warnAlone(&stats.errors, aloneName, warnName, catLine,
                    catName, RAs, Decl1875, Decls,
                    PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
            }
        }

        /* la almacenamos para futuras identificaciones */
        storeStar(&countOA, MAXOASTAR, "OA", oaRef, oaX, oaY, oaZ, oaMag, oeltzenRef, x, y, z, vmag);
        writeCatalogFile(catalogStream, catName, x, y, z, vmag);
    }
    fclose(stream);
    fclose(unidentifiedStream);
    fclose(catalogStream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available OA stars = %d\n", countOA);
    printf("Stars from OA identified with PPM = %d, GSC-PPM = %d, CD = %d and CPD = %d\n", stats.countDist, stats.countGSC, stats.countCD, stats.countCPD);
    printRSMEDist(&stats);
    printRSMEMag(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * readLalande - lee catalogo de Lalande
 */
void readLalande() {
    char buffer[1024], cell[256], catName[20];

    printf("\n***************************************\n");
    printf("Perform comparison between Lalande and PPM...\n");

    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("lalande", &crossPPMStream, &crossSAOStream, &crossHDStream);

    CrossStats stats;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(EPOCH_LAL, false);

    /* leemos catalogo Lalande */
    FILE *stream = fopen("cat/lalande.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read lalande.txt");
		exit(1);
    }

    while (fgets(buffer, 1023, stream) != NULL) {
        /* descarta la estrella si no hay RA o declinación */
        readField(buffer, cell, 19, 1);
        if (cell[0] == ' ') continue;
        readField(buffer, cell, 32, 1);
        if (cell[0] == ' ') continue;

		/* lee numeración */
		readField(buffer, cell, 1, 5);
		int catRef = atoi(cell);
        if (catRef > 0) snprintf(catName, 20, "Lal %d", catRef);

		/* lee magnitud */
		float vmag = readMagIntHalf(buffer, 15, 2, 17);

		/* lee ascension recta B1800.0 (si existiese) */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 19, &RAh, &RAm, &RAs);

		/* lee declinacion B1800.0 (en realidad NPD) */
        int Decld, Declm, Decls;
        double Decl = 90.0 - readDeclField(buffer, 32, 3, 35, 37, 3, 10.0, &Decld, &Declm, &Decls);

		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z, minDistance;
		int ppmIndex;
		sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_LAL_PPM, NULL,
            catRef > 0 ? catName : NULL,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;

        /* si no encuentra una PPM cercana, prueba con GSC */
        bool gscFound = tryGSC(ppmFound, RA, Decl, EPOCH_LAL,
            catRef > 0 ? crossPPMStream : NULL, catName, vmag, &stats);

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_LAL, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        if (!ppmFound && !gscFound && nearestPPMDistance > MAX_DIST_PPM_FAR) {
            char aloneName[20];
            snprintf(aloneName, 20, "Lal %d", catRef);
            warnAlonePPMGSC(&stats.errors, aloneName, aloneName,
                RAs, Decl1875, Decls, PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* la almacenamos para futuras identificaciones */
        storeStar(&countLal, MAXLALSTAR, "Lalande", lalRef, lalX, lalY, lalZ, lalMag, catRef, x, y, z, vmag);
    }
    fclose(stream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);

    printf("Available Lalande stars = %d\n", countLal);
    printf("Stars from Lalande identified with PPM = %d, GSC-PPM = %d\n", stats.countDist, stats.countGSC);
    printRSMEDist(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * readSOM - lee los "Standards of Magnitude" de la Uranometria Argentina
 * (cat/UA_standards.csv, epoca 1875.0) y revisa las referencias al
 * catalogo de Lalande: posicion y magnitud (tolerancia 0.25)
 * (ya se debe haber leido el catalogo Lalande; sus coordenadas
 * almacenadas ya estan en 1875.0, no hace falta transformar)
 */
void readSOM() {
    char buffer[1024], catName[20];
    char catLine[64];

    printf("\n***************************************\n");
    printf("Check references between UA Standards of Magnitude and Lalande...\n");

    int checkLal = 0;
    int errors = 0;

    /* leemos los Standards of Magnitude */
    FILE *stream = fopen("cat/UA_standards.csv", "rt");
    if (stream == NULL) {
        perror("Cannot read UA_standards.csv");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
        char field[17][32];
        double RA, Decl;
        if (!parseSOMLine(buffer, field, catName, catLine, &RA, &Decl)) continue;

        /* omite si no hay referencia a Lalande, o si no es
           confiable (parentesis) */
        if (field[10][0] == 0 || field[10][0] == '(') continue;
        int numRefCat = atoi(field[10]);

        /* calcula rectangulares (Lalande ya esta en 1875.0) */
        double x, y, z;
        sph2rec(RA, Decl, &x, &y, &z);

        for (int i = 0; i < countLal; i++) {
            if (lalRef[i] != numRefCat) continue;
            double dist = 3600.0 * calcAngularDistance(x, y, z, lalX[i], lalY[i], lalZ[i]);
            if (dist > MAX_DIST_CROSS) {
                printf("%d) Warning: %s is FAR from Lal %d (dist = %.1f arcsec).\n",
                    ++errors,
                    catName,
                    numRefCat,
                    dist);
                printf("     Register %s: %s\n", catName, catLine);
            } else checkLal++;

            /* revisa la magnitud, si ambas estan disponibles */
            if (field[11][0] != 0 && lalMag[i] > 0.0) {
                double somMag = atof(field[11]);
                if (fabs(somMag - lalMag[i]) > 0.5) {
                    printf("%d) Warning: %s reports Lal mag %.1f but Lal %d has mag %.1f.\n",
                        ++errors,
                        catName,
                        somMag,
                        numRefCat,
                        lalMag[i]);
                }
            }
        }
    }
	fclose(stream);

    printf("SOM stars properly identified with Lalande = %d\n", checkLal);
    printf("Errors logged = %d\n", errors);
}

/*
 * readStone - lee y cruza catalogo de Stone
 * (ya se deben haber leidos los catalogs CD, CPD y Brisbane)
 */
void readStone() {
    char buffer[1024], cell[256], catName[20], cdName[20], warnName[20], lacName[20];

    /* usamos catalogo CD */
    struct DMstar_struct *CDstar = getDMStruct();

    printf("\n***************************************\n");
    printf("Perform comparison between Stone and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_lacaille_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_lacaille_cpd.csv");
    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("lacaille", &crossPPMStream, &crossSAOStream, &crossHDStream);
    FILE *crossStonePPMStream, *crossStoneSAOStream, *crossStoneHDStream;
    openCrossSet("stone", &crossStonePPMStream, &crossStoneSAOStream, &crossStoneHDStream);
    FILE *catalogStream = openCatalogFile("likelihood/cat1875/lacaille.csv");

    CrossStats stats;
    int checkBri = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(EPOCH_ST, false);

    /* inicializamos el mapa Taylor->Stone (sin referencia) */
    for (int i = 0; i < MAXTAYLORSTAR; i++) stTayRef[i] = -1;

    /* leemos catalogo Stone */
    FILE *stream = fopen("cat/stone1.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read stone1.txt");
		exit(1);
    }
	FILE *stream2 = fopen("cat/stone2.txt", "rt");
    if (stream2 == NULL) {
        perror("Cannot read stone2.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		char buffer2[1024];
		fgets(buffer2, 1023, stream2);

		/* lee numeración catálogo Stone y Lacaille */
		readField(buffer, cell, 8, 5);
		int stoneRef = atoi(cell);
		readField(buffer, cell, 14, 4);
		int lacailleRef = cell[0] == '&' ? 0 : atoi(cell);

        snprintf(warnName, 20, "St %d", stoneRef);
        if (lacailleRef > 0) snprintf(lacName, 20, "L %d", lacailleRef);

		/* lee magnitud */
		readField(buffer, cell, 24, 1);
		float vmag = atof(cell);
		readField(buffer, cell, 25, 1);
		vmag += cell[0] == '&' ? 0 : atof(cell)/10.0;

		/* lee ascension recta B1880.0 */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 33, &RAh, &RAm, &RAs);

		/* lee declinacion B1880.0 (en realidad NPD) */
        int Decld, Declm, Decls;
        double Decl = 90.0 - readDeclField(buffer2, 21, 3, 24, 26, 4, 100.0, &Decld, &Declm, &Decls);

		/* busca la PPM mas cercana y genera los cruzamientos (Lacaille y Stone) */
		double x, y, z, minDistance;
		int ppmIndex;
		sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_PPM, warnName, NULL,
            NULL, NULL, NULL, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;
        if (ppmFound) {
            if (lacailleRef > 0) {
			    writePPMCrossEntry(crossPPMStream, crossSAOStream, crossHDStream, lacName, &PPMstar[ppmIndex], 0.0, minDistance);
            }
            snprintf(catName, 20, "St %d", stoneRef);
            writePPMCrossEntry(crossStonePPMStream, crossStoneSAOStream, crossStoneHDStream, catName, &PPMstar[ppmIndex], vmag, minDistance);
		}

        /* si no encuentra una PPM cercana, prueba con GSC */
        bool gscFound = tryGSC(ppmFound, RA, Decl, EPOCH_ST,
            lacailleRef > 0 ? crossPPMStream : NULL, lacName, 0.0, &stats);

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_ST, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        /* busca la CPD y la CD mas cercanas y genera los cruzamientos */
        bool cpdFound = false;
        if (Decl1875 <= -18.0)
            cpdFound = crossWithCPD(x, y, z, Decl1875, 0.0, MAX_DIST_ST_CPD, crossCPDStream,
                lacailleRef > 0 ? lacName : NULL, &stats, NULL, NULL);

        bool cdFound = false;
        if (Decl1875 <= -22.0) {
            int cdIndex;
            double cdDist;
            cdFound = crossWithCD(x, y, z, Decl1875, vmag, warnName, NULL, NULL, &stats, &cdIndex, &cdDist);
            if (cdFound && lacailleRef > 0) {
	    	    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
		        writeCrossEntry(crossCDStream, lacName, cdName, 0.0, cdDist);
            }
		}

        if (!ppmFound && !cdFound && !cpdFound && !gscFound) {
            char aloneName[32];
            if (lacailleRef > 0) snprintf(aloneName, 32, "St %d (L %d)", stoneRef, lacailleRef);
            else snprintf(aloneName, 32, "St %d", stoneRef);
            warnAlone(&stats.errors, aloneName, NULL, NULL,
                catName, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* revisa identificacion cruzada con Brisbane */
        readField(buffer, cell, 65, 4);
        int brisRefStone = cell[0] == '&' ? 0 : atoi(cell);
        if (brisRefStone > 0) {
            checkCrossRef(warnName, NULL, "B", x, y, z, brisRefStone, &briList, true, -1, &checkBri, &stats.errors);
        }

        /* la almacenamos para futuras identificaciones */
        if (lacailleRef > 0) {
            storeStar(&countLac, MAXLACSTAR, "Lac", lacRef, lacX, lacY, lacZ, lacMag, lacailleRef, x, y, z, vmag);
            writeCatalogFile(catalogStream, lacName, x, y, z, vmag);
        }

        /* guarda la referencia cruzada a Taylor */
        readField(buffer2, cell, 49, 5);
        if (cell[0] != '&') {
            int taylorRefStone = atoi(cell);
            if (taylorRefStone > 0 && taylorRefStone < MAXTAYLORSTAR) {
                stTayRef[taylorRefStone] = countSt;
            }
        }
        storeStar(&countSt, MAXSTSTAR, "Stone", stRef, stX, stY, stZ, NULL, stoneRef, x, y, z, 0.0);
    }
    fclose(stream2);
    fclose(stream);
    fclose(catalogStream);
    closeCrossSet(crossStonePPMStream, crossStoneSAOStream, crossStoneHDStream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available Stone stars = %d, and Lacaille = %d\n", countSt, countLac);
    printf("Stone stars properly identified with Brisbane = %d\n", checkBri);
    printf("Stars from Stone identified with PPM = %d, GSC-PPM = %d, CD = %d and CPD = %d\n", stats.countDist, stats.countGSC, stats.countCD, stats.countCPD);
    printRSMEDist(&stats);
    printRSMEMag(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * readBrisbane - lee catalogo de Brisbane
 */
void readBrisbane() {
    char buffer[1024], cell[256], catName[20];

    printf("\n***************************************\n");
    printf("Perform comparison between Brisbane and PPM...\n");

    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("bri", &crossPPMStream, &crossSAOStream, &crossHDStream);
    FILE *catalogStream = openCatalogFile("likelihood/cat1875/brisbane.csv");

    CrossStats stats;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(EPOCH_BRI, false);

    /* leemos catalogo Brisbane */
    FILE *stream = fopen("cat/brisbane.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read brisbane.txt");
        exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
        /* lee numeración Brisbane */
        readField(buffer, cell, 1, 4);
        int brisRef = atoi(cell);
        snprintf(catName, 20, "B %d", brisRef);

        /* descarta nebulas */
        readField(buffer, cell, 20, 1);
        if (cell[0] == 'N') continue;

        /* descarta si faltan los segundos de RA o SPD */
        readField(buffer, cell, 29, 1);
        if (cell[0] == ' ') continue;
        readField(buffer, cell, 55, 1);
        if (cell[0] == ' ') continue;

        /* lee magnitud (mismo formato que Weiss) */
        readField(buffer, cell, 21, 2);
        float vmag = 0.0;
        if (cell[0] != ' ') {
            int magInt = atoi(cell);
            vmag = (float) magInt;
            readField(buffer, cell, 23, 2);
            if (cell[0] != ' ') {
                int magFrac = atoi(cell);
                if (magFrac == magInt + 1) {
                    vmag += 0.5;
                } else {
                    printf("Warning: B %d has unexpected value %d in magnitude columns 23-24 (expected %d or blank).\n",
                        brisRef, magFrac, magInt + 1);
                }
            }
        }

        /* lee ascension recta B1825.0 */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 25, &RAh, &RAm, &RAs);

        /* lee SPD B1825.0 (distancia polar sur) y la convertimos a declinación */
        int Decld, Declm, Decls;
        double Decl = readDeclField(buffer, 50, 3, 53, 55, 3, 10.0, &Decld, &Declm, &Decls) - 90.0;

		/* busca la PPM mas cercana y genera el cruzamiento */
        double x, y, z, minDistance;
        int ppmIndex;
        sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_BRI_PPM, NULL, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;

        /* si no encuentra una PPM cercana, prueba con GSC */
        bool gscFound = tryGSC(ppmFound, RA, Decl, EPOCH_BRI, crossPPMStream, catName, vmag, &stats);

        /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_BRI, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        if (!ppmFound && !gscFound && nearestPPMDistance > MAX_DIST_PPM_FAR) {
            warnAlonePPMGSC(&stats.errors, catName, catName,
                RAs, Decl1875, Decls, PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* la almacenamos para futuras identificaciones */
        storeStar(&countBri, MAXBRISTAR, "Brisbane", briRef, briX, briY, briZ, briMag, brisRef, x, y, z, vmag);
        writeCatalogFile(catalogStream, catName, x, y, z, vmag);
    }
    fclose(stream);
    fclose(catalogStream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);

    printf("Available Brisbane stars = %d\n", countBri);
    printf("Stars from Brisbane identified with PPM = %d, GSC-PPM = %d\n", stats.countDist, stats.countGSC);
    printRSMEDist(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * readTaylor - lee y cruza catalogo de Taylor
 * (ya se deben haber leidos catalogos Stone, Brisbane y GC)
 */
void readTaylor() {
    char buffer[1024], cell[256], catName[20];

    printf("\n***************************************\n");
    printf("Perform comparison between Taylor and PPM...\n");

    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("taylor", &crossPPMStream, &crossSAOStream, &crossHDStream);
    FILE *catalogStream = openCatalogFile("likelihood/cat1875/taylor.csv");

    CrossStats stats;
    int checkLac = 0;
    int checkBri = 0;
    int checkGC = 0;
    int checkSt = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(EPOCH_TAYLOR, false);

    /* leemos catalogo Taylor */
    FILE *stream = fopen("cat/taylor1.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read taylor1.txt");
		exit(1);
    }
	FILE *stream2 = fopen("cat/taylor2.txt", "rt");
    if (stream2 == NULL) {
        perror("Cannot read taylor2.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		char buffer2[1024];
		fgets(buffer2, 1023, stream2);

		/* lee numeración */
		readField(buffer, cell, 7, 5);
		int taylorRef = atoi(cell);
        snprintf(catName, 20, "T %d", taylorRef);

		/* lee magnitud */
		float vmag = readMagIntHalf(buffer, 42, 2, 44);

		/* lee ascension recta B1835.0 */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 46, &RAh, &RAm, &RAs);

		/* lee declinacion B1835.0 */
        int Decld, Declm, Decls;
        double Decl = readDeclField(buffer2, 8, 2, 10, 12, 4, 100.0, &Decld, &Declm, &Decls);
		readField(buffer2, cell, 7, 1);
		if (cell[0] == '-') Decl = -Decl;

		/* busca la PPM mas cercana y genera el cruzamiento */
        double x, y, z, minDistance;
        int ppmIndex;
        sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_TAY_PPM, NULL, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;

        /* si no encuentra una PPM cercana, prueba con GSC */
        bool gscFound = tryGSC(ppmFound, RA, Decl, EPOCH_TAYLOR, crossPPMStream, catName, vmag, &stats);

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_TAYLOR, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        if (!ppmFound && !gscFound && nearestPPMDistance > MAX_DIST_PPM_FAR) {
            warnAlonePPMGSC(&stats.errors, catName, catName,
                RAs, Decl1875, Decls, PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* lee referencia numerica y referencia a catalogo (que queda en "cell") */
		readField(buffer, cell, 37, 5);
        int numRefCat = atoi(cell);
		readField(buffer, cell, 31, 3);

        if (!strncmp(cell, "LAC", 3))
            checkCrossRef(catName, NULL, "L", x, y, z, numRefCat, &lacList, false, -1, &checkLac, &stats.errors);
        if (!strncmp(cell, "GOU", 3))
            checkCrossRefGCcat(catName, NULL, numRefCat, x, y, z, &checkGC, &stats.errors);
        if (!strncmp(cell, "BRI", 3))
            checkCrossRef(catName, NULL, "B", x, y, z, numRefCat, &briList, true, -1, &checkBri, &stats.errors);

        /* revisa identificacion cruzada con Stone (acceso O(1) via stTayRef) */
        if (taylorRef > 0 && taylorRef < MAXTAYLORSTAR && stTayRef[taylorRef] >= 0) {
            int i = stTayRef[taylorRef];
            double dist = 3600.0 * calcAngularDistance(x, y, z, stX[i], stY[i], stZ[i]);
            if (dist > MAX_DIST_CROSS) {
                printf("%d) Warning: T %d is FAR from St %d (dist = %.1f arcsec).\n",
                    ++stats.errors,
                    taylorRef,
                    stRef[i],
                    dist);
            } else checkSt++;
        }

        /* la almacenamos para futuras identificaciones */
        storeStar(&countTaylor, MAXTAYLORSTAR, "Taylor", tayRef, tayX, tayY, tayZ, tayMag, taylorRef, x, y, z, vmag);
        writeCatalogFile(catalogStream, catName, x, y, z, vmag);
    }
    fclose(stream2);
    fclose(stream);
    fclose(catalogStream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);

    printf("Available Taylor stars = %d\n", countTaylor);
    printf("Taylor stars properly identified with Lacaille = %d\n", checkLac);
    printf("Taylor stars properly identified with Brisbane = %d\n", checkBri);
    printf("Taylor stars properly identified with GC = %d\n", checkGC);
    printf("Taylor stars properly identified with Stone = %d\n", checkSt);
    printf("Stars from Taylor identified with PPM = %d, GSC-PPM = %d\n", stats.countDist, stats.countGSC);
    printRSMEDist(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * readUSNO - lee y cruza catalogo de Yarnall-Frisby
 * (ya se deben haber leidos los catalogs CD, CPD, Weiss y Stone)
 * tambien revisa referencias cruzadas a OA, Lacaille y Brisbane
 */
void readUSNO() {
    char buffer[1024], cell[256], catName[20];
    char catLine[64];

    printf("\n***************************************\n");
    printf("Perform comparison between USNO and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_usno_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_usno_cpd.csv");
    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("usno", &crossPPMStream, &crossSAOStream, &crossHDStream);
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/usno_unidentified.csv");
    FILE *catalogStream = openCatalogFile("likelihood/cat1875/usno.csv");

    CrossStats stats;
    int checkLac = 0;
    int checkTaylor = 0;
    int checkLal = 0;
    int checkOA = 0;
    int checkBri = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(EPOCH_USNO, false);

    /* leemos catalogo USNO */
    FILE *stream = fopen("cat/usno.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read usno.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
        /* omite la estrella si es una doble */
        readField(buffer, cell, 6, 1);
        if (cell[0] != ' ') continue;

        /* descarta si faltan los segundos de RA o Decl */
        readField(buffer, cell, 44, 1);
        if (cell[0] == ' ') continue;
        readField(buffer, cell, 65, 1);
        if (cell[0] == ' ') continue;

		/* lee ascension recta B1860.0 */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 40, &RAh, &RAm, &RAs);

		/* lee declinacion B1860.0 */
        int Decld, Declm, Decls;
        double Decl = readDeclField(buffer, 61, 2, 63, 65, 3, 10.0, &Decld, &Declm, &Decls);
		readField(buffer, cell, 60, 1);
		if (cell[0] == '-') Decl = -Decl;

        formatCatLine(catLine, RAh, RAm, RAs, cell[0], Decld, Declm, Decls);

		/* lee numeracion */
		readField(buffer, cell, 1, 5);
		int numRef = atoi(cell);

        snprintf(catName, 20, "U %d", numRef);

		/* lee magnitud, si es 0 indica ausencia de magnitud */
		readField(buffer, cell, 36, 3);
		float vmag = atof(cell)/10.0;
        if (vmag < __FLT_EPSILON__) {
            vmag = 0.0;
        }

		/* busca la PPM mas cercana y genera el cruzamiento */
        double x, y, z, minDistance;
        int ppmIndex;
        sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_USNO_PPM, catName, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;

        /* si no encuentra una PPM cercana, prueba con GSC */
        bool gscFound = tryGSC(ppmFound, RA, Decl, EPOCH_USNO, crossPPMStream, catName, vmag, &stats);

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_USNO, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        /* busca la CPD y la CD mas cercanas y genera los cruzamientos */
        bool cpdFound = false;
        if (Decl1875 <= -18.0)
            cpdFound = crossWithCPD(x, y, z, Decl1875, vmag, MAX_DIST_USNO_CPD, crossCPDStream, catName, &stats, NULL, NULL);

        bool cdFound = false;
        if (Decl1875 <= -22.0)
            cdFound = crossWithCD(x, y, z, Decl1875, vmag, catName, crossCDStream, catName, &stats, NULL, NULL);

        /* lee referencia numerica y referencia a catalogo (que queda en "cell") */
		readField(buffer, cell, 30, 5);
        int numRefCat = atoi(cell);
		readField(buffer, cell, 19, 4);

        if (!strncmp(cell, "LAC ", 4))
            checkCrossRef(catName, catLine, "L", x, y, z, numRefCat, &lacList, false, -1, &checkLac, &stats.errors);
        if (!strncmp(cell, "TAY ", 4))
            checkCrossRef(catName, catLine, "T", x, y, z, numRefCat, &tayList, false, -1, &checkTaylor, &stats.errors);
        if (!strncmp(cell, "LAL ", 4))
            checkCrossRef(catName, catLine, "Lal", x, y, z, numRefCat, &lalList, false, -1, &checkLal, &stats.errors);
        if (!strncmp(cell, "OARS", 4))
            checkCrossRef(catName, catLine, "OA", x, y, z, numRefCat, &oaList, false, -1, &checkOA, &stats.errors);
        if (!strncmp(cell, "BRI ", 4))
            checkCrossRef(catName, catLine, "B", x, y, z, numRefCat, &briList, true, -1, &checkBri, &stats.errors);

        if (!ppmFound && !cdFound && !cpdFound) {
            if (gscFound) {
                writeUnidentified(unidentifiedStream, catName, x, y, z);
            } else {
                warnAlone(&stats.errors, catName, catName, catLine,
                    catName, RAs, Decl1875, Decls,
                    PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
            }
        }

        /* la almacenamos para futuras identificaciones */
        storeStar(&countUsno, MAXUSNOSTAR, "USNO", usnoRef, usnoX, usnoY, usnoZ, usnoMag, numRef, x, y, z, vmag);
        writeCatalogFile(catalogStream, catName, x, y, z, vmag);
    }
    fclose(stream);
    fclose(unidentifiedStream);
    fclose(catalogStream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available USNO stars = %d\n", countUsno);
    printf("Stars from USNO with Lacaille = %d, Brisbane = %d, Lalande = %d, Taylor = %d and OA = %d\n",
        checkLac, checkBri, checkLal, checkTaylor, checkOA);
    printf("Stars from USNO identified with PPM = %d, GSC-PPM = %d, CD = %d and CPD = %d\n", stats.countDist, stats.countGSC, stats.countCD, stats.countCPD);
    printRSMEDist(&stats);
    printRSMEMag(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * readGCScanned - revisa referencias cruzadas de las páginas escaneadas de GC
 * (ya se deben haber leidos los catalogs OA, Lalande, Brisbane, Stone, Lacaille,
 *  GC, Taylor y Yarnall)
 * Lee todos los archivos scans/rnao14_page*.csv y, para cada estrella (una fila),
 * calcula la distancia entre la estrella GC (1a columna; la 2a columna -magnitud-
 * se ignora) y la estrella de referencia (3a columna) en los catalogs
 * Oeltzen-Argelander (OA.), Lalande (Ll.), Brisbane (B.), Stone (St.),
 * Lacaille (L.), Taylor (T.), Yarnall/USNO (Y.) y Circumpolar List (CL.).
 * Otras designaciones se ignoran.
 */
void readGCScanned() {
    char buffer[1024], ref[64], catgName[20];

    printf("\n***************************************\n");
    printf("Perform cross-checking of scanned GC pages...\n");

    int errors = 0;
    int checkLac = 0, checkLal = 0, checkOA = 0, checkTaylor = 0;
    int checkUSNO = 0, checkBri = 0, checkSt = 0, checkCL = 0;
    int countStars = 0, countRefs = 0;
    int countZCsaved = 0;

    /* coordenadas (1875.0) de las estrellas GC ya leídas */
    struct GCstar_struct *GCstar = getGCStruct();

    /* leemos catalogo PPM a la época 1875.0 (la del Zone Catalog de Gould) para
     * poder identificar las estrellas ZC referidas en las páginas escaneadas */
    preparePPM(1875.0, false);

    /* recorremos las páginas escaneadas */
    for (int page = 1; page <= 616; page++) {
        char filename[64];
        snprintf(filename, 64, "scans/rnao14_page%d.csv", page);
        FILE *stream = fopen(filename, "rt");
        if (stream == NULL) continue;

        while (fgets(buffer, 1023, stream) != NULL) {
            int gcRef;
            if (!parseGCScanLine(buffer, &gcRef, ref)) continue;

            snprintf(catgName, 20, "GC %d", gcRef);

            /* coordenadas (1875.0) de la estrella GC */
            double x, y, z;
            int gcIndex = -1;
            if (!getGCStarData(gcRef, &gcIndex, &x, &y, &z)) continue;
            countStars++;

            /* despachamos según el catálogo referido en la 3a columna */
            if (!strncmp(ref, "Ll.", 3)) {
                countRefs++;
                checkCrossRef(catgName, NULL, "Lal", x, y, z, atoi(&ref[3]), &lalList, false, gcIndex, &checkLal, &errors);
            } else if (!strncmp(ref, "OA.", 3)) {
                countRefs++;
                checkCrossRef(catgName, NULL, "OA", x, y, z, atoi(&ref[3]), &oaList, false, gcIndex, &checkOA, &errors);
            } else if (!strncmp(ref, "St.", 3)) {
                countRefs++;
                checkCrossRef(catgName, NULL, "St", x, y, z, atoi(&ref[3]), &stList, false, gcIndex, &checkSt, &errors);
            } else if (!strncmp(ref, "L.", 2)) {
                countRefs++;
                checkCrossRef(catgName, NULL, "L", x, y, z, atoi(&ref[2]), &lacList, false, gcIndex, &checkLac, &errors);
            } else if (!strncmp(ref, "B.", 2)) {
                countRefs++;
                checkCrossRef(catgName, NULL, "B", x, y, z, atoi(&ref[2]), &briList, true, gcIndex, &checkBri, &errors);
            } else if (!strncmp(ref, "T.", 2)) {
                countRefs++;
                checkCrossRef(catgName, NULL, "T", x, y, z, atoi(&ref[2]), &tayList, false, gcIndex, &checkTaylor, &errors);
            } else if (!strncmp(ref, "Y.", 2)) {
                countRefs++;
                checkYarnallRef(catgName, NULL, atoi(&ref[2]), x, y, z, true, gcIndex, &usnoList, &checkUSNO, &errors);
            } else if (!strncmp(ref, "CL.", 3)) {
                countRefs++;
                checkCrossRef(catgName, NULL, "CL", x, y, z, atoi(&ref[3]), &clList, false, gcIndex, &checkCL, &errors);
            } else if (!strncmp(ref, "ZC.", 3)) {
                /* guardamos la estrella en el Zone Catalog de Gould. */
                int numRefCat = atoi(&ref[3]);
                int RAh = GCstar[gcIndex].RAh;
                int RAs = GCstar[gcIndex].RAs;
                int Decls = GCstar[gcIndex].Decls;
                double Decl1875 = GCstar[gcIndex].Decl1875;

                /* busca la PPM mas cercana (a 1875.0) */
                bool ppmFound = false;
                int ppmIndex = -1;
                double minDistance = HUGE_NUMBER;
                findPPMByCoordinates(x, y, z, Decl1875, &ppmIndex, &minDistance);
                double nearestPPMDistance = minDistance;
                if (minDistance < MAX_DIST_PPM) ppmFound = true;

                /* si no hay PPM cercana, prueba con GSC */
                bool gscFound = false;
                const char *gscId = NULL;
                double nearestGSCDistance = HUGE_NUMBER;
                if (!ppmFound) {
                    gscFound = findGSCStar(GCstar[gcIndex].RA1875, Decl1875, 1875.0, MAX_DIST_GSC);
                    if (gscFound) {
                        gscId = getGSCId();
                        nearestGSCDistance = getDist();
                    }
                }

                /* busca la CPD mas cercana (solo dentro de la faja CPD, Decl <= -18) */
                bool cpdFound = false;
                int cpdIndex = -1;
                double nearestCPDDistance = HUGE_NUMBER;
                if (GCstar[gcIndex].Decld >= 18) {
                    minDistance = HUGE_NUMBER;
                    findCPDByCoordinates(x, y, z, Decl1875, &cpdIndex, &minDistance);
                    nearestCPDDistance = minDistance;
                    if (minDistance < MAX_DIST_CPD) cpdFound = true;
                }

                /* busca la CD mas cercana (solo dentro de la faja CD, Decl <= -22) */
                bool cdFound = false;
                int cdIndex = -1;
                double nearestCDDistance = HUGE_NUMBER;
                if (GCstar[gcIndex].Decld >= 22) {
                    minDistance = HUGE_NUMBER;
                    findDMByCoordinates(x, y, z, Decl1875, &cdIndex, &minDistance);
                    nearestCDDistance = minDistance;
                    if (minDistance < MAX_DIST_CD) cdFound = true;
                }

                saveZC(RAh, RAs, Decls, numRefCat, x, y, z,
                    ppmFound, ppmIndex, nearestPPMDistance,
                    gscFound, gscId, nearestGSCDistance,
                    cdFound, cdIndex, nearestCDDistance,
                    cpdFound, cpdIndex, nearestCPDDistance);
                countZCsaved++;
            }
            /* otras designaciones (WB, UA, F, P, G, M, Melb.I, nombres, ...) se ignoran */
        }
        fclose(stream);
    }

    printf("Scanned GC rows with a reference = %d; references to checked catalogs (OA/Ll/B/St/L/T/Y/CL) = %d\n",
        countStars, countRefs);
    printf("Cross-checks OK: Lacaille = %d, Brisbane = %d, Lalande = %d, Taylor = %d, OA = %d, Stone = %d, USNO = %d, CL = %d\n",
        checkLac, checkBri, checkLal, checkTaylor, checkOA, checkSt, checkUSNO, checkCL);
    printf("ZC (Gould's Zone Catalog) stars saved = %d\n", countZCsaved);
    printf("Errors logged = %d\n", errors);
}

/*
 * checkUACatalogRef - compara la referencia SAO/HD reportada por la UA con la
 * de la estrella PPM asociada, reporta discrepancias y escribe el cruzamiento
 */
static void checkUACatalogRef(char *buffer, int col, bool isHD, char *catName, double vmag,
        double x, double y, double z, char *ppmName, int ppmIndex, double nearestPPMDistance,
        double minDistance, FILE *stream, int *errors) {
    char cell[256], refName[20];
    struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();
    const char *label = isHD ? "HD" : "SAO";

    int catRef = isHD ? PPMstar[ppmIndex].hdRef : PPMstar[ppmIndex].saoRef;
    readFieldSanitized(buffer, cell, col, 6);
    int refFromUA = atoi(cell);
    if (catRef <= 0 || refFromUA <= 0) return;

    bool mismatch = false;
    if (refFromUA < catRef - 1 || refFromUA > catRef + 1) {
        printf("%d) Warning: %s references mismatch for %s: %d != %d.\n",
            ++(*errors),
            label,
            catName,
            catRef,
            refFromUA);
        mismatch = true;
    } else if (refFromUA != catRef) {
        printf("**) Note: %s references mismatch for %s: %d != %d.\n",
            label,
            catName,
            catRef,
            refFromUA);
        mismatch = true;
    }

    // Find cause of mismatch
    bool foundInPPM = false;
    if (mismatch) {
        for (int i = 0; i < PPMstars; i++) {
            int otherRef = isHD ? PPMstar[i].hdRef : PPMstar[i].saoRef;
            if (otherRef == refFromUA) {
                double dist = 3600.0 * calcAngularDistance(x, y, z, PPMstar[i].x, PPMstar[i].y, PPMstar[i].z);
                printf("     Match suggested by UA (vmag=%.1f): PPM %d (vmag=%.1f) at distance %.1f arcsec, instead of %s (vmag=%.1f) at %.1f arcsec\n",
                    vmag,
                    PPMstar[i].ppmRef,
                    PPMstar[i].vmag,
                    dist,
                    ppmName,
                    PPMstar[ppmIndex].vmag,
                    nearestPPMDistance);
                foundInPPM = true;
                break;
            }
        }
        if (!foundInPPM) {
            printf("     No match found in PPM for %s %d referenced by UA.\n", label, refFromUA);
        }
    }

    if (!mismatch || !foundInPPM) {
        // Only save if both matches or if the mismatch comes from the fact that the star is not in PPM
        snprintf(refName, 20, "%s %d", label, refFromUA);
        writeCrossEntry(stream, catName, refName, vmag, minDistance);
    }
}

/*
 * readUA - lee y cruza Uranometria Argentina
 * (ya se deben haber leidos los catalogs CD, CPD, Weiss, Brisbane, Stone y Yarnall)
 * tambien revisa referencias cruzadas; no tiene en cuenta magnitudes
 */
void readUA() {
    char buffer[1024], cell[256], catName[20], ppmName[20];

    printf("\n***************************************\n");
    printf("Perform comparison between UA and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_ua_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_ua_cpd.csv");
    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("ua", &crossPPMStream, &crossSAOStream, &crossHDStream);
    FILE *catalogStream = openCatalogFile("likelihood/cat1875/ua.csv");

    CrossStats stats;
    int checkLac = 0;
    int checkLal = 0;
    int checkOA = 0;
    int checkTaylor = 0;
    int checkUSNO = 0;
    int checkBri = 0;
    int countUA = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(EPOCH_UA, false);

    /* leemos catalogo UA */
    FILE *stream = fopen("cat/ua.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read ua.txt");
		exit(1);
    }
    char serpens = 'a';
    UAstar_struct ua;
    while (fgets(buffer, 1023, stream) != NULL) {
        if (!parseUALine(buffer, &serpens, &ua)) continue;

        /* nombre para los cruzamientos: Bayer, Flamsteed o Gould */
        if (ua.existsRef) {
            /* si está disponible Bayer, usa esa designacion */
            readField(buffer, cell, 15, 8);
            if (cell[0] != ' ' && !(cell[0] == 'N' && cell[1] == 'G' && cell[2] == 'C')) {
                char bayerRef[9];
                copyWithoutSpaces(bayerRef, cell);
                snprintf(catName, 20, "%s%s", bayerRef, ua.cstRef);
            } else {
                /* si está disponible Flamsteed, la usa */
                readField(buffer, cell, 11, 3);
                if (cell[2] != ' ') {
                    int fRef = atoi(cell);
                    snprintf(catName, 20, "%d %s", fRef, ua.cstRef);
                } else {
                    /* caso contrario, usa la denominación de Gould */
                    snprintf(catName, 20, "%dG %s", ua.gouldRef, ua.cstRef);
                }
            }
        } else {
            snprintf(catName, 11, "Annonymous");
        }

        // Lee magnitud (puede ser un entero, un flotante, una fracción o una palabra)
        double vmag = 0.0;
        readField(buffer, cell, 120, 5);
        if (strstr(cell, "var.") != NULL || strstr(cell, "cum.") != NULL || strstr(cell, "neb.") != NULL) {
            vmag = 0.0;
        } else {
            char *fracPtr = strstr(cell, " 1/");
            if (fracPtr == NULL) {
                fracPtr = strstr(cell, " 3/");
            }
            if (fracPtr != NULL) {
                // Tiene fracción: interpreta la parte entera
                char intPart[10];
                int len = fracPtr - cell;
                strncpy(intPart, cell, len);
                intPart[len] = '\0';
                vmag = atof(intPart);

                // Agrega la fracción
                if (strstr(fracPtr, "1/2") != NULL) {
                    vmag += 0.5;
                } else if (strstr(fracPtr, "1/4") != NULL) {
                    vmag += 0.2;
                } else if (strstr(fracPtr, "3/4") != NULL) {
                    vmag += 0.8;
                }
            } else {
                // No tiene fracción: interpreta como número normal
                vmag = atof(cell);
            }
        }

		/* busca la PPM mas cercana y genera el cruzamiento; tambien
           revisa las referencias SAO y HD reportadas por la UA */
        double x, y, z, minDistance;
        int ppmIndex;
        sph2rec(ua.RA, ua.Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, ua.Decl, vmag, MAX_DIST_UA_PPM, NULL, NULL,
            NULL, NULL, NULL, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;
        if (ppmFound && ua.existsRef) {
		    snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
		    writeCrossEntry(crossPPMStream, catName, ppmName, vmag, minDistance);
            checkUACatalogRef(buffer, 74, false, catName, vmag, x, y, z, ppmName, ppmIndex,
                nearestPPMDistance, minDistance, crossSAOStream, &stats.errors);
            checkUACatalogRef(buffer, 66, true, catName, vmag, x, y, z, ppmName, ppmIndex,
                nearestPPMDistance, minDistance, crossHDStream, &stats.errors);
		}

        /* busca la CPD y la CD mas cercanas y genera los cruzamientos */
        bool cpdFound = false;
        if (ua.Decl <= -18.0)
            cpdFound = crossWithCPD(x, y, z, ua.Decl, vmag, MAX_DIST_CPD, crossCPDStream,
                ua.existsRef ? catName : NULL, &stats, NULL, NULL);

        bool cdFound = false;
        if (ua.Decl <= -22.0)
            cdFound = crossWithCD(x, y, z, ua.Decl, vmag, NULL, crossCDStream,
                ua.existsRef ? catName : NULL, &stats, NULL, NULL);

        /* lee referencias a catalogos. */
        char subcell[3][18];
        int count = splitUARefs(buffer, subcell);

        for (int refs = 0; refs < count; refs++) {
            if (!strncmp(subcell[refs], "L.", 2) && subcell[refs][2] != '(')
                checkCrossRef(ua.catgName, ua.catLine, "L", x, y, z, atoi(&subcell[refs][2]), &lacList, false, -1, &checkLac, &stats.errors);
            if (!strncmp(subcell[refs], "B.", 2))
                checkCrossRef(ua.catgName, ua.catLine, "B", x, y, z, atoi(&subcell[refs][2]), &briList, true, -1, &checkBri, &stats.errors);
            if (!strncmp(subcell[refs], "T.", 2))
                checkCrossRef(ua.catgName, ua.catLine, "T", x, y, z, atoi(&subcell[refs][2]), &tayList, false, -1, &checkTaylor, &stats.errors);
            if (!strncmp(subcell[refs], "Ll.", 3) && subcell[refs][3] != '(')
                checkCrossRef(ua.catgName, ua.catLine, "Lal", x, y, z, atoi(&subcell[refs][3]), &lalList, false, -1, &checkLal, &stats.errors);
            if (!strncmp(subcell[refs], "OA.", 3))
                checkCrossRef(ua.catgName, ua.catLine, "OA", x, y, z, atoi(&subcell[refs][3]), &oaList, false, -1, &checkOA, &stats.errors);
            if (!strncmp(subcell[refs], "Y.", 2))
                checkYarnallRef(ua.catgName, ua.catLine, atoi(&subcell[refs][2]), x, y, z, false, -1, &usnoList, &checkUSNO, &stats.errors);
        }

        if (!ppmFound && !cdFound && !cpdFound) {
            printf("%d) Warning: %s is ALONE (no PPM / CD / CPD star near it).\n",
                ++stats.errors,
                ua.catgName);
            printf("     Register %s: %s\n", ua.catgName, ua.catLine);
            readField(buffer, cell, 132, 3);
            bool cumulus = !strncmp(cell, "cum", 3);
            bool nebula = !strncmp(cell, "neb", 3);
            logCauses(catName, true,
                cumulus, nebula, ua.RAs, ua.Decl, 1,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        if (ua.existsRef) {
            writeCatalogFile(catalogStream, catName, x, y, z, vmag);
            countUA++;
        }
    }
    fclose(stream);
    fclose(catalogStream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available UA stars = %d\n", countUA);
    printf("Stars from UA with Lacaille = %d, Brisbane = %d, Lalande = %d, Taylor = %d, OA = %d and USNO = %d\n",
        checkLac, checkBri, checkLal, checkTaylor, checkOA, checkUSNO);
    printf("Stars from UA identified with PPM = %d, CD = %d and CPD = %d\n", stats.countDist, stats.countCD, stats.countCPD);
    printRSMEDist(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * readCL - lee catalogo de la lista de estrellas circumpolares (Circumpolar List)
 * Cada fila de cat/CL.csv es una estrella: numero, designacion, magnitud,
 * ascension recta (h, m, s) y declinacion (grados, minutos) para 1872.0.
 */
void readCL() {
    char buffer[1024], name[32], magStr[16];

    printf("\n***************************************\n");
    printf("Read RNAO2 Circumpolar List (CL)...\n");

    FILE *stream = fopen("cat/CL.csv", "rt");
    if (stream == NULL) {
        perror("Cannot read CL.csv");
        exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
        int clNum, RAh, RAm, Decld;
        double RAs, Declm;

        /* lee la fila; la cabecera y las lineas vacias no matchean y se saltan */
        int n = sscanf(buffer, "%d,%31[^,],%15[^,],%d,%d,%lf,%d,%lf",
            &clNum, name, magStr, &RAh, &RAm, &RAs, &Decld, &Declm);
        if (n != 8) continue;

        /* ascension recta 1872.0 (h, m, s) a grados */
        double RA = ((double) RAh) + ((double) RAm)/60.0 + RAs/3600.0;
        RA *= 15.0; /* conversion horas a grados */

        /* declinacion 1872.0 (siempre negativa en esta lista) */
        double Decl = (double) Decld;
        Decl = Decl < 0.0 ? Decl - Declm/60.0 : Decl + Declm/60.0;

        /* convierte coordenadas de 1872.0 a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_CL, 1875.0, &RA1875, &Decl1875);
        double x, y, z;
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        /* la almacenamos para futuras identificaciones */
        storeStar(&countCL, MAXCLSTAR, "CL", clRef, clX, clY, clZ, NULL, clNum, x, y, z, 0.0);
    }
    fclose(stream);

    printf("Available CL stars = %d\n", countCL);
}

/*
 * readThome - lee y cruza catalogo de los Resultados 15
 * (ya se deben haber leidos los catalogs CD, CPD, GC, Yarnall, Brisbane y Stone)
 * tambien revisa referencias cruzadas a GC, OA, Yarnall, Lacaille, Brisbane y Stone
 * tambien genera identificaciones cruzadas para aquellas estrellas ZC
 */
void readThome(double epoch, const char *filename, int correction) {
    char buffer[1024], cell[256], srcName[24];
    char catLine[64];

    printf("\n***************************************\n");
    printf("Perform comparison between Thome (Resultados ONA 15) and PPM/CD/CPD...\n");

    CrossStats stats;
    int checkLac = 0;
    int checkLal = 0;
    int checkTaylor = 0;
    int checkOA = 0;
    int checkStone = 0;
    int checkGC = 0;
    int checkUSNO = 0;
    int checkBri = 0;
    int checkCL = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(epoch, false);

    /* leemos catalogo Thome */
    FILE *stream = fopen(filename, "rt");
    if (stream == NULL) {
        perror("Cannot read Thome catalog");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		/* lee ascension recta */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 20 + correction, &RAh, &RAm, &RAs);

		/* lee declinacion (incorpora signo negativo; en nuestro caso, siempre) */
        int Decld, Declm, Decls;
        double Decl = -readDeclField(buffer, 45 + correction, 2, 47 + correction, 49 + correction, 3, 10.0,
            &Decld, &Declm, &Decls);

        formatCatLine(catLine, RAh, RAm, RAs, '-', Decld, Declm, Decls);

		/* lee numeración catálogo y magnitud */
		readField(buffer, cell, 1, 4 + correction);
		int numRef = atoi(cell);
		readField(buffer, cell, 13 + correction, 3);
        if (cell[2] == ' ') cell[2] = '0';
		float vmag = atof(cell)/10.0;

        snprintf(srcName, 24, "T %.0f %d", epoch, numRef);

		/* busca la PPM mas cercana */
        double x, y, z, minDistance;
        int ppmIndex;
        sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_PPM, srcName, NULL,
            NULL, NULL, NULL, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;

        /* si no encuentra una PPM cercana, prueba con GSC */
        bool gscFound = false;
        const char *gscId = NULL;
        double nearestGSCDistance = HUGE_NUMBER;
        if (!ppmFound) {
            gscFound = findGSCStar(RA, Decl, epoch, MAX_DIST_GSC);
            if (gscFound) {
                gscId = getGSCId();
                nearestGSCDistance = getDist();
                stats.countGSC++;
            }
        }

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(epoch, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        /* busca la CPD y la CD mas cercanas */
        int cpdIndex = -1;
        double nearestCPDDistance;
        bool cpdFound = crossWithCPD(x, y, z, Decl1875, vmag, MAX_DIST_TH_CPD, NULL, NULL,
            &stats, &cpdIndex, &nearestCPDDistance);

        int cdIndex = -1;
        double nearestCDDistance;
        bool cdFound = crossWithCD(x, y, z, Decl1875, vmag, srcName, NULL, NULL,
            &stats, &cdIndex, &nearestCDDistance);

        if (!ppmFound && !cdFound && !cpdFound && !gscFound) {
            warnAlone(&stats.errors, srcName, srcName, catLine,
                NULL, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* lee referencia numerica y referencia a catalogo (que queda en "cell") */
		readField(buffer, cell, 68 + correction, 5);
        int numRefCat = atoi(cell);
		readField(buffer, cell, 66 + correction, 2);

        if (!strncmp(cell, "L ", 2))
            checkCrossRef(srcName, catLine, "L", x, y, z, numRefCat, &lacList, false, -1, &checkLac, &stats.errors);
        if (!strncmp(cell, "T ", 2))
            checkCrossRef(srcName, catLine, "T", x, y, z, numRefCat, &tayList, false, -1, &checkTaylor, &stats.errors);
        if (!strncmp(cell, "LL", 2))
            checkCrossRef(srcName, catLine, "Lal", x, y, z, numRefCat, &lalList, false, -1, &checkLal, &stats.errors);
        if (!strncmp(cell, "ST", 2))
            checkCrossRef(srcName, catLine, "St", x, y, z, numRefCat, &stList, false, -1, &checkStone, &stats.errors);
        if (!strncmp(cell, "CL", 2))
            checkCrossRef(srcName, catLine, "CL", x, y, z, numRefCat, &clList, false, -1, &checkCL, &stats.errors);
        if (!strncmp(cell, "OA", 2))
            checkCrossRef(srcName, catLine, "OA", x, y, z, numRefCat, &oaList, false, -1, &checkOA, &stats.errors);
        if (!strncmp(cell, "GC", 2))
            checkCrossRefGCcat(srcName, catLine, numRefCat, x, y, z, &checkGC, &stats.errors);
        if (!strncmp(cell, "B ", 2))
            checkCrossRef(srcName, catLine, "B", x, y, z, numRefCat, &briList, false, -1, &checkBri, &stats.errors);
        if (!strncmp(cell, "Y ", 2))
            checkYarnallRef(srcName, catLine, numRefCat, x, y, z, false, -1, &usnoList, &checkUSNO, &stats.errors);

        if (!strncmp(cell, "ZC", 2)) {
            // save Gould's Zone Catalog
            int zcRAh = (int) (floor(RA1875 / 15.0) + __FLT_EPSILON__);
            saveZC(zcRAh, RAs, Decls, numRefCat, x, y, z,
                ppmFound, ppmIndex, nearestPPMDistance,
                gscFound, gscId, nearestGSCDistance,
                cdFound, cdIndex, nearestCDDistance,
                cpdFound, cpdIndex, nearestCPDDistance);
        }
    }
    fclose(stream);

    printf("Stars from Thome %.0f identified with PPM = %d, CD = %d and CPD = %d\n", epoch, stats.countDist, stats.countCD, stats.countCPD);
    printf("Stars from Thome %.0f with Lacaille = %d, Lalande = %d, OA = %d, Taylor = %d, Stone = %d, GC = %d, USNO = %d, Brisbane = %d and CL = %d\n",
        epoch, checkLac, checkLal, checkOA, checkTaylor, checkStone, checkGC, checkUSNO, checkBri, checkCL);
    printRSMEDist(&stats);
    printRSMEMag(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * readGilliss - lee y cruza catalogo de Gilliss
 * (ya se deben haber leidos los catalogs CD, CPD, GC y Stone)
 * tambien revisa referencias cruzadas a GC, Stone, Brisbane y Lacaille
 * tambien genera identificaciones cruzadas para aquellas estrellas ZC
 */
void readGilliss() {
    char buffer[1024], cell[256], catName[20];
    char catLine[64];

    printf("\n***************************************\n");
    printf("Perform comparison between Gilliss and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_gilliss_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_gilliss_cpd.csv");
    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("gilliss", &crossPPMStream, &crossSAOStream, &crossHDStream);
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/gilliss_unidentified.csv");
    FILE *catalogStream = openCatalogFile("likelihood/cat1875/gilliss.csv");

    CrossStats stats;
    int checkLac = 0;
    int checkGC = 0;
    int checkStone = 0;
    int checkBri = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(EPOCH_GILLISS, false);

    /* leemos catalogo Gilliss */
    FILE *stream = fopen("cat/gilliss.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read gilliss.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
        /* descarta dobles */
        readField(buffer, cell, 6, 1);
        if (cell[0] != ' ') {
            continue;
        }

		/* lee ascension recta B1850.0 */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 33, &RAh, &RAm, &RAs);

		/* lee declinacion B1850.0 (incorpora signo negativo; en nuestro caso, siempre) */
        int Decld, Declm, Decls;
        double Decl = -readDeclField(buffer, 49, 2, 51, 53, 3, 10.0, &Decld, &Declm, &Decls);

        formatCatLine(catLine, RAh, RAm, RAs, '-', Decld, Declm, Decls);

		/* lee numeración catálogo Gilliss */
		readField(buffer, cell, 1, 5);
		int giRef = atoi(cell);
        snprintf(catName, 20, "G %d", giRef);

        /* si esta disponible, tambien lee precesiones y chequea (usamos constantes de Struve) */
        readFieldSanitized(buffer, cell, 41, 7);
        checkPrecessionRA(atof(cell) / 1000.0, catName, catLine, RA, Decl, 3.07177, 1.3371, 0.0299, &stats.errors);
        readFieldSanitized(buffer, cell, 56, 5);
        checkPrecessionDecl(atof(cell) / 100.0, catName, catLine, RA, 20.0564, 0.099, &stats.errors);

        /* lee magnitud */
		readField(buffer, cell, 29, 3);
		float vmag = atof(cell)/10.0;

		/* busca la PPM mas cercana y genera el cruzamiento */
        double x, y, z, minDistance;
        int ppmIndex;
        sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_PPM, catName, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);
        double nearestPPMDistance = minDistance;

        /* si no encuentra una PPM cercana, prueba con GSC */
        const char *gscId = NULL;
        double nearestGSCDistance = HUGE_NUMBER;
        bool gscFound = tryGSC(ppmFound, RA, Decl, EPOCH_GILLISS, crossPPMStream, catName, vmag, &stats);
        if (gscFound) {
            gscId = getGSCId();
            nearestGSCDistance = getDist();
        }

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_GILLISS, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        /* la almacenamos para futuras identificaciones */
        storeStar(&countGil, MAXGILLISSSTAR, "Gilliss", gilRef, gilX, gilY, gilZ, gilMag, giRef, x, y, z, vmag);
        writeCatalogFile(catalogStream, catName, x, y, z, vmag);

        /* busca la CPD y la CD mas cercanas y genera los cruzamientos */
        int cpdIndex = -1;
        double nearestCPDDistance;
        bool cpdFound = crossWithCPD(x, y, z, Decl1875, vmag, MAX_DIST_CPD, crossCPDStream, catName,
            &stats, &cpdIndex, &nearestCPDDistance);

        int cdIndex = -1;
        double nearestCDDistance;
        bool cdFound = crossWithCD(x, y, z, Decl1875, vmag, catName, crossCDStream, catName,
            &stats, &cdIndex, &nearestCDDistance);

        if (!ppmFound && !cdFound && !cpdFound) {
            if (gscFound) {
                writeUnidentified(unidentifiedStream, catName, x, y, z);
            } else {
                warnAlone(&stats.errors, catName, catName, catLine,
                    catName, RAs, Decl1875, Decls,
                    PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
            }
        }

        /* lee referencia numerica y referencia a catalogo (que queda en "cell") */
		readField(buffer, cell, 20, 5);
        int numRefCat = atoi(cell);
		readField(buffer, cell, 15, 3);

        if (!strncmp(cell, "LAC", 3))
            checkCrossRef(catName, catLine, "L", x, y, z, numRefCat, &lacList, false, -1, &checkLac, &stats.errors);
        if (!strncmp(cell, "STO", 3))
            checkCrossRef(catName, catLine, "St", x, y, z, numRefCat, &stList, false, -1, &checkStone, &stats.errors);
        if (!strncmp(cell, "GOU", 3))
            checkCrossRefGCcat(catName, catLine, numRefCat, x, y, z, &checkGC, &stats.errors);
        if (!strncmp(cell, "BRI", 3))
            checkCrossRef(catName, catLine, "B", x, y, z, numRefCat, &briList, false, -1, &checkBri, &stats.errors);

        if (!strncmp(cell, "GZC", 3)) {
            // save Gould's Zone Catalog
            readField(buffer, cell, 18, 2);
            int zcRAh = atoi(cell);
            int originalRAh = (int) floor(RA1875/15.0 + __FLT_EPSILON__);
            if (zcRAh != originalRAh) {
                printf("%d) Warning: G %d has a different RA zone (%d) than identified ZC (%d).\n",
                    ++stats.errors,
                    giRef,
                    zcRAh,
                    originalRAh);
            }
            saveZC(zcRAh, RAs, Decls, numRefCat, x, y, z,
                ppmFound, ppmIndex, nearestPPMDistance,
                gscFound, gscId, nearestGSCDistance,
                cdFound, cdIndex, nearestCDDistance,
                cpdFound, cpdIndex, nearestCPDDistance);
        }
    }
    fclose(stream);
    fclose(unidentifiedStream);
    fclose(catalogStream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Stars from Gilliss identified with PPM = %d, GSC-PPM = %d, CD = %d and CPD = %d\n", stats.countDist, stats.countGSC, stats.countCD, stats.countCPD);
    printf("Stars from Gilliss with Lacaille = %d, Stone = %d, GC = %d and Brisbane = %d\n",
        checkLac, checkStone, checkGC, checkBri);
    printRSMEDist(&stats);
    printRSMEMag(&stats);
    printf("Errors logged = %d\n", stats.errors);
}

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("CROSS_SOUTH - Compare several catalogs.\n");
    printf("Made in 2025 by Daniel Severin.\n");

    /* leemos catalogo CD */
    readDM(CURATED ? "cat/cd_curated.txt" : "cat/cd.txt");

    /* leemos catalogo CPD */
    readCPD(false, false);

    /* leemos y cruzamos segundo catálogo argentino */
    readGC2();

    /* leemos y cruzamos Weiss / OA */
    readWeiss();
    makeDoubles(countOA, oaRef, oaX, oaY, oaZ, oaMag, "OA", "results/doubles/oa.csv");

    /* leemos y cruzamos Lalande (solo contra PPM) */
    readLalande();

    /* revisamos las referencias a Lalande de los Standards of Magnitude */
    readSOM();

    /* leemos y cruzamos Brisbane */
    readBrisbane();

    /* leemos y cruzamos Stone / Lacaille */
    readStone();

    /* leemos catalogo GC */
	readGC();

    /* leemos y cruzamos Taylor */
    readTaylor();
    makeDoubles(countTaylor, tayRef, tayX, tayY, tayZ, tayMag, "T", "results/doubles/taylor.csv");

    /* leemos, cruzamos y revisamos identificaciones de Yarnall */
    readUSNO();
    makeDoubles(countUsno, usnoRef, usnoX, usnoY, usnoZ, usnoMag, "U", "results/doubles/usno.csv");

    /* leemos, cruzamos y revisamos Uranometria Argentina */
    readUA();

    crossPPMZCStream = openCrossFile("results/cross/cross_zc_ppm.csv");
    crossCDZCStream = openCrossFile("results/cross/cross_zc_cd.csv");
    crossCPDZCStream = openCrossFile("results/cross/cross_zc_cpd.csv");
    crossSAOZCStream = openCrossFile("results/cross/cross_zc_sao.csv");
    crossHDZCStream = openCrossFile("results/cross/cross_zc_hd.csv");
    unidentifiedZCStream = openUnidentifiedFile("results/cross/zc_unidentified.csv");

    /* leemos la lista de estrellas circumpolares (Circumpolar List) */
    readCL();

    /* leemos, cruzamos y revisamos identificaciones de Thome */
    /* tambien generamos identificaciones de Gould's Zone Catalog */
    readThome(1881.0, "cat/thome1881.txt", 0);
    readThome(1882.0, "cat/thome1882.txt", 0);
    readThome(1883.0, "cat/thome1883.txt", -1);
    readThome(1884.0, "cat/thome1884.txt", -1);

    /* leemos, cruzamos y revisamos identificaciones de Gilliss */
    /* tambien generamos identificaciones de Gould's Zone Catalog */
    readGilliss();
    makeDoubles(countGil, gilRef, gilX, gilY, gilZ, gilMag, "G", "results/doubles/gilliss.csv");

    /* hacemos cross-checkings de páginas escaneadas de GC
     * tambien generamos identificaciones de Gould's Zone Catalog */
    readGCScanned();

    fclose(unidentifiedZCStream);
	fclose(crossPPMZCStream);
	fclose(crossCPDZCStream);
	fclose(crossCDZCStream);
    printf("\nNumber of ZC stars registered = %d\n", countZC);
    printf("Stars from ZC identified with PPM = %d, GSC-PPM = %d, CD = %d and CPD = %d\n", countPPMZC, countGSCZC, countCDZC, countCPDZC);
    return 0;
}
