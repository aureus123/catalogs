
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

#define MAXOASTAR 19000
#define MAXLALSTAR 47400
#define MAXLACSTAR 10000
#define MAXSTSTAR 12500
#define MAXTAYLORSTAR 11100
#define MAXUSNOSTAR 11000
#define MAXZCSTAR 15000
#define EPOCH_GC2 1900.0
#define EPOCH_OA 1850.0
#define EPOCH_LAL 1800.0
#define EPOCH_ST 1880.0
#define EPOCH_TAYLOR 1835.0
#define EPOCH_USNO 1860.0
#define EPOCH_UA 1875.0
#define EPOCH_GILLISS 1850.0
#define MAX_DIST_PPM_FAR 90.0
#define MAX_DIST_OA_PPM 45.0
#define MAX_DIST_LAL_PPM 45.0
#define MAX_DIST_TAY_PPM 40.0
#define MAX_DIST_OA_CPD 45.0
#define MAX_DIST_ST_CPD 30.0
#define MAX_DIST_USNO_CPD 30.0
#define MAX_DIST_TH_CPD 30.0
#define MAX_DIST_CROSS 60.0
#define MAX_DIST_ZC_ZC 15.0
#define CURATED true // true if curated CD catalog should be used

// Here, we save 1875.0 coordinates of OA stars in rectangular form
double oaX[MAXOASTAR], oaY[MAXOASTAR], oaZ[MAXOASTAR];
int oaRef[MAXOASTAR];
int countOA = 0;

// Here, we save 1875.0 coordinates of Lalande stars in rectangular form
double lalX[MAXLALSTAR], lalY[MAXLALSTAR], lalZ[MAXLALSTAR];
int lalRef[MAXLALSTAR];
int countLal = 0;

// Here, we save 1875.0 coordinates of Lacaille stars in rectangular form
double lacX[MAXLACSTAR], lacY[MAXLACSTAR], lacZ[MAXLACSTAR];
int lacRef[MAXLACSTAR];
int countLac = 0;

// Here, we save 1875.0 coordinates of Taylor stars in rectangular form
double tayX[MAXLALSTAR], tayY[MAXLALSTAR], tayZ[MAXLALSTAR];
int tayRef[MAXLALSTAR];
int countTaylor = 0;

// Here, we save 1875.0 coordinates of Stone stars in rectangular form
double stX[MAXSTSTAR], stY[MAXSTSTAR], stZ[MAXSTSTAR];
int stRef[MAXSTSTAR];
int countSt = 0;

// Here, we save 1875.0 coordinates of USNO stars in rectangular form
double usnoX[MAXLACSTAR], usnoY[MAXLACSTAR], usnoZ[MAXLACSTAR];
int usnoRef[MAXLACSTAR];
int countUsno = 0;

// Here, we save 1875.0 coordinates of Gould's ZC stars in rectangular form
double zcX[MAXZCSTAR], zcY[MAXZCSTAR], zcZ[MAXZCSTAR];
int zcHour[MAXZCSTAR], zcNum[MAXZCSTAR];
int countZC = 0;
int countPPMZC = 0, countCDZC = 0, countCPDZC = 0;
FILE *crossPPMZCStream;
FILE *crossCDZCStream;
FILE *crossCPDZCStream;
FILE *unidentifiedZCStream;

void saveZC(int RAh, int RAs, int Decls,
        int numRef, double x, double y, double z,
        bool ppmFound, int ppmIndex, double minPPMDistance,
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
        snprintf(catName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
        writeCrossEntry(crossPPMZCStream, zcName, catName, minPPMDistance);
    }
    if (cdFound) {
        countCDZC++;
        struct DMstar_struct *CDstar = getDMStruct();
        snprintf(catName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
        writeCrossEntry(crossCDZCStream, zcName, catName, minCDDistance);
    }
    if (cpdFound) {
        countCPDZC++;
        struct CPDstar_struct *CPDstar = getCPDStruct();
        snprintf(catName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
        writeCrossEntry(crossCPDZCStream, zcName, catName, minCPDDistance);
    }

    if (!ppmFound && !cdFound && !cpdFound) {
        /* note: some of the parameters are fake in order to avoid log.
         * here, we just want to save the unidentified star. */
        printf("**) %s is also ALONE.\n", zcName);
        logCauses(zcName, unidentifiedZCStream, x, y, z,
            false, false, 10.0, RAs, -50.0, Decls, -1, 0.0);
    }

    /* add the new star */
    if (countZC >= MAXZCSTAR) {
        printf("Error: too many ZC stars.\n");
        exit(1);
    }
    zcX[countZC] = x;
    zcY[countZC] = y;
    zcZ[countZC] = z;
    zcHour[countZC] = RAh;
    zcNum[countZC] = numRef;
    countZC++;
}

/*
 * readGC2 - lee y cruza segundo catálogo argentino
 * (ya se deben haber leidos los catalogs CD y CPD)
 */
void readGC2() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];

    /* usamos catalogs CD y CPD */
    struct DMstar_struct *CDstar = getDMStruct();
    struct CPDstar_struct *CPDstar = getCPDStruct();

    printf("\n***************************************\n");
    printf("Perform comparison between GC2 and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_gc2_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_gc2_cpd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross/cross_gc2_ppm.csv");
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/gc2_unidentified.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countGC = 0;
    int countCD = 0;
    int countCPD = 0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_GC2);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

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
		readFieldSanitized(buffer, cell, 21, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 23, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 25, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1900.0 */
		readFieldSanitized(buffer, cell, 45, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 47, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 49, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */

		/* lee magnitud */
		readField(buffer, cell, 14, 2);
		float vmag = atof(cell);
		readFieldSanitized(buffer, cell, 16, 1);
        vmag += atof(cell)/10.0;

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_PPM) {
			float ppmVmag = PPMstar[ppmIndex].vmag;
			if (vmag > __FLT_EPSILON__ && ppmVmag > __FLT_EPSILON__) {
				// Note: no fit is performed to convert scales of magnitudes
				float delta = fabs(vmag - ppmVmag);
				if (delta < MAX_MAGNITUDE) {
                    akkuDeltaError += delta * delta;
                    countDelta++;
                } else {
					printf("%d) Warning: G2 %d (vmag = %.1f) has a near PPM %d star (vmag = %.1f) with a different magnitude (delta = %.1f).\n",
                        ++errors,
						gcRef,
                        vmag,
                        PPMstar[ppmIndex].ppmRef,
						ppmVmag,
						delta);
				}
			}

            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;

			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_GC2, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        if (Decl1875 <= -18.0) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            findCPDByCoordinates(x, y, z, Decl1875, &cpdIndex, &minDistance);
			if (minDistance < MAX_DIST_CPD) {
                countCPD++;
                cpdFound = true;

			    snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
			    writeCrossEntry(crossCPDStream, catName, cdName, minDistance);
            }
        }

        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
		if (Decl1875 <= -22.0) {
            int cdIndex = -1;
            minDistance = HUGE_NUMBER;
            findDMByCoordinates(x, y, z, Decl1875, &cdIndex, &minDistance);
			if (minDistance < MAX_DIST_CD) {
			    float cdVmag = CDstar[cdIndex].vmag;
			    if (vmag > __FLT_EPSILON__ && cdVmag < 29.9) {
				    float delta = fabs(vmag - cdVmag);
				    if (delta >= MAX_MAGNITUDE) {
					    printf("%d) Warning: G2 %d (vmag = %.1f) has a near CD with dif. magnitude (delta = %.1f). Check dpl.\n",
                            ++errors,
                            gcRef,
                            vmag,
                            delta);
					    writeRegister(cdIndex, false);
				    }
			    }

                countCD++;
                cdFound = true;

			    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
			    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
            }
		}

        if (!ppmFound && !cdFound && !cpdFound) {
            printf("%d) Warning: G2 %d is ALONE (no PPM or CD or CPD star near it).\n",
                ++errors,
                gcRef);
            logCauses(catName, unidentifiedStream, x, y, z,
                false, false, vmag, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }
        countGC++;
    }
    fclose(stream);
    fclose(unidentifiedStream);
	fclose(crossPPMStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available G2 stars = %d\n", countGC);
    printf("Stars from G2 identified with PPM = %d, CD = %d and CPD = %d\n", countDist, countCD, countCPD);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)countDelta),
        countDelta);
    printf("Errors logged = %d\n", errors);
}

/*
 * readWeiss - lee y cruza catalogo de Weiss
 * (ya se deben haber leidos los catalogs CD y CPD)
 */
void readWeiss() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];
    char catLine[64];

    /* usamos catalogs CD y CPD */
    struct DMstar_struct *CDstar = getDMStruct();
    struct CPDstar_struct *CPDstar = getCPDStruct();

    printf("\n***************************************\n");
    printf("Perform comparison between Weiss and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_oa_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_oa_cpd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross/cross_oa_ppm.csv");
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/oa_unidentified.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countCD = 0;
    int countCPD = 0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_OA);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

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

		/* lee ascension recta B1850.0 */
		readFieldSanitized(buffer, cell, 13, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 15, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 17, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1850.0 */
		readFieldSanitized(buffer, cell, 22, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 24, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 26, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */
        
        snprintf(catLine, 64, "%02dh %02dm %02ds%02d -%02d°%02d'%02d\"%01d",
            RAh, RAm, RAs / 100, RAs % 100, Decld, Declm, Decls / 10, Decls % 10);

		/* lee magnitud */
		readField(buffer, cell, 9, 1);
		float vmag = atof(cell);
		readField(buffer, cell, 11, 1);
		vmag += atof(cell)/10.0;

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_OA_PPM) {
			float ppmVmag = PPMstar[ppmIndex].vmag;
			if (vmag > __FLT_EPSILON__ && ppmVmag > __FLT_EPSILON__) {
				// Note: no fit is performed to convert scales of magnitudes
				float delta = fabs(vmag - ppmVmag);
				if (delta < MAX_MAGNITUDE) {
                    akkuDeltaError += delta * delta;
                    countDelta++;
                } else {
					printf("%d) Warning: W %d (vmag = %.1f) has a near PPM %d star (vmag = %.1f) with a different magnitude (delta = %.1f).\n",
                        ++errors,
						weissRef,
                        vmag,
                        PPMstar[ppmIndex].ppmRef,
						ppmVmag,
						delta);	
				}
			}

            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;

			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_OA, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        if (Decl1875 <= -18.0) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            findCPDByCoordinates(x, y, z, Decl1875, &cpdIndex, &minDistance);
			if (minDistance < MAX_DIST_OA_CPD) {
                countCPD++;
                cpdFound = true;

			    snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
			    writeCrossEntry(crossCPDStream, catName, cdName, minDistance);
            }
        }
		
        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
		if (Decl1875 <= -22.0) {
            int cdIndex = -1;
            minDistance = HUGE_NUMBER;
            findDMByCoordinates(x, y, z, Decl1875, &cdIndex, &minDistance);
			if (minDistance < MAX_DIST_CD) {
			    float cdVmag = CDstar[cdIndex].vmag;
			    if (vmag > __FLT_EPSILON__ && cdVmag < 29.9) {
				    float delta = fabs(vmag - cdVmag);
				    if (delta >= MAX_MAGNITUDE) {
					    printf("%d) Warning: W %d (vmag = %.1f) has a near CD with dif. magnitude (delta = %.1f). Check dpl.\n",
                            ++errors,
                            weissRef,
                            vmag,
                            delta);
					    writeRegister(cdIndex, false);
				    }
			    }

                countCD++;
                cdFound = true;

			    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
			    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
            }
		}

        if (!ppmFound && !cdFound && !cpdFound) {
            printf("%d) Warning: W %d (OA %d) is ALONE (no PPM or CD or CPD star near it).\n",
                ++errors,
                weissRef,
                oeltzenRef);
            printf("     Register W %d: %s\n", weissRef, catLine);    
            logCauses(catName, unidentifiedStream, x, y, z,
                false, false, vmag, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* la almacenamos para futuras identificaciones */
        if (countOA >= MAXOASTAR) {
            printf("Error: too many OA stars.\n");
            exit(1);
        }
        oaX[countOA] = x;
        oaY[countOA] = y;
        oaZ[countOA] = z;
        oaRef[countOA] = oeltzenRef;
        countOA++;
    }
    fclose(stream);
    fclose(unidentifiedStream);
	fclose(crossPPMStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available OA stars = %d\n", countOA);
    printf("Stars from OA identified with PPM = %d, CD = %d and CPD = %d\n", countDist, countCD, countCPD);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)countDelta),
        countDelta);
    printf("Errors logged = %d\n", errors);
}

/*
 * readLalande - lee catalogo de Lalande
 */
void readLalande() {
    char buffer[1024], cell[256], catName[20], ppmName[20];

    printf("\n***************************************\n");
    printf("Perform comparison between Lalande and PPM...\n");

    FILE *crossPPMStream = openCrossFile("results/cross/cross_lalande_ppm.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_LAL);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo Lalande */
    FILE *stream = fopen("cat/lalande.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read lalande.txt");
		exit(1);
    }

    double RA = 0.0;
    double Decl = 0.0;
    while (fgets(buffer, 1023, stream) != NULL) {
        /* descarta la estrella si no hay RA o declinación */
        readField(buffer, cell, 19, 1);
        if (cell[0] == ' ') continue;
        readField(buffer, cell, 32, 1);
        if (cell[0] == ' ') continue;

		/* lee numeración */
		readField(buffer, cell, 1, 5);
		int catRef = atoi(cell);

		/* lee ascension recta B1800.0 (si existiese) */
        readField(buffer, cell, 19, 1);
		readFieldSanitized(buffer, cell, 19, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 21, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 23, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1800.0 (en realidad NPD) */
		readFieldSanitized(buffer, cell, 32, 3);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 35, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 37, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		Decl = 90.0 - Decl; /* convertimos NPD a declinación */

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_LAL_PPM) {
            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;

            if (catRef > 0) {
			    snprintf(catName, 20, "Lal %d", catRef);
			    snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			    writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
            }
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_LAL, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        if (!ppmFound && nearestPPMDistance > MAX_DIST_PPM_FAR) {
            printf("%d) Warning: Lal %d is ALONE (no PPM star near it).\n",
                ++errors,
                catRef);
            logCauses(catName, nullptr, x, y, z,
                false, false, 1.0, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* la almacenamos para futuras identificaciones */
        if (countLal >= MAXLALSTAR) {
            printf("Error: too many Lalande stars.\n");
            exit(1);
        }
        lalX[countLal] = x;
        lalY[countLal] = y;
        lalZ[countLal] = z;
        lalRef[countLal] = catRef;
        countLal++;
    }
    fclose(stream);
	fclose(crossPPMStream);

    printf("Available Lalande stars = %d\n", countLal);
    printf("Stars from Lalande identified with PPM = %d\n", countDist);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("Errors logged = %d\n", errors);
}

/*
 * readStone - lee y cruza catalogo de Stone
 * (ya se deben haber leidos los catalogs CD y CPD)
 */
void readStone() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];

    /* usamos catalogs CD y CPD */
    struct DMstar_struct *CDstar = getDMStruct();
    struct CPDstar_struct *CPDstar = getCPDStruct();

    printf("\n***************************************\n");
    printf("Perform comparison between Stone and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_lacaille_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_lacaille_cpd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross/cross_lacaille_ppm.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countCD = 0;
    int countCPD = 0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_ST);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo Stone */
    FILE *stream = fopen("cat/stone1.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read stone1.txt");
		exit(1);
    }
	FILE *stream2 = fopen("cat/stone2.txt", "rt");
    if (stream == NULL) {
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
        readField(buffer, cell, 13, 1);
        bool secondStar = cell[0] == '1';

		/* lee magnitud */
		readField(buffer, cell, 24, 1);
		float vmag = atof(cell);
		readField(buffer, cell, 25, 1);
		vmag += cell[0] == '&' ? 0 : atof(cell)/10.0;

		/* lee ascension recta B1880.0 */
		readFieldSanitized(buffer, cell, 33, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 35, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 37, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1880.0 (en realidad NPD) */
		readFieldSanitized(buffer2, cell, 21, 3);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer2, cell, 24, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer2, cell, 26, 4);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/100.0)/3600.0;
		Decl = 90.0 - Decl; /* convertimos NPD a declinación */

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_PPM) {
			float ppmVmag = PPMstar[ppmIndex].vmag;
			if (vmag > __FLT_EPSILON__ && ppmVmag > __FLT_EPSILON__) {
				// Note: no fit is performed to convert scales of magnitudes
				float delta = fabs(vmag - ppmVmag);
				if (delta < MAX_MAGNITUDE) {
                    akkuDeltaError += delta * delta;
                    countDelta++;
                } else {
					printf("%d) Warning: St %d (vmag = %.1f) has a near PPM %d star (vmag = %.1f) with a different magnitude (delta = %.1f).\n",
                        ++errors,
						stoneRef,
                        vmag,
                        PPMstar[ppmIndex].ppmRef,
						ppmVmag,
						delta);	
				}
			}

            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;

            if (lacailleRef > 0) {
			    snprintf(catName, 20, "L %d", lacailleRef);
			    snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			    writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
            }
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_ST, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        if (Decl1875 <= -18.0) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            findCPDByCoordinates(x, y, z, Decl1875, &cpdIndex, &minDistance);
			if (minDistance < MAX_DIST_ST_CPD) {
                countCPD++;
                cpdFound = true;

                if (lacailleRef > 0) {
			        snprintf(catName, 20, "L %d", lacailleRef);
			        snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
			        writeCrossEntry(crossCPDStream, catName, cdName, minDistance);
                }
            }
        }
		
        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
		if (Decl1875 <= -22.0) {
            int cdIndex = -1;
            minDistance = HUGE_NUMBER;
            findDMByCoordinates(x, y, z, Decl1875, &cdIndex, &minDistance);
			if (minDistance < MAX_DIST_CD) {
			    float cdVmag = CDstar[cdIndex].vmag;
			    if (vmag > __FLT_EPSILON__ && cdVmag < 29.9) {
				    float delta = fabs(vmag - cdVmag);
				    if (delta >= MAX_MAGNITUDE) {
					    printf("%d) Warning: St %d (vmag = %.1f) has a near CD with dif. magnitude (delta = %.1f). Check dpl.\n",
                            ++errors,
                            stoneRef,
                            vmag,
                            delta);
					    writeRegister(cdIndex, false);
				    }
			    }

                countCD++;
                cdFound = true;

                if (lacailleRef > 0) {
    			    snprintf(catName, 20, "L %d", lacailleRef);
	    		    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
		    	    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
                }
            }
		}

        if (!ppmFound && !cdFound && !cpdFound) {
            if (lacailleRef > 0) {
                printf("%d) Warning: St %d (L %d) is ALONE (no PPM or CD or CPD star near it).\n",
                    ++errors,
                    stoneRef,
                    lacailleRef);
            } else {
                printf("%d) Warning: St %d is ALONE (no PPM or CD or CPD star near it).\n",
                    ++errors,
                    stoneRef);
            }
            logCauses(catName, nullptr, x, y, z,
                false, false, vmag, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* la almacenamos para futuras identificaciones */
        if (lacailleRef > 0) {
            if (countLac >= MAXLACSTAR) {
                printf("Error: too many Lac stars.\n");
                exit(1);
            }
            lacX[countLac] = x;
            lacY[countLac] = y;
            lacZ[countLac] = z;
            lacRef[countLac] = lacailleRef;
            countLac++;
        }

        if (countSt >= MAXSTSTAR && !secondStar) {
            printf("Error: too many Stone stars.\n");
            exit(1);
        }
        stX[countSt] = x;
        stY[countSt] = y;
        stZ[countSt] = z;
        stRef[countSt] = stoneRef;
        countSt++; 
    }
    fclose(stream2);
    fclose(stream);
	fclose(crossPPMStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available Stone stars = %d, and Lacaille = %d\n", countSt, countLac);
    printf("Stars from Stone identified with PPM = %d, CD = %d and CPD = %d\n", countDist, countCD, countCPD);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)countDelta),
        countDelta);
    printf("Errors logged = %d\n", errors);
}

/*
 * readTaylor - lee y cruza catalogo de Taylor
 * (ya se deben haber leidos catalogos Stone y GC)
 */
void readTaylor() {
    char buffer[1024], cell[256], catName[20], ppmName[20];

    /* usamos catalogo GC */
    struct GCstar_struct *GCstar = getGCStruct();
    int GCstars = getGCStars();

    printf("\n***************************************\n");
    printf("Perform comparison between Taylor and PPM...\n");

	FILE *crossPPMStream = openCrossFile("results/cross/cross_taylor_ppm.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int checkLac = 0;
    int checkGC = 0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_TAYLOR);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo Taylor */
    FILE *stream = fopen("cat/taylor1.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read taylor1.txt");
		exit(1);
    }
	FILE *stream2 = fopen("cat/taylor2.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read taylor2.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		char buffer2[1024];
		fgets(buffer2, 1023, stream2);

		/* lee numeración */
		readField(buffer, cell, 7, 5);
		int taylorRef = atoi(cell);

		/* lee ascension recta B1835.0 */
		readFieldSanitized(buffer, cell, 46, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 48, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 50, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1835.0 */
		readFieldSanitized(buffer2, cell, 8, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer2, cell, 10, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer2, cell, 12, 4);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/100.0)/3600.0;
		readField(buffer2, cell, 7, 1);
		if (cell[0] == '-') Decl = -Decl;

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_TAY_PPM) {
            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;

			snprintf(catName, 20, "T %d", taylorRef);
			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_TAYLOR, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        if (!ppmFound && nearestPPMDistance > MAX_DIST_PPM_FAR) {
            printf("%d) Warning: T %d is ALONE (no PPM star near it).\n",
                ++errors,
                taylorRef);
            logCauses(catName, nullptr, x, y, z,
                false, false, 1.0, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* lee referencia numerica y referencia a catalogo (que queda en "cell") */
		readField(buffer, cell, 37, 5);
        int numRefCat = atoi(cell);
		readField(buffer, cell, 31, 3);

        if (!strncmp(cell, "LAC", 3)) {
            for (int i = 0; i < countLac; i++) {
                if (lacRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, lacX[i], lacY[i], lacZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: T %d is FAR from L %d (dist = %.1f arcsec).\n",
                        ++errors,
                        taylorRef,
                        numRefCat,
                        dist);
                } else checkLac++;
            }
        }

        if (!strncmp(cell, "GOU", 3)) {
            for (int i = 0; i < GCstars; i++) {
                if (GCstar[i].gcRef != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, GCstar[i].x, GCstar[i].y, GCstar[i].z);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: T %d is FAR from GC %d (dist = %.1f arcsec). Check if it is a 1/2 star in GC.\n",
                        ++errors,
                        taylorRef,
                        numRefCat,
                        dist);
                } else checkGC++;
            }
        }

        /* la almacenamos para futuras identificaciones */
        if (countTaylor >= MAXTAYLORSTAR) {
            printf("Error: too many Taylor stars.\n");
            exit(1);
        }
        tayX[countTaylor] = x;
        tayY[countTaylor] = y;
        tayZ[countTaylor] = z;
        tayRef[countTaylor] = taylorRef;
        countTaylor++;
    }
    fclose(stream2);
    fclose(stream);
	fclose(crossPPMStream);

    printf("Available Taylor stars = %d\n", countTaylor);
    printf("Taylor stars properly identified with Lacaille = %d\n", checkLac);
    printf("Taylor stars properly identified with GC = %d\n", checkGC);
    printf("Stars from Taylor identified with PPM = %d\n", countDist);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("Errors logged = %d\n", errors);
}

/*
 * readUSNO - lee y cruza catalogo de Yarnall-Frisby
 * (ya se deben haber leidos los catalogs CD, CPD, Weiss y Stone)
 * tambien revisa referencias cruzadas a OA y Lacaille
 */
void readUSNO() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];
    char catLine[64];

    /* usamos catalogs CD y CPD */
    struct DMstar_struct *CDstar = getDMStruct();
    struct CPDstar_struct *CPDstar = getCPDStruct();

    printf("\n***************************************\n");
    printf("Perform comparison between USNO and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_usno_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_usno_cpd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross/cross_usno_ppm.csv");
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/usno_unidentified.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countCD = 0;
    int countCPD = 0;
    int checkLac = 0;
    int checkLal = 0;
    int checkOA = 0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_USNO);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

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

		/* lee ascension recta B1860.0 */
		readFieldSanitized(buffer, cell, 40, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 42, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 44, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1860.0 */
		readFieldSanitized(buffer, cell, 61, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 63, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 65, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		readField(buffer, cell, 60, 1);
		if (cell[0] == '-') Decl = -Decl;

        snprintf(catLine, 64, "%02dh %02dm %02ds%02d %c%02d°%02d'%02d\"%01d",
            RAh, RAm, RAs / 100, RAs % 100, cell[0], Decld, Declm, Decls / 10, Decls % 10);

		/* lee numeracion */
		readField(buffer, cell, 1, 5);
		int numRef = atoi(cell);

        snprintf(catName, 20, "U %d", numRef);

		/* lee magnitud */
		readField(buffer, cell, 37, 2);
		float vmag = atof(cell)/10.0;

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_PPM) {
			float ppmVmag = PPMstar[ppmIndex].vmag;
			if (vmag > __FLT_EPSILON__ && ppmVmag > __FLT_EPSILON__) {
				// Note: no fit is performed to convert scales of magnitudes
				float delta = fabs(vmag - ppmVmag);
				if (delta < MAX_MAGNITUDE) {
                    akkuDeltaError += delta * delta;
                    countDelta++;
                } else {
					printf("%d) Warning: U %d (vmag = %.1f) has a near PPM %d star (vmag = %.1f) with a different magnitude (delta = %.1f).\n",
                        ++errors,
						numRef,
                        vmag,
                        PPMstar[ppmIndex].ppmRef,
						ppmVmag,
						delta);	
				}
			}

            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;

			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_USNO, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        if (Decl1875 <= -18.0) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            findCPDByCoordinates(x, y, z, Decl1875, &cpdIndex, &minDistance);
			if (minDistance < MAX_DIST_USNO_CPD) {
                countCPD++;
                cpdFound = true;

			    snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
			    writeCrossEntry(crossCPDStream, catName, cdName, minDistance);
            }
        }
		
        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
		if (Decl1875 <= -22.0) {
            int cdIndex = -1;
            minDistance = HUGE_NUMBER;
            findDMByCoordinates(x, y, z, Decl1875, &cdIndex, &minDistance);
			if (minDistance < MAX_DIST_CD) {
			    float cdVmag = CDstar[cdIndex].vmag;
			    if (vmag > __FLT_EPSILON__ && cdVmag < 29.9) {
				    float delta = fabs(vmag - cdVmag);
				    if (delta >= MAX_MAGNITUDE) {
					    printf("%d) Warning: U %d (vmag = %.1f) has a near CD with dif. magnitude (delta = %.1f). Check dpl.\n",
                            ++errors,
                            numRef,
                            vmag,
                            delta);
					    writeRegister(cdIndex, false);
				    }
			    }

                countCD++;
                cdFound = true;

			    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
			    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
            }
		}

        /* lee referencia numerica y referencia a catalogo (que queda en "cell") */
		readField(buffer, cell, 30, 5);
        int numRefCat = atoi(cell);
		readField(buffer, cell, 19, 4);

        if (!strncmp(cell, "LAC ", 4)) {
            for (int i = 0; i < countLac; i++) {
                if (lacRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, lacX[i], lacY[i], lacZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: U %d is FAR from L %d (dist = %.1f arcsec).\n",
                        ++errors,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register U %d: %s\n", numRef, catLine);    
                } else checkLac++;
            }
        }

        if (!strncmp(cell, "LAL ", 4)) {
            for (int i = 0; i < countLal; i++) {
                if (lalRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, lalX[i], lalY[i], lalZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: U %d is FAR from Lal %d (dist = %.1f arcsec).\n",
                        ++errors,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register U %d: %s\n", numRef, catLine);    
                } else checkLal++;
            }
        }

        if (!strncmp(cell, "OARS", 4)) {
            for (int i = 0; i < countOA; i++) {
                if (oaRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, oaX[i], oaY[i], oaZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: U %d is FAR from OA %d (dist = %.1f arcsec).\n",
                        ++errors,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register U %d: %s\n", numRef, catLine);    
                } else checkOA++;
            }
        }

        /* Alone stars that are from Praesepe or Pleiades are ignored but recognized as unidentified. */
        if (!ppmFound && !cdFound && !cpdFound) {
            readField(buffer, cell, 15, 3);
            bool praesepe = !strncmp(cell, "PRA", 3);
            bool pleiades = !strncmp(cell, "PLE", 3);
            if (praesepe) {
                printf("**) U %s is from Praesepe and is ALONE.\n", catName);
            } else if (pleiades) {
                printf("**) U %s is from Pleiades and is ALONE.\n", catName);
            } else {
                printf("%d) Warning: U %d is ALONE (no PPM or CD or CPD star near it).\n",
                    ++errors,
                    numRef);
                printf("     Register U %d: %s\n", numRef, catLine);    
            }
            logCauses(catName, unidentifiedStream, x, y, z,
                false, false, vmag, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        if (countUsno >= MAXUSNOSTAR) {
            printf("Error: too many USNO stars.\n");
            exit(1);
        }
        usnoX[countUsno] = x;
        usnoY[countUsno] = y;
        usnoZ[countUsno] = z;
        usnoRef[countUsno] = numRef;
        countUsno++;
    }
    fclose(stream);
    fclose(unidentifiedStream);
	fclose(crossPPMStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available USNO stars = %d\n", countUsno);
    printf("Stars from USNO with Lacaille = %d, Lalande = %d and OA = %d\n",
        checkLac, checkLal, checkOA);
    printf("Stars from USNO identified with PPM = %d, CD = %d and CPD = %d\n", countDist, countCD, countCPD);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)countDelta),
        countDelta);
    printf("Errors logged = %d\n", errors);
}

/*
 * readUA - lee y cruza Uranometria Argentina
 * (ya se deben haber leidos los catalogs CD, CPD, Weiss, Stone y Yarnall)
 * tambien revisa referencias cruzadas; no tiene en cuenta magnitudes
 */
void readUA() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];
    char catLine[64];

    /* usamos catalogs CD y CPD */
    struct DMstar_struct *CDstar = getDMStruct();
    struct CPDstar_struct *CPDstar = getCPDStruct();

    printf("\n***************************************\n");
    printf("Perform comparison between UA and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_ua_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_ua_cpd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross/cross_ua_ppm.csv");
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/ua_unidentified.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countCD = 0;
    int countCPD = 0;
    int errors = 0;
    int checkLac = 0;
    int checkLal = 0;
    int checkOA = 0;
    int checkTaylor = 0;
    int checkUSNO = 0;
    int countUA = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_UA);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo UA */
    FILE *stream = fopen("cat/ua.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read ua.txt");
		exit(1);
    }
    char serpens = 'a';
    while (fgets(buffer, 1023, stream) != NULL) {
        /* no leemos estrellas sin designación Gould ni sin coordenadas */
		readFieldSanitized(buffer, cell, 1, 1);
        if (cell[0] != 'G') continue;
		readField(buffer, cell, 101, 1);
        if (cell[0] == ' ') continue;

		/* lee ascension recta B1875.0 */
		readFieldSanitized(buffer, cell, 100, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 103, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 106, 2);
		int RAs = atoi(cell);
        RA += ((double) RAs)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1875.0 */
		readFieldSanitized(buffer, cell, 111, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 114, 4);
		double Declm = atof(cell);
        Decl += Declm/60.0;
		readField(buffer, cell, 110, 1);
		if (cell[0] == '-') Decl = -Decl;

        snprintf(catLine, 64, "%02dh %02dm %02ds %c%02d°%02.1f'",
            RAh, RAm, RAs, cell[0], Decld, Declm);

		/* lee numeración de Gould y constelación */
		readField(buffer, cell, 3, 3);
		int gouldRef = atoi(cell);
        char cstRef[5];
		readField(buffer, cell, 7, 3);
        copyWithoutSpaces(cstRef, cell);
        if (cstRef[0] == 'S' && cstRef[1] == 'e' && cstRef[2] == 'r') {
            /* Serpens tiene parte (a) y (b) */
            cstRef[3] = serpens;
            cstRef[4] = 0;
            /* si es la ultima estrella de (a), actualiza a (b) */
            if (gouldRef == 49) {
                serpens = 'b';
            }
        }

        /* si está disponible Bayer, usa esa designacion */
        readField(buffer, cell, 15, 8);
        if (cell[0] != ' ') {
            char bayerRef[9];
            copyWithoutSpaces(bayerRef, cell);
            snprintf(catName, 20, "%s%s", bayerRef, cstRef);
        } else {
            /* si está disponible Flamsteed, la usa */
            readField(buffer, cell, 11, 3);
            if (cell[2] != ' ') {
                int fRef = atoi(cell);
                snprintf(catName, 20, "%d %s", fRef, cstRef);
            } else {
                /* caso contrario, usa la denominación de Gould */
                snprintf(catName, 20, "%dG %s", gouldRef, cstRef);
            }
        }

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_PPM) {
            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;

			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		}

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        if (Decl <= -18.0) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            findCPDByCoordinates(x, y, z, Decl, &cpdIndex, &minDistance);
			if (minDistance < MAX_DIST_CPD) {
                countCPD++;
                cpdFound = true;

			    snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
			    writeCrossEntry(crossCPDStream, catName, cdName, minDistance);
            }
        }
		
        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
		if (Decl <= -22.0) {
            int cdIndex = -1;
            minDistance = HUGE_NUMBER;
            findDMByCoordinates(x, y, z, Decl, &cdIndex, &minDistance);
			if (minDistance < MAX_DIST_CD) {
                countCD++;
                cdFound = true;

			    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
			    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
            }
		}

        /* lee referencia a catalogo */
        char desigRefCat[4];
		readField(buffer, desigRefCat, 82, 3);

        if (!strncmp(desigRefCat, "L.", 2) && desigRefCat[2] != '(') {
            readField(buffer, cell, 84, 4);
            char numStr[5];
            sscanf(cell, "%4[^, ]", numStr);
            int numRefCat = atoi(numStr);
            for (int i = 0; i < countLac; i++) {
                if (lacRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, lacX[i], lacY[i], lacZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: %dG %s is FAR from L %d (dist = %.1f arcsec).\n",
                        ++errors,
                        gouldRef,
                        cstRef,
                        numRefCat,
                        dist);
                    printf("     Register %dG %s: %s\n", gouldRef, cstRef, catLine);    
                } else checkLac++;
            }
        }

        if (!strncmp(desigRefCat, "T.", 2)) {
            readField(buffer, cell, 84, 5);
            char numStr[6];
            sscanf(cell, "%5[^, ]", numStr);
            int numRefCat = atoi(numStr);
            for (int i = 0; i < countTaylor; i++) {
                if (tayRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, tayX[i], tayY[i], tayZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: %dG %s is FAR from T %d (dist = %.1f arcsec).\n",
                        ++errors,
                        gouldRef,
                        cstRef,
                        numRefCat,
                        dist);
                    printf("     Register %dG %s: %s\n", gouldRef, cstRef, catLine);    
                } else checkTaylor++;
            }
        }

        if (!strncmp(desigRefCat, "Ll.", 3)) {
            readField(buffer, cell, 85, 5);
            char numStr[6];
            sscanf(cell, "%5[^, ]", numStr);
            int numRefCat = atoi(numStr);
            for (int i = 0; i < countLal; i++) {
                if (lalRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, lalX[i], lalY[i], lalZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: %dG %s is FAR from Lal %d (dist = %.1f arcsec).\n",
                        ++errors,
                        gouldRef,
                        cstRef,
                        numRefCat,
                        dist);
                    printf("     Register %dG %s: %s\n", gouldRef, cstRef, catLine);    
                } else checkLal++;
            }
        }

        if (!strncmp(desigRefCat, "OA.", 3)) {
            readField(buffer, cell, 85, 5);
            char numStr[6];
            sscanf(cell, "%5[^, ]", numStr);
            int numRefCat = atoi(numStr);
            for (int i = 0; i < countOA; i++) {
                if (oaRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, oaX[i], oaY[i], oaZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: %dG %s is FAR from OA %d (dist = %.1f arcsec).\n",
                        ++errors,
                        gouldRef,
                        cstRef,
                        numRefCat,
                        dist);
                    printf("     Register %dG %s: %s\n", gouldRef, cstRef, catLine);    
                } else checkOA++;
            }
        }

        if (!strncmp(desigRefCat, "Y.", 2)) {
            readField(buffer, cell, 84, 5);
            char numStr[6];
            sscanf(cell, "%5[^, ]", numStr);
            int numRefCat = atoi(numStr);
            /* Here, it is different from other cross-identifications.
             * Difference between catalogs Yarnall (USNO 2nd Edition) and Yarnall-Frisby
             * (USNO 3rd Edition) are in enumeration, which is similar but not equal. */
            minDistance = HUGE_NUMBER;
            int usnoIndex = -1;
            for (int i = 0; i < countUsno; i++) {
                double dist = 3600.0 * calcAngularDistance(x, y, z, usnoX[i], usnoY[i], usnoZ[i]);
                if (minDistance > dist) {
                    usnoIndex = i;
                    minDistance = dist;
                }
            }
            if (minDistance > MAX_DIST_CROSS) {
                printf("%d) Warning: %dG %s is FAR from Y %d / U %d (dist = %.1f arcsec).\n",
                    ++errors,
                    gouldRef,
                    cstRef,
                    numRefCat,
                    usnoRef[usnoIndex],
                    minDistance);
                printf("     Register %dG %s: %s\n", gouldRef, cstRef, catLine);    
            } else checkUSNO++;
        }
        
        if (!ppmFound && !cdFound && !cpdFound) {
            printf("%d) Warning: %dG %s is ALONE (no PPM or CD or CPD star near it).\n",
                ++errors,
                gouldRef,
                cstRef);
            printf("     Register %dG %s: %s\n", gouldRef, cstRef, catLine);    
            readField(buffer, cell, 132, 3);
            bool cumulus = !strncmp(cell, "cum", 3);
            bool nebula = !strncmp(cell, "neb", 3);
            logCauses(catName, unidentifiedStream, x, y, z,
                cumulus, nebula, 1.0, RAs, Decl, 1,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }
        countUA++;
    }
    fclose(stream);
    fclose(unidentifiedStream);
	fclose(crossPPMStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Available UA stars = %d\n", countUA);
    printf("Stars from UA with Lacaille = %d, Lalande = %d, Taylor = %d, OA = %d and USNO = %d\n",
        checkLac, checkLal, checkTaylor, checkOA, checkUSNO);
    printf("Stars from UA identified with PPM = %d, CD = %d and CPD = %d\n", countDist, countCD, countCPD);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("Errors logged = %d\n", errors);
}

/*
 * readThome - lee y cruza catalogo de los Resultados 15
 * (ya se deben haber leidos los catalogs CD, CPD, GC, Yarnall y Stone)
 * tambien revisa referencias cruzadas a GC, OA, Yarnall, Lacaille y Stone
 * tambien genera identificaciones cruzadas para aquellas estrellas ZC
 */
void readThome(double epoch, const char *filename, int correction) {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];
    char catLine[64];

    /* usamos catalogs CD, CPD y GC */
    struct DMstar_struct *CDstar = getDMStruct();
    struct CPDstar_struct *CPDstar = getCPDStruct();
    struct GCstar_struct *GCstar = getGCStruct();
    int GCstars = getGCStars();

    printf("\n***************************************\n");
    printf("Perform comparison between Thome (Resultados ONA 15) and PPM/CD/CPD...\n");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countCD = 0;
    int countCPD = 0;
    int checkLac = 0;
    int checkLal = 0;
    int checkTaylor = 0;
    int checkOA = 0;
    int checkStone = 0;
    int checkGC = 0;
    int checkUSNO = 0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, epoch);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo Thome */
    FILE *stream = fopen(filename, "rt");
    if (stream == NULL) {
        perror("Cannot read Thome catalog");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		/* lee ascension recta */
		readFieldSanitized(buffer, cell, 20+correction, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 22+correction, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 24+correction, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion */
		readFieldSanitized(buffer, cell, 45+correction, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 47+correction, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 49+correction, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */

        snprintf(catLine, 64, "%02dh %02dm %02ds%02d -%02d°%02d'%02d\"%01d",
            RAh, RAm, RAs / 100, RAs % 100, Decld, Declm, Decls / 10, Decls % 10);

		/* lee numeración catálogo y magnitud */
		readField(buffer, cell, 1, 4+correction);
		int numRef = atoi(cell);
		readField(buffer, cell, 13+correction, 3);
        if (cell[2] == ' ') cell[2] = '0';
		float vmag = atof(cell)/10.0;

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_PPM) {
			float ppmVmag = PPMstar[ppmIndex].vmag;
			if (vmag > __FLT_EPSILON__ && ppmVmag > __FLT_EPSILON__) {
				// Note: no fit is performed to convert scales of magnitudes
				float delta = fabs(vmag - ppmVmag);
				if (delta < MAX_MAGNITUDE) {
                    akkuDeltaError += delta * delta;
                    countDelta++;
                } else {
					printf("%d) Warning: T %.0f %d (vmag = %.1f) has a near PPM %d star (vmag = %.1f) with a different magnitude (delta = %.1f).\n",
                        ++errors,
                        epoch,
						numRef,
                        vmag,
                        PPMstar[ppmIndex].ppmRef,
						ppmVmag,
						delta);	
				}
			}

            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(epoch, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        int cpdIndex = -1;
        minDistance = HUGE_NUMBER;
        findCPDByCoordinates(x, y, z, Decl1875, &cpdIndex, &minDistance);
        double nearestCPDDistance = minDistance;
		if (minDistance < MAX_DIST_TH_CPD) {
            countCPD++;
            cpdFound = true;
        }
		
        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
        int cdIndex = -1;
        minDistance = HUGE_NUMBER;
        findDMByCoordinates(x, y, z, Decl1875, &cdIndex, &minDistance);
        double nearestCDDistance = minDistance;
		if (minDistance < MAX_DIST_CD) {
		    float cdVmag = CDstar[cdIndex].vmag;
		    if (vmag > __FLT_EPSILON__ && cdVmag < 29.9) {
			    float delta = fabs(vmag - cdVmag);
			    if (delta >= MAX_MAGNITUDE) {
				    printf("%d) Warning: T %.0f %d (vmag = %.1f) has a near CD with dif. magnitude (delta = %.1f). Check dpl.\n",
                        ++errors,
                        epoch,
                        numRef,
                        vmag,
                        delta);
				    writeRegister(cdIndex, false);
			    }
		    }

            countCD++;
            cdFound = true;
		}

        if (!ppmFound && !cdFound && !cpdFound) {
            printf("%d) Warning: T %.0f %d is ALONE (no PPM or CD or CPD star near it).\n",
                ++errors,
                epoch,
                numRef);
            printf("     Register T %.0f %d: %s\n", epoch, numRef, catLine);
            logCauses(catName, nullptr, x, y, z,
                false, false, vmag, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* lee referencia numerica y referencia a catalogo (que queda en "cell") */
		readField(buffer, cell, 68+correction, 5);
        int numRefCat = atoi(cell);
		readField(buffer, cell, 66+correction, 2);

        if (!strncmp(cell, "L ", 2)) {
            for (int i = 0; i < countLac; i++) {
                if (lacRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, lacX[i], lacY[i], lacZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: T %.0f %d is FAR from L %d (dist = %.1f arcsec).\n",
                        ++errors,
                        epoch,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register T %.0f %d: %s\n", epoch, numRef, catLine);
                } else checkLac++;
            }
        }

        if (!strncmp(cell, "T ", 2)) {
            for (int i = 0; i < countTaylor; i++) {
                if (tayRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, tayX[i], tayY[i], tayZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: T %.0f %d is FAR from T %d (dist = %.1f arcsec).\n",
                        ++errors,
                        epoch,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register T %.0f %d: %s\n", epoch, numRef, catLine);
                } else checkTaylor++;
            }
        }

        if (!strncmp(cell, "LL", 2)) {
            for (int i = 0; i < countLal; i++) {
                if (lalRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, lalX[i], lalY[i], lalZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: T %.0f %d is FAR from Lal %d (dist = %.1f arcsec).\n",
                        ++errors,
                        epoch,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register T %.0f %d: %s\n", epoch, numRef, catLine);
                } else checkLal++;
            }
        }

        if (!strncmp(cell, "ST", 2)) {
            for (int i = 0; i < countSt; i++) {
                if (stRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, stX[i], stY[i], stZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: T %.0f %d is FAR from St %d (dist = %.1f arcsec).\n",
                        ++errors,
                        epoch,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register T %.0f %d: %s\n", epoch, numRef, catLine);
                } else checkStone++;
            }
        }

        if (!strncmp(cell, "OA", 2)) {
            for (int i = 0; i < countOA; i++) {
                if (oaRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, oaX[i], oaY[i], oaZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: T %.0f %d is FAR from OA %d (dist = %.1f arcsec).\n",
                        ++errors,
                        epoch,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register T %.0f %d: %s\n", epoch, numRef, catLine);
                } else checkOA++;
            }
        }

        if (!strncmp(cell, "GC", 2)) {
            for (int i = 0; i < GCstars; i++) {
                if (GCstar[i].gcRef != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, GCstar[i].x, GCstar[i].y, GCstar[i].z);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: T %.0f %d is FAR from GC %d (dist = %.1f arcsec). Check if it is a 1/2 star in GC.\n",
                        ++errors,
                        epoch,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register T %.0f %d: %s\n", epoch, numRef, catLine);
                } else checkGC++;
            }
        }

        if (!strncmp(cell, "Y ", 2)) {
            /* Here, it is different from other cross-identifications.
             * Difference between catalogs Yarnall (USNO 2nd Edition) and Yarnall-Frisby
             * (USNO 3rd Edition) are in enumeration, which is similar but not equal. */
            minDistance = HUGE_NUMBER;
            int usnoIndex = -1;
            for (int i = 0; i < countUsno; i++) {
                double dist = 3600.0 * calcAngularDistance(x, y, z, usnoX[i], usnoY[i], usnoZ[i]);
                if (minDistance > dist) {
                    usnoIndex = i;
                    minDistance = dist;
                }
            }
            if (minDistance > MAX_DIST_CROSS) {
                printf("%d) Warning: T %.0f %d is FAR from Y %d / U %d (dist = %.1f arcsec).\n",
                    ++errors,
                    epoch,
                    numRef,
                    numRefCat,
                    usnoRef[usnoIndex],
                    minDistance);
                printf("     Register T %.0f %d: %s\n", epoch, numRef, catLine);
            } else checkUSNO++;
        }

        if (!strncmp(cell, "ZC", 2)) {
            // save Gould's Zone Catalog
            RAh = (int) (floor(RA1875 / 15.0) + __FLT_EPSILON__);
            saveZC(RAh, RAs, Decls, numRefCat, x, y, z,
                ppmFound, ppmIndex, nearestPPMDistance,
                cdFound, cdIndex, nearestCDDistance,
                cpdFound, cpdIndex, nearestCPDDistance);
        }
    }
    fclose(stream);

    printf("Stars from Thome %.0f identified with PPM = %d, CD = %d and CPD = %d\n", epoch, countDist, countCD, countCPD);
    printf("Stars from Thome %.0f with Lacaille = %d, Lalande = %d, OA = %d, Taylor = %d, Stone = %d, GC = %d and USNO = %d\n",
        epoch, checkLac, checkLal, checkOA, checkTaylor, checkStone, checkGC, checkUSNO);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)countDelta),
        countDelta);
    printf("Errors logged = %d\n", errors);
}

/*
 * readGilliss - lee y cruza catalogo de Gilliss
 * (ya se deben haber leidos los catalogs CD, CPD, GC y Stone)
 * tambien revisa referencias cruzadas a GC, Stone y Lacaille
 * tambien genera identificaciones cruzadas para aquellas estrellas ZC
 */
void readGilliss() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];
    char catLine[64];

    /* usamos catalogs CD, CPD y GC */
    struct DMstar_struct *CDstar = getDMStruct();
    struct CPDstar_struct *CPDstar = getCPDStruct();
    struct GCstar_struct *GCstar = getGCStruct();
    int GCstars = getGCStars();

    printf("\n***************************************\n");
    printf("Perform comparison between Gilliss and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_gilliss_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_gilliss_cpd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross/cross_gilliss_ppm.csv");
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/gilliss_unidentified.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countCD = 0;
    int countCPD = 0;
    int checkLac = 0;
    int checkGC = 0;
    int checkStone = 0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_GILLISS);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

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
		readFieldSanitized(buffer, cell, 33, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 35, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 37, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1850.0 */
		readFieldSanitized(buffer, cell, 49, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 51, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 53, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */
        
        snprintf(catLine, 64, "%02dh %02dm %02ds%02d -%02d°%02d'%02d\"%01d",
            RAh, RAm, RAs / 100, RAs % 100, Decld, Declm, Decls / 10, Decls % 10);

		/* lee numeración catálogo Gilliss */
		readField(buffer, cell, 1, 5);
		int giRef = atoi(cell);

        /* si esta disponible, tambien lee precesiones y chequea (usamos constantes de Struve) */
        readFieldSanitized(buffer, cell, 41, 7);
        double preRA = atof(cell) / 1000.0;
        if (fabs(preRA) > EPS) {
            double realPreRA = 3.07177 + 1.3371 * dsin(RA) * dtan(Decl);
            double diff = fabs(preRA - realPreRA);
            if (diff > 0.0299) {
                printf("%d) Warning: G %d reports %.3f on precession RA but it should be %.3f (diff=%.3f).\n",
                    ++errors,
                    giRef,
                    preRA,
                    realPreRA,
                    diff);
                printf("     Register G %d: %s\n", giRef, catLine);    
            }
        }

        readFieldSanitized(buffer, cell, 56, 5);
        double preDecl = atof(cell) / 100.0;
        if (fabs(preDecl) > EPS) {
            double realPreDecl = 20.0564 * dcos(RA);
            double diff = fabs(preDecl - realPreDecl);
            if (diff > 0.099) {
                printf("%d) Warning: G %d reports %.2f on precession DECL but it should be %.2f (diff=%.2f).\n",
                    ++errors,
                    giRef,
                    preDecl,
                    realPreDecl,
                    diff);
                printf("     Register G %d: %s\n", giRef, catLine);    
            }
        }

        /* lee magnitud */
		readField(buffer, cell, 29, 3);
		float vmag = atof(cell)/10.0;

        snprintf(catName, 20, "G %d", giRef);

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_PPM) {
			float ppmVmag = PPMstar[ppmIndex].vmag;
			if (vmag > __FLT_EPSILON__ && ppmVmag > __FLT_EPSILON__) {
				// Note: no fit is performed to convert scales of magnitudes
				float delta = fabs(vmag - ppmVmag);
				if (delta < MAX_MAGNITUDE) {
                    akkuDeltaError += delta * delta;
                    countDelta++;
                } else {
					printf("%d) Warning: G %d (vmag = %.1f) has a near PPM %d star (vmag = %.1f) with a different magnitude (delta = %.1f).\n",
                        ++errors,
						giRef,
                        vmag,
                        PPMstar[ppmIndex].ppmRef,
						ppmVmag,
						delta);	
				}
			}

            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;

			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(EPOCH_GILLISS, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        int cpdIndex = -1;
        minDistance = HUGE_NUMBER;
        findCPDByCoordinates(x, y, z, Decl1875, &cpdIndex, &minDistance);
        double nearestCPDDistance = minDistance;
		if (minDistance < MAX_DIST_CPD) {
            countCPD++;
            cpdFound = true;

		    snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
		    writeCrossEntry(crossCPDStream, catName, cdName, minDistance);
        }
		
        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
        int cdIndex = -1;
        minDistance = HUGE_NUMBER;
        findDMByCoordinates(x, y, z, Decl1875, &cdIndex, &minDistance);
        double nearestCDDistance = minDistance;
		if (minDistance < MAX_DIST_CD) {
		    float cdVmag = CDstar[cdIndex].vmag;
		    if (vmag > __FLT_EPSILON__ && cdVmag < 29.9) {
			    float delta = fabs(vmag - cdVmag);
			    if (delta >= MAX_MAGNITUDE) {
				    printf("%d) Warning: G %d (vmag = %.1f) has a near CD with dif. magnitude (delta = %.1f). Check dpl.\n",
                        ++errors,
                        giRef,
                        vmag,
                        delta);
				    writeRegister(cdIndex, false);
			    }
		    }

            countCD++;
            cdFound = true;

		    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
		    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
		}

        if (!ppmFound && !cdFound && !cpdFound) {
            printf("%d) Warning: G %d is ALONE (no PPM or CD or CPD star near it).\n",
                ++errors,
                giRef);
            printf("     Register G %d: %s\n", giRef, catLine);
            logCauses(catName, unidentifiedStream, x, y, z,
                false, false, vmag, RAs, Decl1875, Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
        }

        /* lee referencia numerica y referencia a catalogo (que queda en "cell") */
		readField(buffer, cell, 20, 5);
        int numRefCat = atoi(cell);
		readField(buffer, cell, 15, 3);

        if (!strncmp(cell, "LAC", 3)) {
            for (int i = 0; i < countLac; i++) {
                if (lacRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, lacX[i], lacY[i], lacZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: G %d is FAR from L %d (dist = %.1f arcsec).\n",
                        ++errors,
                        giRef,
                        numRefCat,
                        dist);
                    printf("     Register G %d: %s\n", giRef, catLine);
                } else checkLac++;
            }
        }

        if (!strncmp(cell, "STO", 3)) {
            for (int i = 0; i < countSt; i++) {
                if (stRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, stX[i], stY[i], stZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: G %d is FAR from St %d (dist = %.1f arcsec).\n",
                        ++errors,
                        giRef,
                        numRefCat,
                        dist);
                    printf("     Register G %d: %s\n", giRef, catLine);
                } else checkStone++;
            }
        }

        if (!strncmp(cell, "GOU", 3)) {
            for (int i = 0; i < GCstars; i++) {
                if (GCstar[i].gcRef != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, GCstar[i].x, GCstar[i].y, GCstar[i].z);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: G %d is FAR from GC %d (dist = %.1f arcsec). Check if it is a 1/2 star in GC.\n",
                        ++errors,
                        giRef,
                        numRefCat,
                        dist);
                    printf("     Register G %d: %s\n", giRef, catLine);
                } else checkGC++;
            }
        }

        if (!strncmp(cell, "GZC", 3)) {
            // save Gould's Zone Catalog
            readField(buffer, cell, 18, 2);
            RAh = atoi(cell);
            int originalRAh = (int) floor(RA1875/15.0 + __FLT_EPSILON__);
            if (RAh != originalRAh) {
                printf("%d) Warning: G %d has a different RA zone (%d) than identified ZC (%d).\n",
                    ++errors,
                    giRef,
                    RAh,
                    originalRAh);
            }
            saveZC(RAh, RAs, Decls, numRefCat, x, y, z,
                ppmFound, ppmIndex, nearestPPMDistance,
                cdFound, cdIndex, nearestCDDistance,
                cpdFound, cpdIndex, nearestCPDDistance);
        }
    }
    fclose(stream);
    fclose(unidentifiedStream);
	fclose(crossPPMStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Stars from Gilliss identified with PPM = %d, CD = %d and CPD = %d\n", countDist, countCD, countCPD);
    printf("Stars from Gilliss with Lacaille = %d, Stone = %d and GC = %d\n",
        checkLac, checkStone, checkGC);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)countDelta),
        countDelta);
    printf("Errors logged = %d\n", errors);
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

    /* leemos y cruzamos Lalande (solo contra PPM) */
    readLalande();

    /* leemos y cruzamos Stone / Lacaille */
    readStone();

    /* leemos, cruzamos y revisamos identificaciones de Yarnall */
    readUSNO();

    /* leemos catalogo GC */
	readGC();
	struct GCstar_struct *GCstar = getGCStruct();
	int GCstars = getGCStars();

    /* leemos y cruzamos Taylor */
    readTaylor();

    /* leemos, cruzamos y revisamos Uranometria Argentina */
    readUA();

    /* leemos, cruzamos y revisamos identificaciones de Thome */
    /* tambien generamos Gould's Zone Catalog */
    crossPPMZCStream = openCrossFile("results/cross/cross_zc_ppm.csv");
    crossCDZCStream = openCrossFile("results/cross/cross_zc_cd.csv");
    crossCPDZCStream = openCrossFile("results/cross/cross_zc_cpd.csv");
    unidentifiedZCStream = openUnidentifiedFile("results/cross/zc_unidentified.csv");

    readThome(1881.0, "cat/thome1881.txt", 0);
    readThome(1882.0, "cat/thome1882.txt", 0);
    readThome(1883.0, "cat/thome1883.txt", -1);
    readThome(1884.0, "cat/thome1884.txt", -1);

    /* leemos, cruzamos y revisamos identificaciones de Gilliss */
    /* tambien generamos Gould's Zone Catalog */
    readGilliss();

    fclose(unidentifiedZCStream);
	fclose(crossPPMZCStream);
	fclose(crossCPDZCStream);
	fclose(crossCDZCStream);
    printf("\nNumber of ZC stars registered = %d\n", countZC);
    printf("Stars from ZC identified with PPM = %d, CD = %d and CPD = %d\n", countPPMZC, countCDZC, countCPDZC);
    return 0;
}
