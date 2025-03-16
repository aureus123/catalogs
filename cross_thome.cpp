
/*
 * CROSS_THOME - Compara registros de varios catalogos
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include "read_ppm.h"
#include "read_cpd.h"
#include "read_dm.h"
#include "read_gc.h"
#include "trig.h"
#include "misc.h"

#define MAXOASTAR 19000
#define MAXLACSTAR 19000
#define MAX_DIST_OA_CAT 30.0 // Threshold for OA stars
#define MAX_DIST_ST_CAT 20.0 // Threshold for Stone stars
#define CURATED true // true if curated CD catalog should be used

// Here, we save 1875.0 coordinates of OA stars in rectangular form
double oaX[MAXOASTAR], oaY[MAXOASTAR], oaZ[MAXOASTAR];
int oaRef[MAXOASTAR];
int countOA = 0;

// Here, we save 1875.0 coordinates of Lac stars in rectangular form
double lacX[MAXLACSTAR], lacY[MAXLACSTAR], lacZ[MAXLACSTAR];
int lacRef[MAXLACSTAR];
int countLac = 0;

/*
 * readWeiss - lee y cruza catalogo de Weiss
 * (ya se deben haber leidos los catalogs CD y CPD)
 */
void readWeiss() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];

    /* usamos catalogs CD y CPD */
    struct DMstar_struct *CDstar = getDMStruct();
    int CDstars = getDMStars();
    struct CPDstar_struct *CPDstar = getCPDStruct();
    int CPDstars = getCPDStars();

    printf("\n***************************************\n");
    printf("Perform comparison between Weiss and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross_oa_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross_oa_cpd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross_oa_ppm.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countCD = 0;
    int countCPD = 0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, true, false, 1850.0);
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
		findPPMByCoordinates(x, y, z, &ppmIndex, &minDistance);
		if (minDistance < MAX_DIST_OA_CAT) {
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

			snprintf(catName, 20, "OA %d", oeltzenRef);
			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(1850.0, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        if (Decl1875 <= -18.0) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            for (int i = 0; i < CPDstars; i++) {
			    double dist = 3600.0 * calcAngularDistance(x, y, z, CPDstar[i].x, CPDstar[i].y, CPDstar[i].z);
			    if (minDistance > dist) {
				    cpdIndex = i;
				    minDistance = dist;
			    }
            }
			if (minDistance < MAX_DIST_OA_CAT) {
                countCPD++;
                cpdFound = true;

			    snprintf(catName, 20, "OA %d", oeltzenRef);
			    snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
			    writeCrossEntry(crossCPDStream, catName, cdName, minDistance);
            }
        }
		
        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
		if (Decl1875 <= -22.0) {
            int cdIndex = -1;
            minDistance = HUGE_NUMBER;
            for (int i = 0; i < CDstars; i++) {
			    double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
			    if (minDistance > dist) {
				    cdIndex = i;
				    minDistance = dist;
			    }
            }
			if (minDistance < MAX_DIST_CD) {
			    float cdVmag = CDstar[cdIndex].vmag;
			    if (vmag > __FLT_EPSILON__ && cdVmag < 29.9) {
				    float delta = fabs(vmag - cdVmag);
				    if (delta >= MAX_MAGNITUDE) {
					    printf("%d) Warning: W %d (vmag = %.1f) has a near CD with a different magnitude (delta = %.1f).\n",
                            ++errors,
                            weissRef,
                            vmag,
                            delta);
					    writeRegister(cdIndex, false);
				    }
			    }

                countCD++;
                cdFound = true;

			    snprintf(catName, 20, "OA %d", oeltzenRef);
			    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
			    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
            }
		}

        if (!ppmFound && !cdFound && !cpdFound) {
            printf("%d) Warning: W %d (OA %d) is ALONE (no PPM or CD or CPD star near it).\n",
                ++errors,
                weissRef,
                oeltzenRef);
            if (vmag > 9.9) {
                printf("  Possible cause: dim star.\n");
            }
            if (vmag < 0.1) {
                printf("  Possible cause: no magnitude (cumulus?).\n");
            }
            if (Decld > -18) {
                printf("  Note: no CD/CPD coverage for stars below 18°.\n");
            }
            if (RAs == 0) {
                printf("  Possible cause: lack of RA (s).\n");
            }
            if (Decls == 0) {
                printf("  Possible cause: lack of Decl (s).\n");
            }
        }

        /* la almacenamos para futuras identificaciones */
        if (countOA >= MAXOASTAR) {
            printf("Error: too many OA stars.\n");
            exit(1);
        }
        oaX[countOA] = x;
        oaY[countOA] = y;
        oaZ[countOA] = z;
        oaRef[countOA] = weissRef;
        countOA++;
    }
	fclose(crossPPMStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

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
 * readStone - lee y cruza catalogo de Stone
 * (ya se deben haber leidos los catalogs CD y CPD)
 */
void readStone() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];

    /* usamos catalogs CD y CPD */
    struct DMstar_struct *CDstar = getDMStruct();
    int CDstars = getDMStars();
    struct CPDstar_struct *CPDstar = getCPDStruct();
    int CPDstars = getCPDStars();

    printf("\n***************************************\n");
    printf("Perform comparison between Stone and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross_lacaille_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross_lacaille_cpd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross_lacaille_ppm.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countCD = 0;
    int countCPD = 0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, true, false, 1880.0);
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo Stone (pero solo nos interesa su identificación Lacaille) */
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
		int lacailleRef = cell[0] == '&' ? -1 : atoi(cell);
        if (lacailleRef <= 0) continue;

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
		findPPMByCoordinates(x, y, z, &ppmIndex, &minDistance);
		if (minDistance < MAX_DIST_ST_CAT) {
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

			snprintf(catName, 20, "L %d", lacailleRef);
			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		}

	    /* convierte coordenadas a 1875.0 y calcula rectangulares */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(1880.0, 1875.0, &RA1875, &Decl1875);
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        if (Decl1875 <= -18.0) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            for (int i = 0; i < CPDstars; i++) {
			    double dist = 3600.0 * calcAngularDistance(x, y, z, CPDstar[i].x, CPDstar[i].y, CPDstar[i].z);
			    if (minDistance > dist) {
				    cpdIndex = i;
				    minDistance = dist;
			    }
            }
			if (minDistance < MAX_DIST_ST_CAT) {
                countCPD++;
                cpdFound = true;

			    snprintf(catName, 20, "L %d", lacailleRef);
			    snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
			    writeCrossEntry(crossCPDStream, catName, cdName, minDistance);
            }
        }
		
        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
		if (Decl1875 <= -22.0) {
            int cdIndex = -1;
            minDistance = HUGE_NUMBER;
            for (int i = 0; i < CDstars; i++) {
			    double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
			    if (minDistance > dist) {
				    cdIndex = i;
				    minDistance = dist;
			    }
            }
			if (minDistance < MAX_DIST_CD) {
			    float cdVmag = CDstar[cdIndex].vmag;
			    if (vmag > __FLT_EPSILON__ && cdVmag < 29.9) {
				    float delta = fabs(vmag - cdVmag);
				    if (delta >= MAX_MAGNITUDE) {
					    printf("%d) Warning: St %d (vmag = %.1f) has a near CD with a different magnitude (delta = %.1f).\n",
                            ++errors,
                            stoneRef,
                            vmag,
                            delta);
					    writeRegister(cdIndex, false);
				    }
			    }

                countCD++;
                cdFound = true;

			    snprintf(catName, 20, "L %d", lacailleRef);
			    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
			    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
            }
		}

        if (!ppmFound && !cdFound && !cpdFound) {
            printf("%d) Warning: St %d (L %d) is ALONE (no PPM or CD or CPD star near it).\n",
                ++errors,
                stoneRef,
                lacailleRef);
            if (vmag < 0.1) {
                printf("  Possible cause: no magnitude (cumulus?).\n");
            }
            if (Decld > -18) {
                printf("  Note: no CD/CPD coverage for stars below 18°.\n");
            }
        }

        /* la almacenamos para futuras identificaciones */
        if (countLac >= MAXLACSTAR) {
            printf("Error: too many Stone stars.\n");
            exit(1);
        }
        lacX[countLac] = x;
        lacY[countLac] = y;
        lacZ[countLac] = z;
        lacRef[countLac] = lacailleRef;
        countLac++;
    }
	fclose(crossPPMStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

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
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("CROSS_THOME - Compare several catalogs.\n");
    printf("Made in 2025 by Daniel Severin.\n");

    /* leemos catalogo CD */
    readDM(CURATED ? "cat/cd_curated.txt" : "cat/cd.txt");

    /** leemos catalogo CPD */
    readCPD(false, false);

    /* leemos y cruzamos Weiss / OA */
    readWeiss();

    /* leemos y cruzamos Stone / Lacaille */
    readStone();

    /* leemos catalogo GC */
	//readGC();
	//struct GCstar_struct *GCstar = getGCStruct();
	//int GCstars = getGCStars();

    /* leemos, cruzamos y revisamos identificaciones de Yarnall */
    // readUSNO(); // TODO

    /* leemos, cruzamos y revisamos identificaciones de Gilliss */
    /* tambien generamos Gould's Zone Catalog */
    // readGilliss(); // TODO

    /* leemos, cruzamos y revisamos identificaciones de Thome */
    /* tambien generamos Gould's Zone Catalog */
    // readThome(1881.0, "cat/thome1881.txt"); // TODO
    // readThome(1882.0, "cat/thome1882.txt"); // TODO
    // readThome(1883.0, "cat/thome1883.txt"); // TODO
    // readThome(1884.0, "cat/thome1884.txt"); // TODO

    return 0;
}
