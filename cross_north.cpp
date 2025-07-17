
/*
 * CROSS_NORTH - Compara registros de varios catalogos
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "read_ppm.h"
#include "read_dm.h"
#include "trig.h"
#include "misc.h"

#define MAXWBSTAR 31900
#define MAXOASTAR 26500
#define MAXBACSTAR 8400
#define EPOCH_WB 1825.0
#define EPOCH_OA 1842.0
#define EPOCH_BAC 1850.0
#define EPOCH_USNO 1860.0
#define EPOCH_UA 1875.0
#define MAX_DIST_CROSS 60.0
#define CURATED true // true if curated BD catalog should be used

// Here, we save 1825.0 coordinates of WB stars in rectangular form
double wbX[MAXWBSTAR], wbY[MAXWBSTAR], wbZ[MAXWBSTAR];
int wbRARef[MAXWBSTAR], wbNumRef[MAXWBSTAR];
int countWB = 0;

// Here, we save 1842.0 coordinates of OA stars in rectangular form
double oaX[MAXOASTAR], oaY[MAXOASTAR], oaZ[MAXOASTAR];
int oaRef[MAXOASTAR];
int countOA = 0;

// Here, we save 1850.0 coordinates of BAC stars in rectangular form
double bacX[MAXBACSTAR], bacY[MAXBACSTAR], bacZ[MAXBACSTAR];
int bacRef[MAXBACSTAR];
int countBAC = 0;

/*
 * readWB - lee y cruza catalogo de Weisse
 */
void readWB() {
    char buffer[1024], cell[256];

    printf("\n***************************************\n");
    printf("Perform comparison between Weisse catalog and PPM...\n");

    int countDist = 0;
    double akkuDistError = 0.0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_WB);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo WB */
    FILE *stream = fopen("cat/wb.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read wb.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		/* lee numeración */
		readField(buffer, cell, 14, 2);
		int RARef = atoi(cell);
		readField(buffer, cell, 6, 4);
		int numRef = atoi(cell);

		/* lee ascension recta B1825.0 */
        double RA = (double) RARef; // same as RAh
		readFieldSanitized(buffer, cell, 16, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 18, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1825.0 */
		readFieldSanitized(buffer, cell, 31, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 33, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 35, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
        readField(buffer, cell, 30, 1);
		if (cell[0] == '-') Decl = -Decl;

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
		if (minDistance < MAX_DIST_CROSS) {
            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;
		}

        if (!ppmFound) {
            printf("%d) Warning: WB %dh %d is ALONE (nearest PPM star at %.1f arcsec).\n",
                ++errors,
                RARef,
                numRef,
                minDistance);
        }

        /* la almacenamos para futuras identificaciones */
        if (countWB >= MAXWBSTAR) {
            printf("Error: too many WB stars.\n");
            exit(1);
        }
        wbX[countWB] = x;
        wbY[countWB] = y;
        wbZ[countWB] = z;
        wbRARef[countWB] = RARef;
        wbNumRef[countWB] = numRef;
        countWB++;
    }
	fclose(stream);

    printf("Available WB stars = %d\n", countWB);
    printf("Stars from WB identified with PPM = %d\n", countDist);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("Errors logged = %d\n", errors);
}

/*
 * readOARN - lee y cruza catalogo de Oeltzen-Argelander (North)
 */
void readOARN() {
    char buffer[1024], cell[256];

    printf("\n***************************************\n");
    printf("Perform comparison between Oeltzen-Argelander North catalog and PPM...\n");

    int countDist = 0;
    double akkuDistError = 0.0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, true, EPOCH_OA);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo OARN */
    FILE *stream = fopen("cat/oarn.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read oarn.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		/* lee numeración */
		readField(buffer, cell, 1, 5);
		int oeltzenRef = atoi(cell);

		/* lee ascension recta B1842.0 */
		readFieldSanitized(buffer, cell, 11, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 13, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 15, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1842.0 (siempre positiva) */
		readFieldSanitized(buffer, cell, 20, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 22, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 24, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);
		findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
		if (minDistance < MAX_DIST_CROSS) {
            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;
		}

        if (!ppmFound) {
            printf("%d) Warning: OA %d is ALONE (nearest PPM star at %.1f arcsec).\n",
                ++errors,
                oeltzenRef,
                minDistance);
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

    printf("Available OA stars = %d\n", countOA);
    printf("Stars from OA identified with PPM = %d\n", countDist);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("Errors logged = %d\n", errors);
}

/*
 * readBAC - lee y cruza catalogo de la British Association
 */
void readBAC() {
    char buffer[1024], cell[256];

    printf("\n***************************************\n");
    printf("Perform comparison between BAC catalog and PPM...\n");

    int countDist = 0;
    double akkuDistError = 0.0;
    int errors = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, EPOCH_BAC);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo BAC */
    FILE *stream = fopen("cat/bac.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read bac.txt");
		exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
		/* lee numeración */
		readField(buffer, cell, 1, 4);
		int numRef = atoi(cell);

		/* lee ascension recta B1850.0 */
		readFieldSanitized(buffer, cell, 20, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 22, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 24, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1850.0 (en realidad NPD) */
		readFieldSanitized(buffer, cell, 48, 3);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 51, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 53, 3);
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
		if (minDistance < MAX_DIST_CROSS) {
            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;
		}

        if (!ppmFound) {
            printf("%d) Warning: BAC %d is ALONE (nearest PPM star at %.1f arcsec).\n",
                ++errors,
                numRef,
                minDistance);
        }

        /* la almacenamos para futuras identificaciones */
        if (countBAC >= MAXBACSTAR) {
            printf("Error: too many BAC stars.\n");
            exit(1);
        }
        bacX[countBAC] = x;
        bacY[countBAC] = y;
        bacZ[countBAC] = z;
        bacRef[countBAC] = numRef;
        countBAC++;
    }
	fclose(stream);

    printf("Available BAC stars = %d\n", countBAC);
    printf("Stars from BAC identified with PPM = %d\n", countDist);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("Errors logged = %d\n", errors);
}

/*
 * readUSNO - lee catalogo de Yarnall-Frisby
 * tambien revisa referencias cruzadas a OARN, BAC y BD
 */
void readUSNO() {
    char buffer[1024], cell[256], catName[20];
    char catLine[64];

    /* usamos catalogo BD */
    struct DMstar_struct *BDstar = getDMStruct();

    printf("\n***************************************\n");
    printf("Check references between USNO and BD/BAC/OA...\n");

    int checkDM = 0;
    int checkBAC = 0;
    int checkWB = 0;
    int checkOA = 0;
    int errors = 0;

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

        /* si esta disponible, tambien lee precesiones y chequea (usamos constantes de Struve) */
        readFieldSanitized(buffer, cell, 54, 6);
        double preRA = atof(cell) / 1000.0;
        if (fabs(preRA) > EPS) {
            double realPreRA = 3.07196 + 1.33704 * dsin(RA) * dtan(Decl);
            double diff = fabs(preRA - realPreRA);
            if (diff > 0.0099) {
                printf("%d) Warning: U %d reports %.3f on precession RA but it should be %.3f (diff=%.3f).\n",
                    ++errors,
                    numRef,
                    preRA,
                    realPreRA,
                    diff);
                printf("     Register U %d: %s\n", numRef, catLine);    
            }
        }

        readFieldSanitized(buffer, cell, 74, 5);
        double preDecl = atof(cell) / 100.0;
        if (fabs(preDecl) > EPS) {
            double realPreDecl = 20.05554 * dcos(RA);
            double diff = fabs(preDecl - realPreDecl);
            if (diff > 0.099) {
                printf("%d) Warning: U %d reports %.2f on precession DECL but it should be %.2f (diff=%.2f).\n",
                    ++errors,
                    numRef,
                    preDecl,
                    realPreDecl,
                    diff);
                printf("     Register U %d: %s\n", numRef, catLine);    
            }
        }

	    /* convierte coordenadas a la de OARN y calcula rectangulares */
        double newRA = RA;
        double newDecl = Decl;
		double x, y, z;
        transform(EPOCH_USNO, EPOCH_OA, &newRA, &newDecl);
        sph2rec(newRA, newDecl, &x, &y, &z);

        /* lee referencia a catalogo OARN */
		readField(buffer, cell, 19, 4);
        if (!strncmp(cell, "OARN", 4)) {
    		readField(buffer, cell, 30, 5);
            int numRefCat = atoi(cell);
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

	    /* convierte coordenadas a la de WB y calcula rectangulares */
        newRA = RA;
        newDecl = Decl;
        transform(EPOCH_USNO, EPOCH_WB, &newRA, &newDecl);
        sph2rec(newRA, newDecl, &x, &y, &z);

        /* lee referencia a catalogo WB */
		readField(buffer, cell, 19, 4);
        if (!strncmp(cell, "WEI ", 4)) {
    		readField(buffer, cell, 28, 2);
            int RARef = atoi(cell);
    		readField(buffer, cell, 30, 5);
            int numRefCat = atoi(cell);
            for (int i = 0; i < countWB; i++) {
                if (wbRARef[i] != RARef) continue;
                if (wbNumRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, wbX[i], wbY[i], wbZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: U %d is FAR from WB %dh %d (dist = %.1f arcsec).\n",
                        ++errors,
                        numRef,
                        RARef,
                        numRefCat,
                        dist);
                    printf("     Register U %d: %s\n", numRef, catLine);    
                } else checkWB++;
            }
        }

	    /* convierte coordenadas a la de BAC y calcula rectangulares */
        newRA = RA;
        newDecl = Decl;
        transform(EPOCH_USNO, EPOCH_BAC, &newRA, &newDecl);
        sph2rec(newRA, newDecl, &x, &y, &z);

        /* lee referencia a catalogo BAC */
		readField(buffer, cell, 19, 4);
        if (!strncmp(cell, "BAC ", 4)) {
    		readField(buffer, cell, 30, 5);
            int numRefCat = atoi(cell);
            for (int i = 0; i < countBAC; i++) {
                if (bacRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, bacX[i], bacY[i], bacZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: U %d is FAR from BAC %d (dist = %.1f arcsec).\n",
                        ++errors,
                        numRef,
                        numRefCat,
                        dist);
                    printf("     Register U %d: %s\n", numRef, catLine);    
                } else checkBAC++;
            }
        }

	    /* convierte coordenadas a la de BD y calcula rectangulares */
        newRA = RA;
        newDecl = Decl;
        transform(EPOCH_USNO, 1855.0, &newRA, &newDecl);
        sph2rec(newRA, newDecl, &x, &y, &z);

        /* lee referencia a catalogo Durchmusterung, también puede ser
           del catálogo B.VI. cuando se dan grados en vez de horas
           (se identifica si tiene un signo + o -). */
		readField(buffer, cell, 19, 9);
        if (!strncmp(cell, "DM  ", 4) ||
                (!strncmp(cell, "BON6", 4) && 
                    (cell[8] == '+' || cell[8] == '-'))) {
		    readField(buffer, cell, 27, 1);
            bool signRef = (cell[0] == '-');
            readField(buffer, cell, 28, 2);
            int declRef = atoi(cell);
            readField(buffer, cell, 30, 5);
            int numRefCat = atoi(cell);

            int index = getDMindex(signRef, declRef, numRefCat);
            if (index == -1) {
                printf("DM not found for declRef = %d, numRef = %d.\n",
                    declRef,
                    numRefCat);
                bye("Cannot find DM star!");
            }
            double dist = 3600.0 * calcAngularDistance(x, y, z, BDstar[index].x, BDstar[index].y, BDstar[index].z);
            if (dist > MAX_DIST_CROSS) {
                printf("%d) Warning: U %d is FAR from BD star (dist = %.1f arcsec).\n",
                    ++errors,
                    numRef,
                    dist);
                printf("     Register U %d: %s\n", numRef, catLine);
                writeRegister(index, false);
            } else checkDM++;
        }
    }
	fclose(stream);
    printf("USNO properly identified with BD = %d\n", checkDM);
    printf("USNO properly identified with WB = %d\n", checkWB);
    printf("USNO properly identified with BAC = %d\n", checkBAC);
    printf("USNO properly identified with OA = %d\n", checkOA);
    printf("Errors logged = %d\n", errors);
}

/*
 * readUA - lee Uranometria Argentina
 * tambien revisa referencias cruzadas a BD
 */
void readUA() {
    char buffer[1024], cell[256], catName[20];
    char catLine[64];

    /* usamos catalogo BD */
    struct DMstar_struct *BDstar = getDMStruct();

    printf("\n***************************************\n");
    printf("Check references between UA and BD...\n");

    int checkDM = 0;
    int checkWB = 0;
    int errors = 0;

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

	    /* convierte coordenadas a la de BD y calcula rectangulares */
        double newRA = RA;
        double newDecl = Decl;
        double x, y, z;
        transform(EPOCH_UA, 1855.0, &newRA, &newDecl);
        sph2rec(newRA, newDecl, &x, &y, &z);

        /* lee referencia a catalogo Durchmusterung */
		readField(buffer, cell, 82, 3);
        if (!strncmp(cell, "DM.", 3)) {
            bool signRef = newDecl < 0.0;
            int declRef = (int) fabs(newDecl);
            readField(buffer, cell, 85, 4);
            char numStr[5];
            sscanf(cell, "%4[^, ]", numStr);
            int numRefCat = atoi(numStr);
            
            int index = getDMindex(signRef, declRef, numRefCat);
            if (index == -1) {
                bye("Cannot find DM star!");
            }
            double dist = 3600.0 * calcAngularDistance(x, y, z, BDstar[index].x, BDstar[index].y, BDstar[index].z);
            if (dist > MAX_DIST_CROSS) {
                printf("%d) Warning: %dG %s is FAR from BD star (dist = %.1f arcsec).\n",
                    ++errors,
                    gouldRef,
                    cstRef,
                    dist);
                printf("     Register %dG %s: %s\n", gouldRef, cstRef, catLine);
                writeRegister(index, false);
            } else checkDM++;
        }

	    /* convierte coordenadas a la de WB y calcula rectangulares */
        newRA = RA;
        newDecl = Decl;
        transform(EPOCH_UA, EPOCH_WB, &newRA, &newDecl);
        sph2rec(newRA, newDecl, &x, &y, &z);

        /* lee referencia a catalogo Weisse */
		readField(buffer, cell, 82, 4);
        if (!strncmp(cell, "WB.", 3) && cell[3] != '(') {
            int RARef = (int) floor(newRA / 15.0);
            readField(buffer, cell, 85, 4);
            char numStr[5];
            sscanf(cell, "%4[^, ]", numStr);
            int numRefCat = atoi(numStr);
            for (int i = 0; i < countWB; i++) {
                if (wbRARef[i] != RARef) continue;
                if (wbNumRef[i] != numRefCat) continue;
                double dist = 3600.0 * calcAngularDistance(x, y, z, wbX[i], wbY[i], wbZ[i]);
                if (dist > MAX_DIST_CROSS) {
                    printf("%d) Warning: %dG %s is FAR from WB %dh %d star (dist = %.1f arcsec).\n",
                        ++errors,
                        gouldRef,
                        cstRef,
                        RARef,
                        numRefCat,
                        dist);
                    printf("     Register %dG %s: %s\n", gouldRef, cstRef, catLine);
                } else checkWB++;
            }
        }
    }
	fclose(stream);
    printf("UA stars properly identified with BD = %d\n", checkDM);
    printf("UA stars properly identified with WB = %d\n", checkWB);
    printf("Errors logged = %d\n", errors);
}

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("CROSS_NORTH - Compare several catalogs.\n");
    printf("Made in 2025 by Daniel Severin.\n");

    /* leemos catalogo CD */
    readDM(CURATED ? "cat/bd_curated.txt" : "cat/bd.txt");

    /* leemos y cruzamos WB */
    readWB();

    /* leemos y cruzamos OA */
    readOARN();

    /* leemos y cruzamos BAC */
    readBAC();

    /* leemos y revisamos identificaciones de Yarnall */
    readUSNO();

    /* leemos, cruzamos y revisamos Uranometria Argentina */
    readUA();
    return 0;
}
