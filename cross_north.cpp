
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

#define MAXOASTAR 26500
#define EPOCH_OA 1842.0
#define EPOCH_USNO 1860.0
#define EPOCH_UA 1875.0
#define MAX_DIST_OA_PPM 45.0
#define MAX_DIST_CROSS 60.0
#define CURATED true // true if curated BD catalog should be used

// Here, we save 1842.0 coordinates of OA stars in rectangular form
double oaX[MAXOASTAR], oaY[MAXOASTAR], oaZ[MAXOASTAR];
int oaRef[MAXOASTAR];
int countOA = 0;

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
		if (minDistance < MAX_DIST_OA_PPM) {
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
 * readUSNO - lee catalogo de Yarnall-Frisby
 * tambien revisa referencias cruzadas a OARN y BD
 */
void readUSNO() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];
    char catLine[64];

    /* usamos catalogo BD */
    struct DMstar_struct *BDstar = getDMStruct();

    printf("\n***************************************\n");
    printf("Check references between USNO and BD/OA...\n");

    int errors = 0;

    /* leemos catalogo USNO */
    FILE *stream = fopen("cat/yarnall.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read yarnall.txt");
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
                }
            }
        }

	    /* convierte coordenadas a la de BD y calcula rectangulares */
        newRA = RA;
        newDecl = Decl;
        transform(EPOCH_USNO, 1855.0, &newRA, &newDecl);
        sph2rec(newRA, newDecl, &x, &y, &z);

        /* lee referencia a catalogo Durchmusterung */
		readField(buffer, cell, 19, 4);
        if (!strncmp(cell, "DM  ", 4)) {
		    readField(buffer, cell, 27, 1);
            bool signRef = (cell[0] == '-');
            readField(buffer, cell, 28, 2);
            int declRef = atoi(cell);
            readField(buffer, cell, 30, 5);
            int numRefCat = atoi(cell);

            bool isFirstVolume = false;
            if (signRef) {
                if (declRef <= 1) {
                    /* primer volumen de BD llega hasta declinación -01 (inclusive) */
                    isFirstVolume = true;
                }
            } else {
                if (declRef <= 19) {
                    /* primer volumen de BD comienza en declinación +19 */
                    isFirstVolume = true;
                }
            }

            if (isFirstVolume) {
                int index = getDMindex(signRef, declRef, numRefCat);
                if (index == -1) {
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
                }
            }
        }
    }
	fclose(stream);
    printf("Errors logged = %d\n", errors);
}

/*
 * readUA - lee Uranometria Argentina
 * tambien revisa referencias cruzadas a BD
 */
void readUA() {
    char buffer[1024], cell[256], catName[20], cdName[20], ppmName[20];
    char catLine[64];

    /* usamos catalogo BD */
    struct DMstar_struct *BDstar = getDMStruct();

    printf("\n***************************************\n");
    printf("Check references between UA and BD...\n");

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
            int numRefCat = atoi(cell);

            bool isFirstVolume = false;
            if (signRef) {
                if (declRef <= 1) {
                    /* primer volumen de BD llega hasta declinación -01 (inclusive) */
                    isFirstVolume = true;
                }
            } else {
                if (declRef <= 19) {
                    /* primer volumen de BD comienza en declinación +19 */
                    isFirstVolume = true;
                }
            }
            if (!isFirstVolume) {
                bye("Found an inexistent DM star!");
            }

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
            }
        }
    }
	fclose(stream);
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

    /* leemos y cruzamos OA */
    readOARN();

    /* leemos y revisamos identificaciones de Yarnall */
    readUSNO();

    /* leemos, cruzamos y revisamos Uranometria Argentina */
    readUA();
    return 0;
}
