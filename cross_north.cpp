
/*
 * CROSS_NORTH - Compara registros de varios catalogos
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "read_ppm.h"
#include "read_dm.h"
#include "read_gc.h"
#include "trig.h"
#include "misc.h"
#include "find_gsc.h"
#include "cross_utils.h"

#define MAXWBSTAR 31900
#define MAXOASTAR 26500
#define MAXBACSTAR 8400
#define EPOCH_WB 1825.0
#define EPOCH_OA 1842.0
#define EPOCH_BAC 1850.0
#define EPOCH_USNO 1860.0
#define EPOCH_UA 1875.0
#define EPOCH_GC 1875.0
#define MAX_DIST_CAT_PPM 90.0
#define MAX_DIST_OA_PPM 45.0
#define CURATED true // true if curated BD catalog should be used

// Here, we save 1825.0 coordinates of WB stars in rectangular form
double wbX[MAXWBSTAR], wbY[MAXWBSTAR], wbZ[MAXWBSTAR];
int wbRARef[MAXWBSTAR], wbNumRef[MAXWBSTAR];
int countWB = 0;
StarList wbList = {&countWB, wbNumRef, wbX, wbY, wbZ};

// Here, we save 1842.0 coordinates of OA stars in rectangular form
double oaX[MAXOASTAR], oaY[MAXOASTAR], oaZ[MAXOASTAR];
int oaRef[MAXOASTAR];
int countOA = 0;
StarList oaList = {&countOA, oaRef, oaX, oaY, oaZ};

// Here, we save 1850.0 coordinates of BAC stars in rectangular form
double bacX[MAXBACSTAR], bacY[MAXBACSTAR], bacZ[MAXBACSTAR];
int bacRef[MAXBACSTAR];
int countBAC = 0;
StarList bacList = {&countBAC, bacRef, bacX, bacY, bacZ};

/*
 * readWB - lee y cruza catalogo de Weisse
 */
void readWB() {
    char buffer[1024], cell[256], catName[20];

    printf("\n***************************************\n");
    printf("Perform comparison between Weisse catalog and PPM...\n");

    CrossStats stats;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    preparePPM(EPOCH_WB, false);

    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("wb", &crossPPMStream, &crossSAOStream, &crossHDStream);

    /* leemos catalogo WB */
    FILE *stream = fopen("cat/wb.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read wb.txt");
		exit(1);
    }
    int lineNum = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
        lineNum++;
		/* lee numeración */
		readField(buffer, cell, 14, 2);
		int RARef = atoi(cell);
		readField(buffer, cell, 6, 4);
		int numRef = atoi(cell);

		/* lee magnitud */
		float vmag = 0.0;
		readField(buffer, cell, 10, 2);
		if (cell[0] != ' ') {
		    int magInt = atoi(cell);
		    vmag = (float) magInt;
		    readField(buffer, cell, 12, 2);
		    if (cell[0] != ' ') {
		        int magFrac = atoi(cell);
		        if (magFrac == magInt + 1) {
		            vmag += 0.3;
		        } else if (magFrac == magInt - 1) {
		            vmag -= 0.3;
		        } else {
		            printf("Error: line %d has unexpected value %d in magnitude columns 12-13 (expected %d, %d or blank). Using integer magnitude.\n",
		                lineNum, magFrac, magInt - 1, magInt + 1);
		        }
		    }
		}

		/* lee ascension recta B1825.0 (la hora es la misma numeración RARef) */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 14, &RAh, &RAm, &RAs);

		/* lee declinacion B1825.0 */
        int Decld, Declm, Decls;
        double Decl = readDeclField(buffer, 31, 2, 33, 35, 3, 10.0, &Decld, &Declm, &Decls);
        readField(buffer, cell, 30, 1);
		if (cell[0] == '-') Decl = -Decl;

		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z, minDistance;
		int ppmIndex;
		sph2rec(RA, Decl, &x, &y, &z);
        snprintf(catName, 20, "WB %dh %d", RARef, numRef);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_CAT_PPM, NULL, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);

        warnIfAloneNorth(ppmFound, minDistance, RA, Decl, EPOCH_WB, catName, &stats.errors);

        /* la almacenamos para futuras identificaciones */
        int i = storeStar(&countWB, MAXWBSTAR, "WB", wbNumRef, wbX, wbY, wbZ, NULL, numRef, x, y, z, 0.0);
        wbRARef[i] = RARef;
    }
	fclose(stream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);

    printf("Available WB stars = %d\n", countWB);
    printf("Stars from WB identified with PPM = %d\n", stats.countDist);
    printRSMEDist(&stats);
    printf("Stars not identified with PPM nor GSC = %d\n", stats.errors);
}

/*
 * readOARN - lee y cruza catalogo de Oeltzen-Argelander (North)
 */
void readOARN() {
    char buffer[1024], cell[256], catName[20];

    printf("\n***************************************\n");
    printf("Perform comparison between Oeltzen-Argelander North catalog and PPM...\n");

    CrossStats stats;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    preparePPM(EPOCH_OA, true);

    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("oarn", &crossPPMStream, &crossSAOStream, &crossHDStream);

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

		/* lee magnitud (columnas 6-9; en blanco si no hay magnitud) */
		readField(buffer, cell, 6, 2);
		int magInt = atoi(cell);
		float vmag = (float) magInt;
		readField(buffer, cell, 8, 2);
		if (cell[0] != ' ') {
		    int magFrac = atoi(cell);
		    if (magFrac == magInt + 1) {
		        vmag += 0.5;
		    } else if (magFrac == magInt - 1) {
		        vmag -= 0.2;
		    } else {
		        printf("Error: OA %d has unexpected value %d in magnitude columns 8-9 (expected %d, %d or blank).\n",
		            oeltzenRef, magFrac, magInt - 1, magInt + 1);
		        exit(1);
		    }
		}

		/* lee ascension recta B1842.0 */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 11, &RAh, &RAm, &RAs);

		/* lee declinacion B1842.0 (siempre positiva) */
        int Decld, Declm, Decls;
        double Decl = readDeclField(buffer, 20, 2, 22, 24, 3, 10.0, &Decld, &Declm, &Decls);

		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z, minDistance;
		int ppmIndex;
		sph2rec(RA, Decl, &x, &y, &z);
        snprintf(catName, 20, "OA %d", oeltzenRef);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_OA_PPM, NULL, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);

        warnIfAloneNorth(ppmFound, minDistance, RA, Decl, EPOCH_OA, catName, &stats.errors);

        /* la almacenamos para futuras identificaciones */
        storeStar(&countOA, MAXOASTAR, "OA", oaRef, oaX, oaY, oaZ, NULL, oeltzenRef, x, y, z, 0.0);
    }
	fclose(stream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);

    printf("Available OA stars = %d\n", countOA);
    printf("Stars from OA identified with PPM = %d\n", stats.countDist);
    printRSMEDist(&stats);
    printf("Stars not identified with PPM nor GSC = %d\n", stats.errors);
}

/*
 * readBAC - lee y cruza catalogo de la British Association
 */
void readBAC() {
    char buffer[1024], cell[256], catName[20];

    printf("\n***************************************\n");
    printf("Perform comparison between BAC catalog and PPM...\n");

    CrossStats stats;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    preparePPM(EPOCH_BAC, false);

    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("bac", &crossPPMStream, &crossSAOStream, &crossHDStream);

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

		/* lee magnitud */
		float vmag = readMagIntHalf(buffer, 17, 1, 18);

		/* lee ascension recta B1850.0 */
        int RAh, RAm, RAs;
        double RA = readRAField(buffer, 20, &RAh, &RAm, &RAs);

		/* lee declinacion B1850.0 (en realidad NPD) */
        int Decld, Declm, Decls;
        double Decl = readDeclField(buffer, 48, 3, 51, 53, 3, 10.0, &Decld, &Declm, &Decls);
		Decl = 90.0 - Decl; /* convertimos NPD a declinación */

		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z, minDistance;
		int ppmIndex;
		sph2rec(RA, Decl, &x, &y, &z);
        snprintf(catName, 20, "BAC %d", numRef);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_CAT_PPM, NULL, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);

        warnIfAloneNorth(ppmFound, minDistance, RA, Decl, EPOCH_BAC, catName, &stats.errors);

        /* la almacenamos para futuras identificaciones */
        storeStar(&countBAC, MAXBACSTAR, "BAC", bacRef, bacX, bacY, bacZ, NULL, numRef, x, y, z, 0.0);
    }
	fclose(stream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);

    printf("Available BAC stars = %d\n", countBAC);
    printf("Stars from BAC identified with PPM = %d\n", stats.countDist);
    printRSMEDist(&stats);
    printf("Stars not identified with PPM nor GSC = %d\n", stats.errors);
}

/*
 * readUSNO - lee catalogo de Yarnall-Frisby
 * tambien revisa referencias cruzadas a OARN, BAC y BD
 */
void readUSNO() {
    char buffer[1024], cell[256], srcName[20];
    char catLine[64];

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
        snprintf(srcName, 20, "U %d", numRef);

        /* si esta disponible, tambien lee precesiones y chequea (usamos constantes de Struve) */
        readFieldSanitized(buffer, cell, 54, 6);
        checkPrecessionRA(atof(cell) / 1000.0, srcName, catLine, RA, Decl, 3.07196, 1.33704, 0.0099, &errors);
        readFieldSanitized(buffer, cell, 74, 5);
        checkPrecessionDecl(atof(cell) / 100.0, srcName, catLine, RA, 20.05554, 0.099, &errors);

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
            checkCrossRef(srcName, catLine, "OA", x, y, z, atoi(cell), &oaList, false, -1, &checkOA, &errors);
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
            checkCrossRefWB(srcName, catLine, false, x, y, z, RARef, atoi(cell), &wbList, wbRARef, -1, &checkWB, &errors);
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
            checkCrossRef(srcName, catLine, "BAC", x, y, z, atoi(cell), &bacList, false, -1, &checkBAC, &errors);
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
            checkBDRef(srcName, catLine, true, signRef, declRef, atoi(cell), x, y, z, &checkDM, &errors);
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
 * readGCScanned - revisa referencias cruzadas de las páginas escaneadas de GC
 * (ya se deben haber leidos los catálogos WB)
 * Lee todos los archivos scans/rnao14_page*.csv y, para cada estrella (una fila),
 * calcula la distancia entre la estrella GC (1a columna; la 2a columna -magnitud-
 * se ignora) y la estrella de referencia (3a columna)
 */
void readGCScanned() {
    char buffer[1024], ref[64], catgName[20];

    printf("\n***************************************\n");
    printf("Perform cross-checking of scanned GC pages...\n");

    int errors = 0;
    int checkWB = 0;
    int countStars = 0, countRefs = 0;

    /* coordenadas (1875.0) de las estrellas GC ya leídas */
    struct GCstar_struct *GCstar = getGCStruct();

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
            if (!strncmp(ref, "WB.", 3)) {
                countRefs++;
                int numRefCat = atoi(&ref[3]);

        	    /* convierte coordenadas a la de WB y calcula rectangulares */
                double newRA = GCstar[gcIndex].RA1875;
                double newDecl = GCstar[gcIndex].Decl1875;
                transform(EPOCH_GC, EPOCH_WB, &newRA, &newDecl);
                sph2rec(newRA, newDecl, &x, &y, &z);

                int RARef = (int) floor(newRA / 15.0);
                checkCrossRefWB(catgName, NULL, true, x, y, z, RARef, numRefCat, &wbList, wbRARef, gcIndex, &checkWB, &errors);
            }
        }
        fclose(stream);
    }

    printf("Scanned GC rows with a reference = %d; references to checked catalogs (WB) = %d\n",
        countStars, countRefs);
    printf("GC properly identified with WB = %d\n", checkWB);
    printf("Errors logged = %d\n", errors);
}

/*
 * readUA - lee Uranometria Argentina
 * tambien revisa referencias cruzadas a BD
 */
void readUA() {
    char buffer[1024];

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
    UAstar_struct ua;
    while (fgets(buffer, 1023, stream) != NULL) {
        if (!parseUALine(buffer, &serpens, &ua)) continue;

	    /* convierte coordenadas a la de BD y calcula rectangulares */
        double newRA = ua.RA;
        double newDecl = ua.Decl;
        double x, y, z;
        transform(EPOCH_UA, 1855.0, &newRA, &newDecl);
        sph2rec(newRA, newDecl, &x, &y, &z);

        /* lee referencias a catalogos. */
        char subcell[3][18];
        int count = splitUARefs(buffer, subcell);

        for (int refs = 0; refs < count; refs++) {
            /* lee referencia a catalogo Durchmusterung */
            if (!strncmp(subcell[refs], "DM.", 3)) {
                int numRefCat = atoi(&subcell[refs][3]);
                bool signRef = newDecl < 0.0;
                int declRef = (int) fabs(newDecl);
                checkBDRef(ua.catgName, ua.catLine, false, signRef, declRef, numRefCat, x, y, z, &checkDM, &errors);
            }

            /* lee referencia a catalogo Weisse */
            if (!strncmp(subcell[refs], "WB.", 3) && subcell[refs][3] != '(') {
                int numRefCat = atoi(&subcell[refs][3]);

        	    /* convierte coordenadas a la de WB y calcula rectangulares */
                newRA = ua.RA;
                newDecl = ua.Decl;
                transform(EPOCH_UA, EPOCH_WB, &newRA, &newDecl);
                sph2rec(newRA, newDecl, &x, &y, &z);

                int RARef = (int) floor(newRA / 15.0);
                checkCrossRefWB(ua.catgName, ua.catLine, true, x, y, z, RARef, numRefCat, &wbList, wbRARef, -1, &checkWB, &errors);
            }
        }
    }
	fclose(stream);
    printf("UA stars properly identified with BD = %d\n", checkDM);
    printf("UA stars properly identified with WB = %d\n", checkWB);
    printf("Errors logged = %d\n", errors);
}

/*
 * readSOM - lee los "Standards of Magnitude" de la Uranometria Argentina
 * (cat/UA_standards.csv, epoca 1875.0) y los cruza con PPM
 * tambien revisa referencias cruzadas a BD (posicion y magnitud) y WB (posicion)
 */
void readSOM() {
    char buffer[1024], catName[20];
    char catLine[64];

    /* usamos catalogo BD */
    struct DMstar_struct *BDstar = getDMStruct();

    printf("\n***************************************\n");
    printf("Perform comparison between UA Standards of Magnitude and PPM...\n");

    CrossStats stats;
    int checkDM = 0;
    int checkWB = 0;

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    preparePPM(EPOCH_UA, false);

    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("som", &crossPPMStream, &crossSAOStream, &crossHDStream);

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

		/* lee magnitud */
        float vmag = (float) atof(field[7]);

		/* busca la PPM mas cercana y genera el cruzamiento */
		double x, y, z, minDistance;
		int ppmIndex;
		sph2rec(RA, Decl, &x, &y, &z);
        bool ppmFound = crossWithPPM(x, y, z, Decl, vmag, MAX_DIST_CAT_PPM, NULL, catName,
            crossPPMStream, crossSAOStream, crossHDStream, &ppmIndex, &minDistance, &stats);

        /* estrellas brillantes: no hace falta consultar GSC, con PPM alcanza */
        if (!ppmFound) {
            printf("%d) Warning: %s is ALONE (nearest PPM star at %.1f arcsec).\n",
                ++stats.errors,
                catName,
                minDistance);
        }

	    /* convierte coordenadas a la de BD y calcula rectangulares */
        double newRA = RA;
        double newDecl = Decl;
        transform(EPOCH_UA, 1855.0, &newRA, &newDecl);
        sph2rec(newRA, newDecl, &x, &y, &z);

        /* revisa la referencia a BD: posicion y magnitud */
        int numRefCat = atoi(field[8]);
        bool signRef = newDecl < 0.0;
        int declRef = (int) fabs(newDecl);
        int index = getDMindex(signRef, declRef, numRefCat);
        if (index == -1) {
            printf("%d) Warning: %s refers to BD %c%d°%d but it does not exist.\n",
                ++stats.errors,
                catName,
                signRef ? '-' : '+',
                declRef,
                numRefCat);
            printf("     Register %s: %s\n", catName, catLine);
        } else {
            double dist = 3600.0 * calcAngularDistance(x, y, z, BDstar[index].x, BDstar[index].y, BDstar[index].z);
            if (dist > MAX_DIST_CROSS) {
                printf("%d) Warning: %s is FAR from BD %c%d°%d star (dist = %.1f arcsec).\n",
                    ++stats.errors,
                    catName,
                    signRef ? '-' : '+',
                    declRef,
                    numRefCat,
                    dist);
                printf("     Register %s: %s\n", catName, catLine);
                writeRegister(index, false);
            } else checkDM++;

            double bdMag = atof(field[9]);
            if (fabs(bdMag - BDstar[index].vmag) > 0.1) {
                printf("%d) Warning: %s reports BD mag %.1f but BD %c%d°%d has mag %.1f.\n",
                    ++stats.errors,
                    catName,
                    bdMag,
                    signRef ? '-' : '+',
                    declRef,
                    numRefCat,
                    BDstar[index].vmag);
            }
        }

        /* revisa la referencia a WB, si existe: posicion
           (parentesis: identificación no confiable, no se revisa) */
        if (field[12][0] != 0 && field[12][0] != '(') {
            numRefCat = atoi(field[12]);

    	    /* convierte coordenadas a la de WB y calcula rectangulares */
            newRA = RA;
            newDecl = Decl;
            transform(EPOCH_UA, EPOCH_WB, &newRA, &newDecl);
            sph2rec(newRA, newDecl, &x, &y, &z);

            /* la serie de Weisse 1846 (cat/wb.txt) solo cubre hasta +15°
               en su propia epoca: mas alla la columna W.Bessel refiere a
               la serie de Weisse 1863, cuya numeracion no esta aqui */
            if (newDecl < 15.0) {
                int RARef = (int) floor(newRA / 15.0);
                checkCrossRefWB(catName, catLine, true, x, y, z, RARef, numRefCat, &wbList, wbRARef, -1, &checkWB, &stats.errors);
            }
        }
    }
	fclose(stream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);

    printf("Stars from SOM identified with PPM = %d\n", stats.countDist);
    printRSMEDist(&stats);
    printf("SOM stars properly identified with BD = %d\n", checkDM);
    printf("SOM stars properly identified with WB = %d\n", checkWB);
    printf("Errors logged = %d\n", stats.errors);
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

    /* leemos, cruzamos y revisamos Standards of Magnitude de la UA */
    readSOM();

    /* leemos y cruzamos OA */
    readOARN();

    /* leemos y cruzamos BAC */
    readBAC();

    /* leemos y revisamos identificaciones de Yarnall */
    readUSNO();

    /* leemos, cruzamos y revisamos Uranometria Argentina */
    readUA();

    /* leemos, cruzamos y revisamos GC */
    readGC();
    readGCScanned();
    return 0;
}
