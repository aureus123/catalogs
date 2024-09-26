
/*
 * COMPARE_GC - Compara registros de CD contra GC (y otros catálogos)
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "trig.h"
#include "misc.h"

#define MAXGCSTAR 30000
#define MAX_DISTANCE 360.0   // 6 minutos de arco
#define MAX_DIST_CATALOGS 180.0  // 3 minutos de arco
#define MAX_MAGNITUDE 1.0
#define HUGE_NUMBER 9999999999

/* Para uso de la libreria WCS: */
#define WCS_B1950 2 /* B1950(FK4) right ascension and declination */
extern "C" void wcsconp(int sys1, int sys2, double eq1, double eq2, double ep1, double ep2,
             double *dtheta, double *dphi, double *ptheta, double *pphi);


struct GCstar_struct {
	bool discard; /* true si debe ser descartada (por doble) */
    int gcRef; /* identificador con numero */
    int obs; /* cantidad de observaciones */
	int RAh, RAm, RAs, Decld, Declm, Decls; /* detalle primera observacion */
	double x, y, z; /* primera obs. en coordenadas rectangulares */
    double RA1875[10], Decl1875[10], epoch[10]; /* coordenadas y epoca de la observacion */
    double vmag; /* magnitud visual */
    int page; /* pagina donde se encuentra */
	/* CD Cross-index */
    int cdIndex; /* indice (no identificador) a la estrella CD más cercana */
	double dist; /* distancia a la estrella CD */
	int cdIndexWithinMag; /* indice a la CD más cercana, pero dentro del rango de magnitud */
	double distWithinMag;
	/* Yarnall Cross-index */
	int yarnallRef; /* identificador a catalogo de la USNO */
	char yarnallCat[29]; /* referencia a otros catalogos */
	double distYarnall; /* distancia a USNO */
	double vmagYarnall; /* magnitud USNO */
	/* Weiss Cross-index */
	int weissRef, oeltzenRef; /* identificador a catalogo de Weiss y OA */
	double distWeiss; /* distancia a Weiss */
	double distOeltzen; /* distancia a OA */
	double vmagWeiss; /* magnitud Weiss */
	/* Stone Cross-index */
	int stoneRef, lacailleRef; /* identificador a catalogo de Stone y Lacaille */
	double distStone; /* distancia a Stone */
	double vmagStone; /* magnitud Stone */
	/* Gillis Cross-index */
	int giRef; /* identificador a catalogo de Gillis */
	char giCat[19]; /* referencia a otros catalogos */
	double distGi; /* distance a Gillis */
	double vmagGi; /* magnitud Gillis */
} GCstar[MAXGCSTAR];

int GCstars;

void writeRegisterGC(int index) {
	printf("     Register GC %d: mag = %.1f, RA = %02dh%02dm%02ds%02d, DE = %02d°%02d'%02d''%01d, obs = %d (pag %d)\n",
		GCstar[index].gcRef,
		GCstar[index].vmag,
		GCstar[index].RAh,
		GCstar[index].RAm,
		GCstar[index].RAs / 100,
		GCstar[index].RAs % 100,
		GCstar[index].Decld,
		GCstar[index].Declm,
		GCstar[index].Decls / 10,
		GCstar[index].Decls % 10,
		GCstar[index].obs,
		GCstar[index].page
	);
	if (GCstar[index].yarnallRef != -1) {
		printf("       corresponds to USNO %d <%s> (mag=%.1f) at %.1f arcsec.\n",
			GCstar[index].yarnallRef,
			GCstar[index].yarnallCat,
			GCstar[index].vmagYarnall,
			GCstar[index].distYarnall
		);
	}
	if (GCstar[index].weissRef != -1) {
		printf("       corresponds to WEI %d or OA %d (mag=%.1f) at %.1f arcsec.\n",
			GCstar[index].weissRef,
			GCstar[index].oeltzenRef,
			GCstar[index].vmagWeiss,
			GCstar[index].distWeiss
		);
	}
	if (GCstar[index].stoneRef != -1) {
		char lacaille[20];
		lacaille[0] = 0;
		if (GCstar[index].lacailleRef != -1) {
			snprintf(lacaille, 20, "or L. %d ", GCstar[index].lacailleRef);
		}
		printf("       corresponds to ST %d %s(mag=%.1f) at %.1f arcsec.\n",
			GCstar[index].stoneRef,
			lacaille,
			GCstar[index].vmagStone,
			GCstar[index].distStone
		);
	}
	if (GCstar[index].giRef != -1) {
		printf("       corresponds to GI %d <%s> (mag=%.1f) at %.1f arcsec.\n",
			GCstar[index].giRef,
			GCstar[index].giCat,
			GCstar[index].vmagGi,
			GCstar[index].distGi
		);
	}
}

void copy(char *dest, char *src) {
	bool previousIsSpace = true;
	int srcPtr = 0;
	int destPtr = 0;
	while (src[srcPtr] != 0) {
		if (src[srcPtr] != ' ') {
			dest[destPtr++] = src[srcPtr];
			previousIsSpace = false;
		} else {
			if (!previousIsSpace) {
				// omit further consecutive spaces
				previousIsSpace = true;
				dest[destPtr++] = src[srcPtr];
			}
		}
		srcPtr++;
	}
	dest[destPtr] = 0;
}

/* lee estrellas del Primer Catalogo Argentino, coordenadas 1875.0
 * supuestamente todas estas estrellas deberian estar incluidas en el CD (excepto las que estan fuera de la faja)
 */
void readGC()
{
    FILE *stream;
    char buffer[1024], cell[256];
	double vmag;

    stream = fopen("cat/gc.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read gc.txt");
		exit(1);
    }

    struct DMstar_struct *CDstar = getDMStruct();
	int CDstars = getDMStars();

    int page = 1;
    int entry = 0;
    GCstars = 0;
	vmag = 0.0;
    while (fgets(buffer, 1023, stream) != NULL) {
		//if (GCstars >= 500) break;
		entry++;
		if ((entry-53) % 70 == 0) {
			page++;
			printf("Progress: stars = %d, page = %d\n", GCstars, page);
		}

		/* lee numeracion */
		readField(buffer, cell, 1, 5);
        int gcRef = atoi(cell);

		/* ver si es cumulo, nebulosa o variable */
		char type = buffer[11-1];
		//if (type == 'C' || type == 'N') continue;
		if (type == 'V') {
	    	vmag = 0.0;
		}
		else {
	    	/* lee magnitud (excepto si son espacios, en cuyo caso la magnitud y variabilidad es de la entrada anterior) */
	    	readField(buffer, cell, 8, 3);
	    	if (cell[0] != ' ') {
				if (cell[2] == ' ') cell[2] = '0';
				vmag = atof(cell)/10.0;
	    	}
		}

		/* lee epoca en que fue hecha la observacion */
		readField(buffer, cell, 12, 4);
        double epoch = (atof(cell)/100.0) + 1800.0;

		/* lee ascension recta B1875.0 */
		readField(buffer, cell, 16, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readField(buffer, cell, 18, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readField(buffer, cell, 20, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1875.0 */
		readField(buffer, cell, 39, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readField(buffer, cell, 41, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readField(buffer, cell, 43, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */
        if (Decl > -22) continue;

        /* calcula coordenadas rectangulares) */
        double x, y, z;
        sph2rec(RA, Decl, &x, &y, &z);

        /* Busca la estrella asociada en DM; en caso de haber más de
		 * una, escoge la de menor distancia */
        int cdIndex = -1;
        double minDistance = HUGE_NUMBER;
		int cdIndexWithinMag = -1;
		double minDistWithinMag = minDistance;
		for (int i = 0; i < CDstars; i++) {
        	double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
        	if (minDistance > dist) {
        		cdIndex = i;
        		minDistance = dist;
        	}
        	double cdVmag = CDstar[i].vmag;
			if (fabs(vmag - cdVmag) <= MAX_MAGNITUDE) {
        		if (minDistWithinMag > dist) {
            		cdIndexWithinMag = i;
            		minDistWithinMag = dist;
          		}
 		  	}
		}

		int previous = GCstars - 1;
		if (previous >= 0 && gcRef == GCstar[previous].gcRef) {
	    	int obs = GCstar[previous].obs;
			if (obs == 10) {
				printf("Max amount (of obs.) reached!\n");
				exit(1);
	    	}
	
	    	/* almacena otra entrada de la misma estrella */
	    	GCstar[previous].RA1875[obs] = RA;
	    	GCstar[previous].Decl1875[obs] = Decl;
	    	GCstar[previous].epoch[obs] = epoch;
	    	GCstar[previous].obs++;

			/* actualiza estrella CD si es más cercana */
			if (minDistance < GCstar[previous].dist) {
				GCstar[previous].cdIndex = cdIndex;
				GCstar[previous].dist = minDistance;
			}
			if (minDistWithinMag < GCstar[previous].distWithinMag) {
				GCstar[previous].cdIndexWithinMag = cdIndexWithinMag;
				GCstar[previous].distWithinMag = minDistWithinMag;
			}
			continue;
		}

	    if (GCstars == MAXGCSTAR) {
			printf("Max amount reached!\n");
			exit(1);
	    }

	    /* almacena la estrella por primera vez */
		GCstar[GCstars].discard = false;
	    GCstar[GCstars].gcRef = gcRef;
	    GCstar[GCstars].obs = 1;
		GCstar[GCstars].RAh = RAh;
		GCstar[GCstars].RAm = RAm;
		GCstar[GCstars].RAs = RAs;
		GCstar[GCstars].Decld = Decld;
		GCstar[GCstars].Declm = Declm;
		GCstar[GCstars].Decls = Decls;
		GCstar[GCstars].x = x;
		GCstar[GCstars].y = y;
		GCstar[GCstars].z = z;
	    GCstar[GCstars].RA1875[0] = RA;
	    GCstar[GCstars].Decl1875[0] = Decl;
	    GCstar[GCstars].epoch[0] = epoch;
	    GCstar[GCstars].vmag = vmag;
	    GCstar[GCstars].page = page;
	    GCstar[GCstars].cdIndex = cdIndex;
	    GCstar[GCstars].dist = minDistance;
		GCstar[GCstars].cdIndexWithinMag = cdIndexWithinMag;
		GCstar[GCstars].distWithinMag = minDistWithinMag;
		GCstar[GCstars].yarnallRef = -1;
		GCstar[GCstars].distYarnall = HUGE_NUMBER;
		GCstar[GCstars].weissRef = -1;
		GCstar[GCstars].oeltzenRef = -1;
		GCstar[GCstars].distWeiss = HUGE_NUMBER;
		GCstar[GCstars].stoneRef = -1;
		GCstar[GCstars].lacailleRef = -1;
		GCstar[GCstars].distStone = HUGE_NUMBER;
		GCstar[GCstars].giRef = -1;
		GCstar[GCstars].distGi = HUGE_NUMBER;

		/* proxima estrella */
		GCstars++;
		//printf("Pos %d: id=%d RA=%.4f Decl=%.4f (%.2f) Vmag=%.1f\n", GCstars, gcRef, RA, Decl, epoch, vmag);
    }
    fclose(stream);

	/* descartamos aquellas dobles cercanas débiles (cuya estrella CD es la misma) */
	int discarded = 0;
	for (int i = 0; i < GCstars; i++) {
		if (GCstar[i].discard) continue;
		int cdIndex1 = GCstar[i].cdIndex;
		for (int j = -5; j <= +5; j++) {
			if (j != 0) {
				if (GCstar[i+j].discard) continue;
				int cdIndex2 = GCstar[i+j].cdIndex;
				if (cdIndex1 == cdIndex2) {
					if (GCstar[i].vmag < GCstar[i+j].vmag) {
						//printf("Discarding GC %d\n", GCstar[i+j].gcRef);
						GCstar[i+j].discard = true;
					} else {
						//printf("Discarding GC %d\n", GCstar[i].gcRef);
						GCstar[i].discard = true;
					}
					discarded++;
					break;
				}
			}
		}
	}
		
	/* añadimos referencia al catálogo de la USNO */
	printf("Reading USNO catalog...\n");
    stream = fopen("cat/yarnall.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read yarnall.txt");
		exit(1);
    }
	int countUSNO = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
		/* lee signo declinación y descarta hemisferio norte */
		readField(buffer, cell, 60, 1);
		if (cell[0] != '-') continue;

		/* lee ascension recta B1860.0 */
		readField(buffer, cell, 40, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readField(buffer, cell, 42, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readField(buffer, cell, 44, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1860.0 */
		readField(buffer, cell, 61, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readField(buffer, cell, 63, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readField(buffer, cell, 65, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */
    
	    /* convierte coordenadas a 1875.0 */
        double RA1875 = RA;
        double Decl1875 = Decl;
        double pmRA = 0.0;
        double pmDecl = 0.0;
        wcsconp(WCS_B1950, WCS_B1950, 1860.0, 1875.0, 1860.0, 1875.0, &RA1875, &Decl1875, &pmRA, &pmDecl);
    	if (Decl1875 > -22) continue;

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RA1875, Decl1875, &x, &y, &z);

		/* lee magnitud y referencia */
		readField(buffer, cell, 37, 2);
		vmag = atof(cell)/10.0;
		readField(buffer, cell, 1, 5);
		int yarnallRef = atoi(cell);
		readField(buffer, cell, 7, 28);
		cell[28] = 0;

		/* barre estrellas de GC para identificarlas con este catálogo */
        double minDistance = HUGE_NUMBER;
		int gcIndex = -1;
		for (int i = 0; i < GCstars; i++) {
			if (GCstar[i].discard) continue;
			double dist = 3600.0 * calcAngularDistance(x, y, z, GCstar[i].x, GCstar[i].y, GCstar[i].z);
			if (dist < minDistance) {
				gcIndex = i;
				minDistance = dist;
			}
		}
		if (minDistance > MAX_DIST_CATALOGS) {
			//printf("Warning: nearest star is GC %d at %.1f arcsec. from USNO %d\n",
			//	GCstar[gcIndex].gcRef,
			//	minDistance,
			//	yarnallRef);
			continue;
		}
		if (minDistance < GCstar[gcIndex].distYarnall) {
			copy(GCstar[gcIndex].yarnallCat, cell);
			GCstar[gcIndex].yarnallRef = yarnallRef;
			GCstar[gcIndex].distYarnall = minDistance;
			GCstar[gcIndex].vmagYarnall = vmag;
			countUSNO++;
		}
	}
    fclose(stream);
	printf("  cross-referenced with %d stars\n", countUSNO);
		
	/* añadimos referencia al catálogo de Weiss */
	printf("Reading Weiss catalog...\n");
    stream = fopen("cat/weiss.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read weiss.txt");
		exit(1);
    }
	int countWeiss = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
		/* descarta otras obs que no sean la primera */
		readField(buffer, cell, 6, 1);
		if (cell[0] != '1') continue;

		/* lee ascension recta B1850.0 */
		readField(buffer, cell, 13, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readField(buffer, cell, 15, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readField(buffer, cell, 17, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1850.0 */
		readField(buffer, cell, 22, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readField(buffer, cell, 24, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readField(buffer, cell, 26, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */
    
	    /* convierte coordenadas a 1875.0 */
        double RA1875 = RA;
        double Decl1875 = Decl;
        double pmRA = 0.0;
        double pmDecl = 0.0;
        wcsconp(WCS_B1950, WCS_B1950, 1850.0, 1875.0, 1850.0, 1875.0, &RA1875, &Decl1875, &pmRA, &pmDecl);
    	if (Decl1875 > -22) continue;

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RA1875, Decl1875, &x, &y, &z);

		/* lee numeración catálogo Weiss y Oeltzen-Argelander */
		readField(buffer, cell, 1, 5);
		int weissRef = atoi(cell);
		readField(buffer, cell, 49, 5);
		int oeltzenRef = atoi(cell);

		/* lee magnitud */
		readField(buffer, cell, 9, 1);
		vmag = atof(cell);
		readField(buffer, cell, 11, 1);
		vmag += atof(cell)/10.0;

		/* barre estrellas de GC para identificarlas con este catálogo */
        double minDistance = HUGE_NUMBER;
		int gcIndex = -1;
		for (int i = 0; i < GCstars; i++) {
			if (GCstar[i].discard) continue;
			double dist = 3600.0 * calcAngularDistance(x, y, z, GCstar[i].x, GCstar[i].y, GCstar[i].z);
			if (dist < minDistance) {
				gcIndex = i;
				minDistance = dist;
			}
		}
		if (minDistance > MAX_DIST_CATALOGS) {
			//printf("Warning: nearest star is GC %d at %.1f arcsec. from WEI %d or OA %d\n",
			//	GCstar[gcIndex].gcRef,
			//	minDistance,
			//	weissRef,
			//	oeltzenRef);
			continue;
		}
		if (minDistance < GCstar[gcIndex].distWeiss) {
			GCstar[gcIndex].weissRef = weissRef;
			GCstar[gcIndex].oeltzenRef = oeltzenRef;
			GCstar[gcIndex].distWeiss = minDistance;
			GCstar[gcIndex].vmagWeiss = vmag;
			countWeiss++;
		}
	}
    fclose(stream);
	printf("  cross-referenced with %d stars\n", countWeiss);

	/* añadimos referencia al catálogo de Stone */
	printf("Reading Stone catalog...\n");
    stream = fopen("cat/stone1.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read stone1.txt");
		exit(1);
    }
	FILE *stream2 = fopen("cat/stone2.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read stone2.txt");
		exit(1);
    }

	int countStone = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
		char buffer2[1024];
		fgets(buffer2, 1023, stream2);

		/* lee ascension recta B1880.0 */
		readField(buffer, cell, 33, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readField(buffer, cell, 35, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readField(buffer, cell, 37, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1880.0 (en realidad NPD) */
		readField(buffer2, cell, 21, 3);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readField(buffer2, cell, 24, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readField(buffer2, cell, 26, 4);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/100.0)/3600.0;
		Decl = 90.0 - Decl; /* convertimos NPD a declinación */

	    /* convierte coordenadas a 1875.0 según cada mean time */
        double RA1875 = RA;
        double Decl1875 = Decl;
        double pmRA = 0.0;
        double pmDecl = 0.0;
        wcsconp(WCS_B1950, WCS_B1950, 1880.0, 1875.0, 1880.0, 1875.0, &RA1875, &Decl1875, &pmRA, &pmDecl);
    	if (Decl1875 > -22) continue;

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RA1875, Decl1875, &x, &y, &z);

		/* lee numeración catálogo Stone y Lacaille */
		readField(buffer, cell, 8, 5);
		int stoneRef = atoi(cell);
		readField(buffer, cell, 14, 4);
		int lacailleRef = cell[0] == '&' ? -1 : atoi(cell);

		/* lee magnitud */
		readField(buffer, cell, 24, 1);
		vmag = atof(cell);
		readField(buffer, cell, 25, 1);
		vmag += cell[0] == '&' ? 0 : atof(cell)/10.0;

		/* barre estrellas de GC para identificarlas con este catálogo */
        double minDistance = HUGE_NUMBER;
		int gcIndex = -1;
		for (int i = 0; i < GCstars; i++) {
			if (GCstar[i].discard) continue;
			double dist = 3600.0 * calcAngularDistance(x, y, z, GCstar[i].x, GCstar[i].y, GCstar[i].z);
			if (dist < minDistance) {
				gcIndex = i;
				minDistance = dist;
			}
		}
		if (minDistance > MAX_DIST_CATALOGS) {
			//printf("Warning: nearest star is GC %d at %.1f arcsec. from ST %d\n",
			//	GCstar[gcIndex].gcRef,
			//	minDistance,
			//	stoneRef);
			continue;
		}
		if (minDistance < GCstar[gcIndex].distStone) {
			GCstar[gcIndex].stoneRef = stoneRef;
			GCstar[gcIndex].lacailleRef = lacailleRef;
			GCstar[gcIndex].distStone = minDistance;
			GCstar[gcIndex].vmagStone = vmag;
			countStone++;
		}
	}
    fclose(stream2);
    fclose(stream);
	printf("  cross-referenced with %d stars\n", countStone);

	/* añadimos referencia al catálogo de Gillis */
	printf("Reading Gillis catalog...\n");
    stream = fopen("cat/gillis.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read gillis.txt");
		exit(1);
    }
	int countGi = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
		/* lee signo declinación y descarta hemisferio norte */
		readField(buffer, cell, 48, 1);
		if (cell[0] != '-') continue;

		/* lee ascension recta B1850.0 */
		readField(buffer, cell, 33, 2);
		int RAh = atoi(cell);
        double RA = (double) RAh;
		readField(buffer, cell, 35, 2);
		int RAm = atoi(cell);
        RA += ((double) RAm)/60.0;
		readField(buffer, cell, 37, 4);
		int RAs = atoi(cell);
        RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1850.0 */
		readField(buffer, cell, 49, 2);
		int Decld = atoi(cell);
        double Decl = (double) Decld;
		readField(buffer, cell, 51, 2);
		int Declm = atoi(cell);
        Decl += ((double) Declm)/60.0;
		readField(buffer, cell, 53, 3);
		int Decls = atoi(cell);
        Decl += (((double) Decls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */
    
	    /* convierte coordenadas a 1875.0 */
        double RA1875 = RA;
        double Decl1875 = Decl;
        double pmRA = 0.0;
        double pmDecl = 0.0;
        wcsconp(WCS_B1950, WCS_B1950, 1850.0, 1875.0, 1850.0, 1875.0, &RA1875, &Decl1875, &pmRA, &pmDecl);
    	if (Decl1875 > -22) continue;

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RA1875, Decl1875, &x, &y, &z);

		/* lee magnitud y referencia */
		readField(buffer, cell, 29, 3);
		vmag = atof(cell)/10.0;
		readField(buffer, cell, 1, 5);
		int giRef = atoi(cell);
		readField(buffer, cell, 7, 18);
		cell[18] = 0;

		/* barre estrellas de GC para identificarlas con este catálogo */
        double minDistance = HUGE_NUMBER;
		int gcIndex = -1;
		for (int i = 0; i < GCstars; i++) {
			if (GCstar[i].discard) continue;
			double dist = 3600.0 * calcAngularDistance(x, y, z, GCstar[i].x, GCstar[i].y, GCstar[i].z);
			if (dist < minDistance) {
				gcIndex = i;
				minDistance = dist;
			}
		}
		if (minDistance > MAX_DIST_CATALOGS) {
			//printf("Warning: nearest star is GC %d at %.1f arcsec. from GI %d\n",
			//	GCstar[gcIndex].gcRef,
			//	minDistance,
			//	giRef);
			continue;
		}
		if (minDistance < GCstar[gcIndex].distGi) {
			copy(GCstar[gcIndex].giCat, cell);
			GCstar[gcIndex].giRef = giRef;
			GCstar[gcIndex].distGi = minDistance;
			GCstar[gcIndex].vmagGi = vmag;
			countGi++;
		}
	}
    fclose(stream);
	printf("  cross-referenced with %d stars\n", countGi);
		
	printf("Stars read from Catalogo General Argentino: %d  (discarded = %d)\n", GCstars, discarded);
}

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("COMPARE_GC - Compare CD and GC catalogs.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    /* leemos catalogo CD */
    readDM("cat/cd.txt");
    struct DMstar_struct *CDstar = getDMStruct();

    /* leemos catalogo GC */
    readGC();

    int maxDistError = 0;
    int magDiffError = 0;
    int indexError = 0;
    int goodStars = 0;
    double akkuDistError = 0.0;
    double akkuDeltaError = 0.0;
    for (int i = 0; i < GCstars; i++) {
		if (GCstar[i].discard) continue;
        int cdIndex = GCstar[i].cdIndex;
        float dist = GCstar[i].dist;
        if (dist > MAX_DISTANCE) {
            // posiciones muy separadas, supera umbral
            maxDistError++;
            indexError++;
            printf("%d) GC %d is ALONE; nearest CD %d°%d separated in %.1f arcsec.\n",
                indexError,
                GCstar[i].gcRef,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                dist);
			writeRegisterGC(i);	
            continue;
        }

        int cdIndexWithinMag = GCstar[i].cdIndexWithinMag;
		double gcVmag = GCstar[i].vmag;
        double cdVmag = CDstar[cdIndexWithinMag].vmag;
        if (gcVmag < 0.00001 || cdVmag > 29.9) continue; // se omiten aquellas estrellas con Vmag=0 o variables

        float distWithinMag = GCstar[i].distWithinMag;
        if (distWithinMag > MAX_DISTANCE) {
            // supera umbral (estrella CD mas cercana de magnitud similar)
            magDiffError++;
            indexError++;
            printf("%d) GC %d is alone within MAG = %.1f; similar CD %d°%d separated in %.1f arcsec although nearest CD %d°%d at %.1f arcsec.\n",
                indexError,
                GCstar[i].gcRef,
                gcVmag,
                CDstar[cdIndexWithinMag].declRef,
                CDstar[cdIndexWithinMag].numRef,
                distWithinMag,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                dist);
			writeRegisterGC(i);	
            writeRegister(cdIndexWithinMag);
            writeRegister(cdIndex);
            continue;
        }
        goodStars++;
        akkuDistError += distWithinMag * distWithinMag;
		double delta = fabs(gcVmag - cdVmag);
        akkuDeltaError += delta * delta;
    }
    printf("Total errors: %d (alone: %d, mag: %d)\n",
        indexError, maxDistError, magDiffError);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)goodStars),
        goodStars);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)goodStars),
        goodStars);
    return 0;
}
