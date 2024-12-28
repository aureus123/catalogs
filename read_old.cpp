
/*
 * READ_OLD - Lee varios catálogos antiguos
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
#include "read_old.h"

#define CD_SEARCH true // true if search nearest CD star over other old catalogs
#define PPM_SEARCH true // true if search nearest PPM star over other old catalogs (very intensive)
struct GCstar_struct GCstar[MAXGCSTAR];

int GCstars;

/*
 * getGCStars - devuelve la cantidad de estrellas de GC leidas
 */
int getGCStars()
{
    return GCstars;
}

/*
 * getGCStruct - devuelve la estructura GC
 */
struct GCstar_struct *getGCStruct()
{
    return &GCstar[0];
}

/*
 * writeRegister - escribe en pantalla un registro de GC
 */
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

/*
 * Tiene dos modos:
 * mode = true:
 *   Lee estrellas del Primer Catalogo Argentino, coordenadas 1875.0
 *   supuestamente todas estas estrellas deberian estar incluidas en el catálogo CD
 *   (excepto las que están fuera de la faja, y algunas de CD marcadas como "dobles")
 *   También lee otros catálogos antiguos (Yarnall, Stone, Weiss, Gilliss).
 * mode = false:
 *   Recibe una coordenada en 1875.0 y genera una estrella ficticia llamada GC 1.
 *   Luego la cruza con los catálogos CD y otros antiguos.
 */
void readGC(bool mode, int fictRAh, int fictRAm, int fictRAs, int fictDecld, int fictDeclm, int fictDecls)
{
    FILE *stream;
    char buffer[1024], cell[256];
	double vmag;
	
    struct DMstar_struct *CDstar = getDMStruct();
	int CDstars = getDMStars();
    int page = 1;
    int entry = 0;
	int discarded = 0;
    GCstars = 0;

	if (PPM_SEARCH) {
		/* leemos catalogo PPM (de ser necesario)
		 * pero sin identificaciones cruzadas */
    	readPPM(false, true, 1875.0);
    	int PPMstars = getPPMStars();
    	struct PPMstar_struct *PPMstar = getPPMStruct();
	}

	// Lee Catálogo General Argentino
	stream = fopen("cat/gc.txt", "rt");
	if (stream == NULL) {
		perror("Cannot read gc.txt");
		exit(1);
	}

	if (mode) {
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
			readFieldSanitized(buffer, cell, 16, 2);
			int RAh = atoi(cell);
			double RA = (double) RAh;
			readFieldSanitized(buffer, cell, 18, 2);
			int RAm = atoi(cell);
			RA += ((double) RAm)/60.0;
			readFieldSanitized(buffer, cell, 20, 4);
			int RAs = atoi(cell);
			RA += (((double) RAs)/100.0)/3600.0;
			RA *= 15.0; /* conversion horas a grados */

			/* lee declinacion B1875.0 */
			readFieldSanitized(buffer, cell, 39, 2);
			int Decld = atoi(cell);
			double Decl = (double) Decld;
			readFieldSanitized(buffer, cell, 41, 2);
			int Declm = atoi(cell);
			Decl += ((double) Declm)/60.0;
			readFieldSanitized(buffer, cell, 43, 3);
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
			GCstar[GCstars].distYarnall = MAX_DISTANCE;
			GCstar[GCstars].weissRef = -1;
			GCstar[GCstars].oeltzenRef = -1;
			GCstar[GCstars].distWeiss = MAX_DISTANCE;
			GCstar[GCstars].stoneRef = -1;
			GCstar[GCstars].lacailleRef = -1;
			GCstar[GCstars].distStone = MAX_DISTANCE;
			GCstar[GCstars].giRef = -1;
			GCstar[GCstars].distGi = MAX_DISTANCE;

			/* proxima estrella */
			GCstars++;
			//printf("Pos %d: id=%d RA=%.4f Decl=%.4f (%.2f) Vmag=%.1f\n", GCstars, gcRef, RA, Decl, epoch, vmag);
		}

		/* descartamos aquellas dobles cercanas débiles (cuya estrella CD es la misma) */
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
	} else {
		// Genera estrella ficticia GC 1

		/* lee ascension recta B1875.0 */
		double RA = (double) fictRAh;
		RA += ((double) fictRAm)/60.0;
		RA += (((double) fictRAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1875.0 */
		double Decl = (double) fictDecld;
		Decl += ((double) fictDeclm)/60.0;
		Decl += (((double) fictDecls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */

		/* calcula coordenadas rectangulares) */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);

		/* Busca la estrella asociada en DM más cercana */
		double minDistance = HUGE_NUMBER;
		int cdIndex;
		findDMByCoordinates(x, y, z, &cdIndex, &minDistance);

		/* ahora busca en el catálogo GC */
		int minGcRef = 0;
		double minGcDistance = MAX_DISTANCE;
		while (fgets(buffer, 1023, stream) != NULL) {
			/* lee numeracion */
			readField(buffer, cell, 1, 5);
			int gcRef = atoi(cell);

			/* lee ascension recta B1875.0 */
			readFieldSanitized(buffer, cell, 16, 2);
			int RAh = atoi(cell);
			double cRA = (double) RAh;
			readFieldSanitized(buffer, cell, 18, 2);
			int RAm = atoi(cell);
			cRA += ((double) RAm)/60.0;
			readFieldSanitized(buffer, cell, 20, 4);
			int RAs = atoi(cell);
			cRA += (((double) RAs)/100.0)/3600.0;
			cRA *= 15.0; /* conversion horas a grados */

			/* lee declinacion B1875.0 */
			readFieldSanitized(buffer, cell, 39, 2);
			int Decld = atoi(cell);
			double cDecl = (double) Decld;
			readFieldSanitized(buffer, cell, 41, 2);
			int Declm = atoi(cell);
			cDecl += ((double) Declm)/60.0;
			readFieldSanitized(buffer, cell, 43, 3);
			int Decls = atoi(cell);
			cDecl += (((double) Decls)/10.0)/3600.0;
			cDecl = -cDecl; /* incorpora signo negativo (en nuestro caso, siempre) */

			/* calcula coordenadas rectangulares y distancia a la ficticia */
			double cx, cy, cz;
			sph2rec(cRA, cDecl, &cx, &cy, &cz);
			double dist = 3600.0 * calcAngularDistance(x, y, z, cx, cy, cz);
			if (minGcDistance > dist) {
				minGcDistance = dist;
				minGcRef = gcRef;
			}
		}
		if (minGcRef > 0) {
			printf("Matched with GC %d at %.1f arcsec!\n", minGcRef, minGcDistance);
		}

		GCstar[GCstars].discard = false;
		GCstar[GCstars].gcRef = 1;
		GCstar[GCstars].obs = 1;
		GCstar[GCstars].RAh = fictRAh;
		GCstar[GCstars].RAm = fictRAm;
		GCstar[GCstars].RAs = fictRAs;
		GCstar[GCstars].Decld = fictDecld;
		GCstar[GCstars].Declm = fictDeclm;
		GCstar[GCstars].Decls = fictDecls;
		GCstar[GCstars].x = x;
		GCstar[GCstars].y = y;
		GCstar[GCstars].z = z;
		GCstar[GCstars].RA1875[0] = RA;
		GCstar[GCstars].Decl1875[0] = Decl;
		GCstar[GCstars].epoch[0] = 1875.0;
		GCstar[GCstars].vmag = 0.0;
		GCstar[GCstars].page = page;
		GCstar[GCstars].cdIndex = cdIndex;
		GCstar[GCstars].dist = minDistance;
		GCstar[GCstars].cdIndexWithinMag = -1;
		GCstar[GCstars].distWithinMag = HUGE_NUMBER;
		GCstar[GCstars].yarnallRef = -1;
		GCstar[GCstars].distYarnall = MAX_DISTANCE;
		GCstar[GCstars].weissRef = -1;
		GCstar[GCstars].oeltzenRef = -1;
		GCstar[GCstars].distWeiss = MAX_DISTANCE;
		GCstar[GCstars].stoneRef = -1;
		GCstar[GCstars].lacailleRef = -1;
		GCstar[GCstars].distStone = MAX_DISTANCE;
		GCstar[GCstars].giRef = -1;
		GCstar[GCstars].distGi = MAX_DISTANCE;
		GCstars++;
	}
	fclose(stream);

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
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */
    
	    /* convierte coordenadas a 1875.0 */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(1860.0, 1875.0, &RA1875, &Decl1875);
    	if (mode && Decl1875 > -22) continue;

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
		if (minDistance < GCstar[gcIndex].distYarnall) {
			copy(GCstar[gcIndex].yarnallCat, cell);
			GCstar[gcIndex].yarnallRef = yarnallRef;
			GCstar[gcIndex].distYarnall = minDistance;
			GCstar[gcIndex].vmagYarnall = vmag;
			countUSNO++;
		}
		if (mode && CD_SEARCH) {
			/* aprovechamos a revisar CD con esta estrella */
			minDistance = HUGE_NUMBER;
			int cdIndex;
			findDMByCoordinates(x, y, z, &cdIndex, &minDistance);
			if (minDistance > MAX_DISTANCE) {
				char yarnallCat[28];
				copy(yarnallCat, cell);
				printf("  USNO %d <%s> is ALONE with mag=%.1f; nearest CD %dº%d separated in %.1f arcsec.\n",
					yarnallRef,
					yarnallCat,
					vmag,
					CDstar[cdIndex].declRef,
					CDstar[cdIndex].numRef,
					minDistance
				);
				printf("    RA %02dh%02dm%02ds%02d DE %02d°%02d'%02d''%01d\n",
					RAh, RAm, RAs / 100, RAs % 100,
					Decld, Declm, Decls / 10, Decls % 10);
			}
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
    
	    /* convierte coordenadas a 1875.0 */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(1850.0, 1875.0, &RA1875, &Decl1875);
    	if (mode && Decl1875 > -22) continue;

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
		if (minDistance < GCstar[gcIndex].distWeiss) {
			GCstar[gcIndex].weissRef = weissRef;
			GCstar[gcIndex].oeltzenRef = oeltzenRef;
			GCstar[gcIndex].distWeiss = minDistance;
			GCstar[gcIndex].vmagWeiss = vmag;
			countWeiss++;
		}
		if (mode && CD_SEARCH) {
			/* aprovechamos a revisar CD con esta estrella */
			minDistance = HUGE_NUMBER;
			int cdIndex;
			findDMByCoordinates(x, y, z, &cdIndex, &minDistance);
			if (minDistance > MAX_DISTANCE) {
				printf("  WEI %d or OA %d is ALONE with mag=%.1f; nearest CD %dº%d separated in %.1f arcsec.\n",
					weissRef,
					oeltzenRef,
					vmag,
					CDstar[cdIndex].declRef,
					CDstar[cdIndex].numRef,
					minDistance
				);
				printf("    RA %02dh%02dm%02ds%02d DE %02d°%02d'%02d''%01d\n",
					RAh, RAm, RAs / 100, RAs % 100,
					Decld, Declm, Decls / 10, Decls % 10);
			}
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

	    /* convierte coordenadas a 1875.0 según cada mean time */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(1880.0, 1875.0, &RA1875, &Decl1875);
    	if (mode && Decl1875 > -22) continue;

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
		if (minDistance < GCstar[gcIndex].distStone) {
			GCstar[gcIndex].stoneRef = stoneRef;
			GCstar[gcIndex].lacailleRef = lacailleRef;
			GCstar[gcIndex].distStone = minDistance;
			GCstar[gcIndex].vmagStone = vmag;
			countStone++;
		}
		if (mode && CD_SEARCH) {
			/* aprovechamos a revisar CD con esta estrella */
			minDistance = HUGE_NUMBER;
			int cdIndex;
			findDMByCoordinates(x, y, z, &cdIndex, &minDistance);
			if (minDistance > MAX_DISTANCE) {
				printf("  ST %d is ALONE with mag=%.1f; nearest CD %dº%d separated in %.1f arcsec.\n",
					stoneRef,
					vmag,
					CDstar[cdIndex].declRef,
					CDstar[cdIndex].numRef,
					minDistance
				);
			}
		}
	}
    fclose(stream2);
    fclose(stream);
	printf("  cross-referenced with %d stars\n", countStone);

	/* añadimos referencia al catálogo de Gilliss */
	printf("Reading Gilliss catalog...\n");
    stream = fopen("cat/gilliss.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read gilliss.txt");
		exit(1);
    }
	int countGi = 0;
	int countGiCD = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
		/* lee signo declinación y descarta hemisferio norte */
		readField(buffer, cell, 48, 1);
		if (cell[0] != '-') continue;

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
    
	    /* convierte coordenadas a 1875.0 */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(1850.0, 1875.0, &RA1875, &Decl1875);
    	if (mode && Decl1875 > -22) continue;

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
		if (minDistance < GCstar[gcIndex].distGi) {
			copy(GCstar[gcIndex].giCat, cell);
			GCstar[gcIndex].giRef = giRef;
			GCstar[gcIndex].distGi = minDistance;
			GCstar[gcIndex].vmagGi = vmag;
			countGi++;
		}
		if (mode && CD_SEARCH) {
			/* aprovechamos a revisar CD con esta estrella
			 * (tambien mencionamos el numero de GZC en caso de haberlo)*/
			char gouldZC[30];
			gouldZC[0] = 0;
			if (cell[9] == 'Z' && cell[10] == 'C') {
				char hour[3], num[6];
				hour[0] = cell[11]; hour[1] = cell[12]; hour[2] = 0;
				int gouldHour = atoi(hour);
				num[0] = cell[13]; num[1] = cell[14]; num[2] = cell[15];
				num[3] = cell[16]; num[4] = cell[17]; num[5] = 0;
				int gouldNum = atoi(num);
				snprintf(gouldZC, 30, " or GZC %dh %d", gouldHour, gouldNum);
			}
			minDistance = HUGE_NUMBER;
			int cdIndex;
			findDMByCoordinates(x, y, z, &cdIndex, &minDistance);
			if (minDistance > MAX_DISTANCE) {
				countGiCD++;
				printf("  %d> GI %d%s is ALONE with mag=%.1f; nearest CD %dº%d separated in %.1f arcsec.\n",
					countGiCD,
					giRef,
					gouldZC,
					vmag,
					CDstar[cdIndex].declRef,
					CDstar[cdIndex].numRef,
					minDistance
				);
				printf("    RA %02dh%02dm%02ds%02d DE %02d°%02d'%02d''%01d\n",
					RAh, RAm, RAs / 100, RAs % 100,
					Decld, Declm, Decls / 10, Decls % 10);
			}
		}
	}
    fclose(stream);
	printf("  cross-referenced with %d stars\n", countGi);

	if (mode) {	
		printf("Stars read from Catalogo General Argentino: %d  (discarded = %d)\n", GCstars, discarded);
	}
}
