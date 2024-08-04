
/*
 * COMPARE_GC - Compara registros de CD contra GC
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
#define MAX_MAGNITUDE 1.0

struct GCstar_struct {
	bool discard; /* true si debe ser descartada (por doble) */
    int gcRef; /* identificador con numero */
    int obs; /* cantidad de observaciones */
	int RAh, RAm, RAs, Decld, Declm, Decls; /* detalle primera observacion */
    double RA1875[10], Decl1875[10], epoch[10]; /* coordenadas y epoca de la observacion */
    double vmag; /* magnitud visual */
    int page; /* pagina donde se encuentra */
    int cdIndex; /* indice (no identificador) a la estrella CD más cercana */
	double dist; /* distancia a la estrella CD */
	int cdIndexWithinMag; /* indice a la CD más cercana, pero dentro del rango de magnitud */
	double distWithinMag;
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
        double minDistance = 9999999999;
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
	    GCstar[GCstars].RA1875[0] = RA;
	    GCstar[GCstars].Decl1875[0] = Decl;
	    GCstar[GCstars].epoch[0] = epoch;
	    GCstar[GCstars].vmag = vmag;
	    GCstar[GCstars].page = page;
	    GCstar[GCstars].cdIndex = cdIndex;
	    GCstar[GCstars].dist = minDistance;
		GCstar[GCstars].cdIndexWithinMag = cdIndexWithinMag;
		GCstar[GCstars].distWithinMag = minDistWithinMag;

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
            printf("%d) GC %d is alone within MAG = %.1f; similar CD %d°%d separated in %.1f arcsec.\n",
                indexError,
                GCstar[i].gcRef,
                gcVmag,
                CDstar[cdIndexWithinMag].declRef,
                CDstar[cdIndexWithinMag].numRef,
                distWithinMag);
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
