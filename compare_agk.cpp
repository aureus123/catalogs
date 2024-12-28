
/*
 * COMPARE_AGK - Compara registros de CD contra AGK
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "trig.h"
#include "misc.h"

#define MAX_DISTANCE 180.0   // 3 minutos de arco
#define MAX_DISTANCE_DROP 3600.0 // 1 grado de arco
#define MAX_MAGNITUDE 1.0
#define MAXAGKSTAR 50000

struct AGKstar_struct {
    int agkRef; /* identificador */
    char letter; /* Cordoba A, B o C */
    double RAeq, Decleq, RA1875, Decl1875, vmag; /* magnitud y posiciones */
    int cdIndex; /* indice a CD */
    int cdIndexWarning; /* otro indice a CD en caso que no coincida declinacion, o -1 */
    double x, y, z;
    double dist;
} AGKstar[MAXAGKSTAR];

int AGKstars;

/*
 * lee estrellas del Astronomische Gesellschaft Katalog, coordenadas segun equinox que las pasa a 1875.0
 * cada entrada se corresponde con una estrella de CD, cuya declinacion varia entre declA y declB
 * se deben dar las posiciones y cantidad de bytes a leer del archivo
 */
void readAGK(
    const char *filename,
    int ref_pos,
    int ref_bytes,
    int vmag_pos,
    int vmag_bytes,
    int RA_pos,
    int Decl_pos,
    int durch_pos,
    int durch_bytes,
	int cpdmark_pos,
    char letter)
{
    FILE *stream;
    char buffer[1024], cell[256];

    int CDstars = getDMStars();
    struct DMstar_struct *CDstar = getDMStruct();

    stream = fopen(filename, "rt");
    if (stream == NULL) {
	    snprintf(buffer, 1024, "Cannot read %s", filename);
        perror(buffer);
	    exit(1);
    }

    printf("Reading AGK %c:\n", letter);
    while (fgets(buffer, 1023, stream) != NULL) {
    	/* lee numeracion */
	    readField(buffer, cell, ref_pos, ref_bytes);
        int agkRef = atoi(cell);

        /* lee otros campos (evita usar estrellas de CPD) */
		readField(buffer, cell, cpdmark_pos, 1);
		if (cell[0] == '-') continue;
        readField(buffer, cell, durch_pos, durch_bytes);
	    if (cell[0] == '&') continue;
	    int numRef = atoi(cell);

    	/* lee magnitud */
	    readField(buffer, cell, vmag_pos, vmag_bytes);
	    double vmag = 0.0;
    	if (cell[0] != '&' && cell[1] != '&') vmag = atof(cell)/10.0;

	    /* lee ascension recta de epoca */
	    readField(buffer, cell, RA_pos, 2);
        double RA = atof(cell);
	    readField(buffer, cell, RA_pos+2, 2);
        RA += atof(cell)/60.0;
	    readField(buffer, cell, RA_pos+4, 4);
        RA += (atof(cell)/100.0)/3600.0;
	    RA *= 15.0; /* conversion horas a grados */

	    /* lee declinacion de epoca */
	    readField(buffer, cell, Decl_pos, 2);
        double Decl = atof(cell);
	    readField(buffer, cell, Decl_pos+2, 2);
        Decl += atof(cell)/60.0;
	    readField(buffer, cell, Decl_pos+4, 3);
        Decl += (atof(cell)/10.0)/3600.0;
	    Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */

	    /* convierte coordenadas a 1875.0 */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(1900.0, 1875.0, &RA1875, &Decl1875);

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RA1875, Decl1875, &x, &y, &z);

	    /* computa la estrella del Durchmusterung correspondiente
         * pero hay que revisar también las declinaciones adyacentes */
	    int cdIndex = -1;
        int cdIndexWarning = -1;
        double minDistance = 9999999999;
	    int declRef = -((int)floor(-Decl1875));
        int declFinal;
        int i = getDMindex(true, declRef, numRef);
        while (i != -1) {
          double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
          if (minDistance > dist) {
            cdIndex = i;
            cdIndexWarning = i;
            minDistance = dist;
          }
          i = CDstar[i].next;
        }
        i = getDMindex(true, declRef - 1, numRef);
        while (i != -1) {
          double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
          if (minDistance > dist) {
            cdIndex = i;
            minDistance = dist;
          }
          i = CDstar[i].next;
        }
        i = getDMindex(true, declRef + 1, numRef);
        while (i != -1) {
          double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
          if (minDistance > dist) {
            cdIndex = i;
            minDistance = dist;
          }
          i = CDstar[i].next;
        }

        if (cdIndex == -1) {
          printf("Star DM %d not found (corresponding to AGK %d). Discarding star.\n", numRef, agkRef);
          continue;
        }    
        if (CDstar[cdIndex].declRef != declRef) {
            printf("Warning: Decl = %d of AGK %d does not coincide with CD %d°%d\n",
                declRef, agkRef, CDstar[cdIndex].declRef, numRef);
        } else {
            cdIndexWarning = -1;
        }

    	/* la almacena en memoria */
		if (AGKstars == MAXAGKSTAR) bye("Maximum amount reached!\n");
    	AGKstar[AGKstars].agkRef = agkRef;
        AGKstar[AGKstars].letter = letter;
	    AGKstar[AGKstars].RAeq = RA;
	    AGKstar[AGKstars].Decleq = Decl;
	    AGKstar[AGKstars].RA1875 = RA1875;
	    AGKstar[AGKstars].Decl1875 = Decl1875;
	    AGKstar[AGKstars].vmag = vmag;
	    AGKstar[AGKstars].cdIndex = cdIndex;
        AGKstar[AGKstars].cdIndexWarning = cdIndexWarning;
	    AGKstar[AGKstars].x = x;
	    AGKstar[AGKstars].y = y;
	    AGKstar[AGKstars].z = z;
	    AGKstar[AGKstars].dist = minDistance;
    	/* proxima estrella */
	    AGKstars++;
    }
    fclose(stream);
}

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("COMPARE_AGK - Compare CD and AGK catalogs.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    char agkName[20];

    /* leemos catalogo CD */
    readDM("cat/cd.txt");
    struct DMstar_struct *CDstar = getDMStruct();

    /* leemos catalogos Cordoba A, B y C */
    AGKstars = 0;
    readAGK("cat/corda.txt", 8, 5, 13, 2, 15, 36, 57, 5, 62, 'A');
    readAGK("cat/cordb.txt", 8, 5, 13, 3, 16, 37, 60, 5, 66, 'B');
    readAGK("cat/cordc.txt", 8, 5, 13, 2, 15, 36, 59, 5, 65, 'C');
    printf("Stars read from AGK: %d\n", AGKstars);

    /* revisamos la identificación cruzada y generamos dos planillas */
    FILE *posStream, *magStream;
    posStream = openPositionFile("results/table_pos_agk.csv");
    magStream = openMagnitudeFile("results/table_mag_agk.csv");

    int dropDistError = 0;
    int maxDistError = 0;
    int magDiffError = 0;
    int indexError = 0;
    int goodStarsPosition = 0;
    double akkuDistError = 0.0;
    int goodStarsMagnitude = 0;
    double akkuDeltaError = 0.0;
    for (int i = 0; i < AGKstars; i++) {
        int cdIndex = AGKstar[i].cdIndex;
        float dist = AGKstar[i].dist;
        if (dist > MAX_DISTANCE_DROP) {
            // posiciones demasiado separadas, posiblemente mala identificacion
            dropDistError++;
            printf("*) CD %d°%d too separated from AGK %c%d in %.1f arcsec.\n",
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                AGKstar[i].letter,
                AGKstar[i].agkRef,
                dist);
            writeRegister(cdIndex, false);
            continue;
        }
        if (dist > MAX_DISTANCE) {
            // posiciones muy separadas, supera umbral
            maxDistError++;
            indexError++;
            printf("%d) CD %d°%d separated from AGK %c%d in %.1f arcsec.\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                AGKstar[i].letter,
                AGKstar[i].agkRef,
                dist);
            writeRegister(cdIndex, true);
            int cdIndexWarning = AGKstar[i].cdIndexWarning;
            if (cdIndexWarning != -1) {
                double distWarning = 3600.0 * calcAngularDistance(
                    AGKstar[i].x, AGKstar[i].y, AGKstar[i].z,
                    CDstar[cdIndexWarning].x, CDstar[cdIndexWarning].y, CDstar[cdIndexWarning].z);
                printf("     Also check CD %d°%d (dist = %.1f arcsec.)\n",
                    CDstar[cdIndexWarning].declRef,
                    CDstar[cdIndexWarning].numRef,
                    distWarning);
                writeRegister(cdIndexWarning, false);
            }
            snprintf(agkName, 20, "AGK %c%d", AGKstar[i].letter, AGKstar[i].agkRef);
            writePositionEntry(posStream,
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                agkName,
                dist);
            continue;
        }
        goodStarsPosition++;
        akkuDistError += dist * dist;

        double agkVmag = AGKstar[i].vmag;
        double cdVmag = CDstar[cdIndex].vmag;
        if (agkVmag < 0.00001 || cdVmag > 29.9) continue; // se omiten aquellas estrellas con Vmag=0 o variables
        double delta = fabs(agkVmag - cdVmag);
        if (delta >= MAX_MAGNITUDE) {
            // diferencia en magnitud visual supera umbral
            magDiffError++;
            indexError++;
            printf("%d) CD %d°%d reports mag=%.1f but it should be mag=%.1f from AGK %c%d: Delta = %.1f.\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                cdVmag,
                agkVmag,
                AGKstar[i].letter,
                AGKstar[i].agkRef,
                delta);
            writeRegister(cdIndex, false);
            snprintf(agkName, 20, "AGK %c%d", AGKstar[i].letter, AGKstar[i].agkRef);
            writeMagnitudeEntry(magStream,
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                agkName,
                delta);
            continue;
        }
        goodStarsMagnitude++;
        akkuDeltaError += delta * delta;
    }
    fclose(magStream);
    fclose(posStream);
    printf("Total errors: %d (position: %d, mag: %d), Drops: %d\n",
        indexError, maxDistError, magDiffError, dropDistError);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)goodStarsPosition),
        goodStarsPosition);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)goodStarsMagnitude),
        goodStarsMagnitude);
    return 0;
}
