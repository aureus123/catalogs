
/*
 * COMPARE_SD - Compara registros de CD contra SD a través de una identificación cruzada
 * Copyright (C) 2024 Daniel E. Severin
 *
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
#define MAX_MAGNITUDE 0.8
#define MAXSDSTAR 10000

/* Para uso de la libreria WCS: */
#define WCS_B1950 2 /* B1950(FK4) right ascension and declination */
extern "C" void wcsconp(int sys1, int sys2, double eq1, double eq2, double ep1, double ep2,
             double *dtheta, double *dphi, double *ptheta, double *pphi);


struct SDstar_struct {
    bool discard; /* true if should not be considered */
    int numRef; /* identificador con numero (siempre decl = -22) */
    float rah, ramin, raseg, decldeg, declmin; /* Ascension recta y declinacion, en partes */
    double RA1855, Decl1855, RA1875, Decl1875, vmag; /* coordenadas y magnitud visual */
    int cdIndex; /* Indice a CD */
    double x, y, z; /* coordenadas rectangulares en circulo unidad */
    double dist; /* distancia angular a su CD asociada (en arcsec) */
} SDstar[MAXSDSTAR];

int SDstars;

/*
 * writeRegisterSD - escribe en pantalla un registro de SD en su formato
 */
void writeRegisterSD(int sdIndex)
{
    struct SDstar_struct *s = &SDstar[sdIndex];

    int numRef = s->numRef;
    int page = 438 + 22 * sdIndex / SDstars;
    printf("     Register SD -22°%d (en %.0fh):  %.1f | %.0fm%.1fs | %.1f'     (pag. %d)\n",
                numRef,
                s->rah,
                s->vmag,
                s->ramin,
                s->raseg,
                s->declmin,
                page);
}

/*
 * lee estrellas del Southern Durchmusterung, solo decl -22 (mismo formato que BD, ver read_bd.cpp).
 */
void readSD() {
    FILE *stream;
    char buffer[1024];
    char cell[256];

    int CDstars = getDMStars();
    struct DMstar_struct *CDstar = getDMStruct();

    stream = fopen("cat/sd.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read sd.txt");
        exit(1);
    }

    SDstars = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
        /* lee la zona de declinacion */
        readField(buffer, cell, 3, 3);
        int declRef = atoi(cell);
        if (declRef != -22) continue;

        /* lee numeracion y caracter suplementario */
        readField(buffer, cell, 6, 5);
        int numRef = atoi(cell);
        char supplRef = buffer[11-1];
        if (supplRef == 'D') continue;
        if (supplRef != ' ') {
            printf("Ommitting star SD %d°%d%c\n", declRef, numRef, supplRef);
            continue;
        }

        /* lee magnitud visual */
        readField(buffer, cell, 12, 4);
        double vmag = atof(cell);
        if (vmag > 12.1) {
            /* vmag no es una magnitud, si no un codigo:
             * 20.0 = neb
             * 30.0 = var
             * 40.0 = nova or nova?
             * 50.0 = cluster */
            if ((vmag > 19.9 && vmag < 20.1) || (vmag > 39.9 && vmag < 40.1) || (vmag > 49.9 && vmag < 50.1)) continue;
            if (vmag > 29.9 && vmag < 30.1) {
                /* estrella variable */
            } else {
                printf("Unknown code: %f, SD %d°%d\n", vmag, declRef, numRef);
                exit(1);
            }
        }

        /* lee ascension recta B1855.0 */
        readField(buffer, cell, 16, 2);
        float rah = atof(cell);
        double RA = rah;
        readField(buffer, cell, 18, 2);
        float ramin = atof(cell);
        RA += ramin/60.0;
        readField(buffer, cell, 20, 4);
        float raseg = atof(cell);
        RA += raseg/3600.0;
        RA *= 15.0; /* conversion horas a grados */

        /* lee declinacion B1855.0 */
        readField(buffer, cell, 25, 2);
        float decldeg = atof(cell);
        double Decl = decldeg;
        readField(buffer, cell, 27, 6);
        float declmin = atof(cell);
        Decl += declmin/60.0;
        Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */

        /* calcula coordenadas rectangulares, pero de 1875.0 */
        double RA1875 = RA;
        double Decl1875 = Decl;
        double pmRA = 0.0;
        double pmDecl = 0.0;
        wcsconp(WCS_B1950, WCS_B1950, 1855.0, 1875.0, 1855.0, 1875.0, &RA1875, &Decl1875, &pmRA, &pmDecl);
        double x, y, z;
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        /* la almacena en memoria */
		if (SDstars == MAXSDSTAR) bye("Maximum amount reached!\n");
        SDstar[SDstars].discard = true; /* todas comienzan descartadas a menos que se cruce */
        SDstar[SDstars].numRef = numRef;
        SDstar[SDstars].rah = rah;
        SDstar[SDstars].ramin = ramin;
        SDstar[SDstars].raseg = raseg;
        SDstar[SDstars].decldeg = decldeg;
        SDstar[SDstars].declmin = declmin;
        SDstar[SDstars].RA1855 = RA;
        SDstar[SDstars].Decl1855 = Decl;
        SDstar[SDstars].RA1875 = RA1875;
        SDstar[SDstars].Decl1875 = Decl1875;
        SDstar[SDstars].vmag = vmag;
        SDstar[SDstars].cdIndex = -1; /* a ser rellenado en la siguiente fase */
        SDstar[SDstars].x = x;
        SDstar[SDstars].y = y;
        SDstar[SDstars].z = z;
        SDstar[SDstars].dist = 0.0; /* a ser rellenado en la siguiente fase */
        SDstars++;
    }
    printf("Stars read from Southern Durchmusterung: %d\n", SDstars);
    fclose(stream);

    /* siguiente fase: leer identificación cruzada del catálogo 4005 */
    stream = fopen("cat/4005.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read 4005.txt");
        exit(1);
    }

    int crossed = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
        /* necesitamos que una sea SD y la otra sea CD */
        readField(buffer, cell, 12, 1);
        bool isSourceSD = cell[0] == ' ';
        bool isSourceCD = cell[0] == '2';
        readField(buffer, cell, 25, 1);
        bool isTargetSD = cell[0] == ' ';
        bool isTargetCD = cell[0] == '2';
        int declRefSource, numRefSource, declRefTarget, numRefTarget;
        if (isSourceSD && isTargetCD) {
            declRefSource = 2;
            numRefSource = 6;
            declRefTarget = 15;
            numRefTarget = 19;
        } else if (isSourceCD && isTargetSD) {
            declRefSource = 15;
            numRefSource = 19;
            declRefTarget = 2;
            numRefTarget = 6;
        } else continue;

        /* lee la zona de declinacion de SD y numero */
        readField(buffer, cell, declRefSource, 3);
        int declRefSD = atoi(cell);
        if (declRefSD != -22) continue;
        readField(buffer, cell, numRefSource, 5);
        int numRefSD = atoi(cell);
        if (numRefSD == 0) bye("Error in SD num!");

        /* lee la estrella de CD asociada, si existe */
        readField(buffer, cell, declRefTarget, 3);
        int declRefCD = atoi(cell);
        readField(buffer, cell, numRefTarget, 5);
        int numRefCD = atoi(cell);
        if (numRefCD == 0) bye("Error in CD num!");

        /* buscar la estrella en el catalogo SD */
        int catIndex = -1;
        for (int i = 0; i < SDstars; i++) {
            if (numRefSD == SDstar[i].numRef) {
                catIndex = i;
                break;
            }
        }
        if (catIndex == -1) {
            printf("Warning: Star SD %d°%d not found\n", declRefSD, numRefSD);
            continue;
        }

        /* obtiene posición en la esfera unidad */
        double x = SDstar[catIndex].x;
        double y = SDstar[catIndex].y;
        double z = SDstar[catIndex].z;

        /* Busca la estrella asociada en CD; en caso de haber más
         * de una, escoger la de menor distancia */
        int cdIndex = -1;
        double minDistance = 9999999999;
        int i = getDMindex(true, declRefCD, numRefCD);
        while (i != -1) {
          double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
          if (minDistance > dist) {
            cdIndex = i;
            minDistance = dist;
          }
          i = CDstar[i].next;
        }
        if (cdIndex == -1) {
          printf("CD %d°%d not found (corresponding to SD %d°%d). Discarding SD star.\n",
            declRefCD, numRefCD, declRefSD, numRefSD);
          continue;
        }

        /* se almacena la lectura de la identificacion */
        SDstar[catIndex].discard = false;
        SDstar[catIndex].cdIndex = cdIndex;
        SDstar[catIndex].dist = minDistance;
        CDstar[cdIndex].catIndex = catIndex;
        crossed++;
    }
    printf("Number of SD stars cross-identified with CD stars: %d\n", crossed);
    fclose(stream);
}

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("COMPARE_SD - Compare CD and SD catalogs.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    /* leemos catalogo CD */
    readDM("cat/cd.txt");
    struct DMstar_struct *CDstar = getDMStruct();

    /* leemos catalogo SD */
    readSD();

    /* revisamos la identificación cruzada y generamos dos planillas */
    FILE *posStream, *magStream;

    posStream = fopen("results/table_pos_sd.csv", "wt");
    if (posStream == NULL) {
        perror("Cannot write in table_pos_sd.csv");
        exit(1);
    }
    fprintf(posStream, "index,decl,num,ref,dist10\n");

    magStream = fopen("results/table_mag_sd.csv", "wt");
    if (magStream == NULL) {
        perror("Cannot write in table_mag_sd.csv");
        exit(1);
    }
    fprintf(magStream, "index,decl,num,ref,delta10\n");

    int dropDistError = 0;
    int maxDistError = 0;
    int magDiffError = 0;
    int indexError = 0;
    int goodStarsPosition = 0;
    double akkuDistError = 0.0;
    int goodStarsMagnitude = 0;
    double akkuDeltaError = 0.0;
    for (int i = 0; i < SDstars; i++) {
        if (SDstar[i].discard) continue;
        int cdIndex = SDstar[i].cdIndex;
        float dist = SDstar[i].dist;
        if (dist > MAX_DISTANCE_DROP) {
            // posiciones demasiado separadas, posiblemente mala identificacion
            dropDistError++;
            printf("*) CD %d°%d too separated from SD -22°%d in %.1f arcsec.\n",
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                SDstar[i].numRef,
                dist);
            writeRegister(cdIndex);
            writeRegisterSD(i);
            continue;
        }
        if (dist > MAX_DISTANCE) {
            // posiciones muy separadas, supera umbral
            maxDistError++;
            indexError++;
            printf("%d) CD %d°%d separated from SD -22°%d in %.1f arcsec.\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                SDstar[i].numRef,
                dist);
            writeRegister(cdIndex);
            writeRegisterSD(i);
            fprintf(posStream, "%d,%d,%d,SD -22°%d,%.0f\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                SDstar[i].numRef,
                10.0 * dist);
            continue;
        }
        goodStarsPosition++;
        akkuDistError += dist * dist;

        double sdVmag = SDstar[i].vmag;
        double cdVmag = CDstar[cdIndex].vmag;
        if (sdVmag > 29.9 || cdVmag > 29.9) continue; // se omiten aquellas estrellas variables
        double delta = fabs(sdVmag - cdVmag);
        if (delta >= MAX_MAGNITUDE) {
            // diferencia en magnitud visual supera umbral
            magDiffError++;
            indexError++;
            printf("%d) CD %d°%d reports mag=%.1f but it should be mag=%.1f from SD -22°%d: Delta = %.1f.\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                cdVmag,
                sdVmag,
                SDstar[i].numRef,
                delta);
            writeRegister(cdIndex);
            writeRegisterSD(i);
            fprintf(magStream, "%d,%d,%d,SD -22°%d,%.0f\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                SDstar[i].numRef,
                10.0 * delta);
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
