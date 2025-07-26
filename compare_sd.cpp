
/*
 * COMPARE_SD - Compara registros de CD contra SD a través de una identificación cruzada
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "read_sd.h"
#include "trig.h"
#include "misc.h"

#define MAX_DISTANCE 180.0   // 3 minutos de arco
#define MAX_DISTANCE_DROP 3600.0 // 1 grado de arco
#define MAX_MAGNITUDE 0.8

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("COMPARE_SD - Compare CD and SD catalogs.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    char sdName[20];

    /* leemos catalogo CD */
    readDM("cat/cd.txt");
    struct DMstar_struct *CDstar = getDMStruct();

    /* leemos catalogo SD */
    readSD(true);
    int SDstars = getSDStars();
    struct SDstar_struct *SDstar = getSDStruct();

    /* revisamos la identificación cruzada y generamos dos planillas */
    FILE *posStream, *magStream;
    posStream = openPositionFile("results/table_pos_sd.csv");
    magStream = openMagnitudeFile("results/table_mag_sd.csv");

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
            writeRegister(cdIndex, false);
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
            writeRegister(cdIndex, true);
            writeRegisterSD(i);
            snprintf(sdName, 20, "SD -22°%d", SDstar[i].numRef);
            writePositionEntry(posStream,
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                sdName,
                dist);
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
            writeRegister(cdIndex, false);
            writeRegisterSD(i);
            snprintf(sdName, 20, "SD -22°%d", SDstar[i].numRef);
            writeMagnitudeEntry(magStream,
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                sdName,
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
