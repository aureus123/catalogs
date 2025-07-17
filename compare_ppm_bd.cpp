
/*
 * COMPARE_PPM_BD - Compara registros de BD contra PPM
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

#define MAX_DISTANCE 180.0   // 3 minutos de arco
#define MAX_MAGNITUDE 0.4
#define CURATED false // true if curated BD catalog should be used

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("COMPARE_PPM_BD - Compare BD and PPM catalogs.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    char ppmName[20];

    /* leemos catalogo BD (primer volumen) */
    readDM(CURATED ? "cat/bd_vol1_curated.txt" : "cat/bd_vol1.txt");
    struct DMstar_struct *BDstar = getDMStruct();

    /* leemos catalogo PPM */
    readPPM(true, true, false, false, 1855.0);
    int PPMstars = getPPMStars();
    struct PPMstar_struct *PPMstar = getPPMStruct();

    /* revisamos la identificación cruzada y generamos dos planillas */
    FILE *posStream, *magStream;
    posStream = openPositionFile("results/table_pos_bd.csv");
    magStream = openMagnitudeFile("results/table_mag_bd.csv");

    int maxDistError = 0;
    int magDiffError = 0;
    int totalErrorsMinusDoubles = 0;
    int indexError = 0;
    int problematic = 0;
    int goodStarsPosition = 0;
    double akkuDistError = 0.0;
    int goodStarsMagnitude = 0;
    double akkuDeltaError = 0.0;
    for (int i = 0; i < PPMstars; i++) {
        if (PPMstar[i].discard == true) continue;
        bool isProblematic = false;
        if (PPMstar[i].problem == 1) {
            // se avisa si la estrella es "problematica"
            isProblematic = true;
            problematic++;
        }
        int bdIndex = PPMstar[i].dmIndex;
        char *bdString = PPMstar[i].dmString;
        double dist = PPMstar[i].dist;
        if (dist > MAX_DISTANCE) {
            // posiciones muy separadas, supera umbral
            maxDistError++;
            indexError++;
            printf("%d) %s separated from PPM %d%s in %.1f arcsec.\n",
                indexError,
                bdString,
                PPMstar[i].ppmRef,
                isProblematic ? " (PROBLEM)" : "",
                dist);
            writeRegister(bdIndex, true);
            if (!revise(i)) totalErrorsMinusDoubles++;
            snprintf(ppmName, 20, "PPM %d", PPMstar[i].ppmRef);
            writePositionEntry(posStream,
                indexError,
                BDstar[bdIndex].declRef,
                BDstar[bdIndex].numRef,
                ppmName,
                dist);
            continue;
        }
        goodStarsPosition++;
        akkuDistError += dist * dist;

        double ppmVmag = PPMstar[i].vmag;
        double bdVmag = BDstar[bdIndex].vmag;
        if (ppmVmag < 0.00001 || bdVmag > 29.9) continue; // se omiten aquellas estrellas con Vmag=0 o variables
        ppmVmag = compVmagToCDmag(BDstar[bdIndex].declRef, ppmVmag);
        double delta = fabs(ppmVmag - bdVmag);
        if (delta >= MAX_MAGNITUDE) {
            // diferencia en magnitud V y visual supera umbral
            magDiffError++;
            indexError++;
            printf("%d) %s reports mag=%.1f but it should be mag=%.1f from PPM %d%s: Delta = %.1f.\n",
                indexError,
                bdString,
                bdVmag,
                ppmVmag,
                PPMstar[i].ppmRef,
                isProblematic ? " (PROBLEM)" : "",
                delta);
            writeRegister(bdIndex, false);
            if (!revise(i)) totalErrorsMinusDoubles++;
            snprintf(ppmName, 20, "PPM %d", PPMstar[i].ppmRef);
            writeMagnitudeEntry(magStream,
                indexError,
                BDstar[bdIndex].declRef,
                BDstar[bdIndex].numRef,
                ppmName,
                delta);
            continue;
        }
        goodStarsMagnitude++;
        akkuDeltaError += delta * delta;
    }
    fclose(magStream);
    fclose(posStream);
    printf("Total errors: %d (position: %d, mag: %d); errors without warning = %d, PPM with problems = %d\n",
        indexError, maxDistError, magDiffError, totalErrorsMinusDoubles,  problematic);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)goodStarsPosition),
        goodStarsPosition);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)goodStarsMagnitude),
        goodStarsMagnitude);
    return 0;
}
