
/*
 * COMPARE_PPM - Compara registros de CD contra PPM
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

#define MAX_DISTANCE 120.0   // 2 minutos de arco
#define MAX_MAGNITUDE 1.5

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("COMPARE_PPM - Compare CD and PPM catalogs.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    if (argc < 2) {
        printf("Usage: compare_ppm file\n");
        printf("    where file can be:\n");
        printf("        cd.txt = Current CD catalog at Vizier\n");
        printf("        1114.txt = NASA-ADC CD catalog (has some errors)\n");
        printf("        cd_vol1.txt = same as cd.txt but only 1st. Volume (Resultados XVI)\n");
        printf("        I88.txt = 1982 CD catalog version (has some errors)\n");
        exit(-1);
    }
    bool allSky = true;
    if (strcmp(argv[1], "cd.txt") && strcmp(argv[1], "1114.txt")) {
        if (strcmp(argv[1], "cd_vol1.txt") && strcmp(argv[1], "I88.txt")) {
            printf("Bad file name. See usage.\n");
            exit(-1);
        }
        allSky = false;
    }

    /* leemos catalogo CD */
    char buffer[64];
    snprintf(buffer, 64, "cat/%s", argv[1]);
    readDM(buffer);
    struct DMstar_struct *CDstar = getDMStruct();

    /* leemos catalogo PPM */
    readPPM(allSky);
    int PPMstars = getPPMStars();
    struct PPMstar_struct *PPMstar = getPPMStruct();

    /* revisamos la identificación cruzada y generamos dos planillas */
    FILE *posStream, *magStream;

    posStream = fopen("results/table_pos_ppm.csv", "wt");
    if (posStream == NULL) {
        perror("Cannot write in table_pos_cd.csv");
        exit(1);
    }
    fprintf(posStream, "index,decl,num,ref,dist10\n");

    magStream = fopen("results/table_mag_ppm.csv", "wt");
    if (magStream == NULL) {
        perror("Cannot write in table_mag_cd.csv");
        exit(1);
    }
    fprintf(magStream, "index,decl,num,ref,delta10\n");

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
        int cdIndex = PPMstar[i].dmIndex;
        float dist = PPMstar[i].dist;
        if (dist > MAX_DISTANCE) {
            // posiciones muy separadas, supera umbral
            maxDistError++;
            indexError++;
            printf("%d) CD %d°%d separated from PPM %d%s in %.1f arcsec.\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                PPMstar[i].ppmRef,
                isProblematic ? " (PROBLEM)" : "",
                dist);
            writeRegister(cdIndex);
            if (!revise(i)) totalErrorsMinusDoubles++;
            fprintf(posStream, "%d,%d,%d,PPM %d,%.0f\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                PPMstar[i].ppmRef,
                10.0 * dist);
            continue;
        }
        goodStarsPosition++;
        akkuDistError += dist * dist;

        double ppmVmag = PPMstar[i].vmag;
        double cdVmag = CDstar[cdIndex].vmag;
        if (ppmVmag < 0.00001 || cdVmag > 29.9) continue; // se omiten aquellas estrellas con Vmag=0 o variables
        ppmVmag = compVmagToCDmag(CDstar[cdIndex].declRef, ppmVmag);
        double delta = fabs(ppmVmag - cdVmag);
        if (delta >= MAX_MAGNITUDE) {
            // diferencia en magnitud V y visual supera umbral
            magDiffError++;
            indexError++;
            printf("%d) CD %d°%d reports mag=%.1f but it should be mag=%.1f from PPM %d%s: Delta = %.1f.\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                cdVmag,
                ppmVmag,
                PPMstar[i].ppmRef,
                isProblematic ? " (PROBLEM)" : "",
                delta);
            writeRegister(cdIndex);
            if (!revise(i)) totalErrorsMinusDoubles++;
            fprintf(magStream, "%d,%d,%d,PPM %d,%.0f\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                PPMstar[i].ppmRef,
                10.0 * delta);
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
