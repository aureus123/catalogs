
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
//#define MAGNITUDE_METHOD

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
        printf("        cd_vol1_curated.txt = curated version of cd_vol1.txt\n");
        printf("        I88.txt = 1982 CD catalog version (has some errors)\n");
        exit(-1);
    }
    bool allSky = true;
    if (strcmp(argv[1], "cd.txt") && strcmp(argv[1], "1114.txt")) {
        if (strcmp(argv[1], "cd_vol1.txt") && strcmp(argv[1], "cd_vol1_curated.txt") && strcmp(argv[1], "I88.txt")) {
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
    readPPM(true, allSky, true, false, 1875.0);
    int PPMstars = getPPMStars();
    struct PPMstar_struct *PPMstar = getPPMStruct();

    /* revisamos la identificación cruzada y generamos planillas */
    FILE *posStream, *magStream, *crossPPMStreamV1, *crossPPMStreamV2;
    FILE *crossPPMStreamV3, *crossPPMStreamV4, *crossPPMStreamV5;
    posStream = openPositionFile("results/table_pos_ppm.csv");
    magStream = openMagnitudeFile("results/table_mag_ppm.csv");
    if (allSky) {
        /* si cubrimos todo el cielo, generamos planillas para cada volumen */
        crossPPMStreamV1 = openCrossFile("results/cross/cross_cd_vol1_ppm.csv");
        crossPPMStreamV2 = openCrossFile("results/cross/cross_cd_vol2_ppm.csv");
        crossPPMStreamV3 = openCrossFile("results/cross/cross_cd_vol3_ppm.csv");
        crossPPMStreamV4 = openCrossFile("results/cross/cross_cd_vol4_ppm.csv");
        crossPPMStreamV5 = openCrossFile("results/cross/cross_cd_vol5_ppm.csv");
    }

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

        char ppmName[20];
        snprintf(ppmName, 20, "PPM %d", PPMstar[i].ppmRef);

        int cdIndex = PPMstar[i].dmIndex;
        char *cdString = PPMstar[i].dmString;
        float dist = PPMstar[i].dist;
        if (dist > MAX_DISTANCE) {
            // posiciones muy separadas, supera umbral
            maxDistError++;
            indexError++;
            printf("%d) %s separated from %s%s in %.1f arcsec.\n",
                indexError,
                cdString,
                ppmName,
                isProblematic ? " (PROBLEM)" : "",
                dist);
            writeRegister(cdIndex, true);
            if (!revise(i)) totalErrorsMinusDoubles++;
            writePositionEntry(posStream,
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                ppmName,
                dist);
            continue;
        }
        goodStarsPosition++;
        akkuDistError += dist * dist;

        double ppmVmag = PPMstar[i].vmag;
        double cdVmag = CDstar[cdIndex].vmag;
        if (ppmVmag < 0.00001 || cdVmag > 29.9) continue; // se omiten aquellas estrellas con Vmag=0 o variables
#ifdef MAGNITUDE_METHOD        
        double convertedCDmag = compCDmagToVmag(CDstar[cdIndex].declRef, cdVmag);
        double delta = fabs(ppmVmag - convertedCDmag);
#else
        double convertedPPMVmag = compVmagToCDmag(CDstar[cdIndex].declRef, ppmVmag);
        double delta = fabs(cdVmag - convertedPPMVmag);
#endif        
        if (delta >= MAX_MAGNITUDE) {
            // diferencia en magnitud V y visual supera umbral
            magDiffError++;
            indexError++;
            printf("%d) %s reports mag=%.1f but %s%s has Vmag=%.1f: Delta = %.1f.\n",
                indexError,
                cdString,
                cdVmag,
                ppmName,
                isProblematic ? " (PROBLEM)" : "",
                ppmVmag,
                delta);
            writeRegister(cdIndex, false);
            if (!revise(i)) totalErrorsMinusDoubles++;
            writeMagnitudeEntry(magStream,
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                ppmName,
                delta);
            continue;
        }

        if (allSky) {
            // Las estrellas sin errores se cruzan con PPM y se almacenan según el volumen
            int decl = -CDstar[cdIndex].declRef;
            FILE *targetStream = NULL;
            if (decl >= 22 && decl <= 31) {
                targetStream = crossPPMStreamV1;
            } else if (decl >= 32 && decl <= 41) {
                targetStream = crossPPMStreamV2;
            } else if (decl >= 42 && decl <= 51) {
                targetStream = crossPPMStreamV3;
            } else if (decl >= 52 && decl <= 61) {
                targetStream = crossPPMStreamV4;
            } else if (decl >= 62 && decl <= 89) {
                targetStream = crossPPMStreamV5;
            }
            if (targetStream != NULL) {
                writeCrossEntry(targetStream, cdString, ppmName, cdVmag, dist);
            }
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
