
/*
 * COMPARE_CPD - Compara registros de CD contra CPD
 * Copyright (C) 2024 Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "read_cpd.h"
#include "trig.h"
#include "misc.h"

#define MAX_DISTANCE 180.0   // 3 minutos de arco
#define MAX_DISTANCE_DROP 3600.0 // 1 grado de arco
#define CROSS_CATALOG true   // true = 4011, false = 4005

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("COMPARE_CPD - Compare CD and CPD catalogs.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    /* leemos catalogo CD */
    readDM("cat/cd.txt");
    struct DMstar_struct *CDstar = getDMStruct();

    /* leemos catalogo CPD */
    readCPD(CROSS_CATALOG);
    int CPDstars = getCPDStars();
    struct CPDstar_struct *CPDstar = getCPDStruct();

    /* revisamos la identificación cruzada y generamos una planilla */
    FILE *posStream, *magStream;

    posStream = fopen("results/table_pos_cpd.csv", "wt");
    if (posStream == NULL) {
        perror("Cannot write in table_pos_cpd.csv");
        exit(1);
    }
    fprintf(posStream, "index,decl,num,ref,dist10\n");

    int dropDistError = 0;
    int indexError = 0;
    int goodStarsPosition = 0;
    double akkuDistError = 0.0;
    for (int i = 0; i < CPDstars; i++) {
        if (CPDstar[i].discard == true) continue;
        int cdIndex = CPDstar[i].dmIndex;
        float dist = CPDstar[i].dist;
        if (dist > MAX_DISTANCE_DROP) {
            // posiciones demasiado separadas, posiblemente mala identificacion
            dropDistError++;
            printf("*) CD %d°%d too separated from CPD %d°%d in %.1f arcsec.\n",
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                CPDstar[i].declRef,
                CPDstar[i].numRef,
                dist);
            writeRegister(cdIndex);
            continue;
        }
        if (dist > MAX_DISTANCE) {
            // posiciones muy separadas, supera umbral
            indexError++;
            printf("%d) CD %d°%d separated from CPD %d°%d in %.1f arcsec.\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                CPDstar[i].declRef,
                CPDstar[i].numRef,
                dist);
            if (CPDstar[i].declRef2 != -1) {
                printf("    (also identified as CD %d°%d)\n", CPDstar[i].declRef2, CPDstar[i].numRef2);
            }
            writeRegister(cdIndex);
            fprintf(posStream, "%d,%d,%d,CPD %d°%d,%.0f\n",
                indexError,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                CPDstar[i].declRef,
                CPDstar[i].numRef,
                10.0 * dist);
            continue;
        }
        goodStarsPosition++;
        akkuDistError += dist * dist;
    }
    fclose(posStream);
    printf("Total errors: %d, Drops = %d\n", indexError, dropDistError);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)goodStarsPosition),
        goodStarsPosition);
    return 0;
}
