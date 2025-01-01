
/*
 * MAG_CD - Genera datos para ajustar magnitudes visuales de CD
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

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("MAG_CD - Compare registers from CD and PPM for quadratic fit.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    /* leemos catalogo CD */
    readDM("cat/cd.txt");
    int CDstars = getDMStars();
    struct DMstar_struct *CDstar = getDMStruct();

    /* leemos catalogo PPM */
    readPPM(true, true, true, false, 1875.0);
    int PPMstars = getPPMStars();
    struct PPMstar_struct *PPMstar = getPPMStruct();

    /* revisamos la identificación cruzada y generamos la matriz */
    FILE *stream;

    stream = fopen("results/matrices.csv", "wt");
    if (stream == NULL) {
        perror("Cannot write matrices.csv");
        exit(1);
    }

    int filas = 0;
    for (int i = 0; i < PPMstars; i++) {
        /* tomamos 1 de cada 3 estrellas */
        if (i%3 > 0) continue;
        
        /* también descartamos problemáticas, con Vmag=0 o muy separadas */
        if (PPMstar[i].problem == 1) continue;
        double ppmVmag = PPMstar[i].vmag;
        if (ppmVmag < 0.00001) continue;
        int cdIndex = PPMstar[i].dmIndex;
        float dist = 3600.0 * calcAngularDistance(PPMstar[i].x, PPMstar[i].y, PPMstar[i].z, CDstar[cdIndex].x, CDstar[cdIndex].y, CDstar[cdIndex].z);
        if (dist > MAX_DISTANCE) continue;

        /* obtenemos el tomo */
        int tomo = 0;
        int declRef = CDstar[cdIndex].declRef;
        if (declRef <= -22 && declRef >= -31) tomo = 1;
        if (declRef <= -32 && declRef >= -41) tomo = 2;
        if (declRef <= -42 && declRef >= -51) tomo = 3;
        if (declRef <= -52 && declRef >= -61) tomo = 4;
        if (declRef <= -62) tomo = 5;

        /* generamos la fila con las magnitudes */
        double cdVmag = CDstar[cdIndex].vmag;
        fprintf(stream, "%d,%.2f,%.2f\n", tomo, ppmVmag, cdVmag);
        filas++;
    }
    fclose(stream);
    printf("Rows generated = %d\n", filas);
    return 0;
}
