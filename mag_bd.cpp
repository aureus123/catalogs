
/*
 * MAG_BD - Genera datos para ajustar magnitudes visuales de BD
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
    printf("MAG_BD - Compare registers from BD and PPM for quadratic fit.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    /* leemos catalogo BD */
    readDM("cat/bd.txt");
    int BDstars = getDMStars();
    struct DMstar_struct *BDstar = getDMStruct();

    /* leemos catalogo PPM */
    readPPM(true, true, false, false, 1855.0);
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
        /* descartamos problemáticas, con Vmag=0 (hay muchas fotográficas en esta franja) o muy separadas */
        if (PPMstar[i].problem == 1) continue;
        double ppmVmag = PPMstar[i].vmag;
        if (ppmVmag < 0.00001) continue;
        int bdIndex = PPMstar[i].dmIndex;
        float dist = 3600.0 * calcAngularDistance(PPMstar[i].x, PPMstar[i].y, PPMstar[i].z, BDstar[bdIndex].x, BDstar[bdIndex].y, BDstar[bdIndex].z);
        if (dist > MAX_DISTANCE) continue;

        /* generamos la fila con las magnitudes */
        double bdVmag = BDstar[bdIndex].vmag;
        fprintf(stream, "1,%.2f,%.2f\n", ppmVmag, bdVmag);
        filas++;
    }
    fclose(stream);
    printf("Rows generated = %d\n", filas);
    return 0;
}
