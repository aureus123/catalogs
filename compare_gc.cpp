
/*
 * COMPARE_GC - Compara registros de CD contra GC (y otros catálogos)
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "read_old.h"
#include "misc.h"

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
	readGC(true, 0, 0, 0, 0, 0, 0);
	struct GCstar_struct *GCstar = getGCStruct();
	int GCstars = getGCStars();

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
            printf("%d) GC %d is alone within MAG = %.1f; similar CD %d°%d separated in %.1f arcsec although nearest CD %d°%d at %.1f arcsec.\n",
                indexError,
                GCstar[i].gcRef,
                gcVmag,
                CDstar[cdIndexWithinMag].declRef,
                CDstar[cdIndexWithinMag].numRef,
                distWithinMag,
                CDstar[cdIndex].declRef,
                CDstar[cdIndex].numRef,
                dist);
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
