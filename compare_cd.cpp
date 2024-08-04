
/*
 * COMPARE_CD - Compara registros de CD contra si mismos
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "trig.h"
#include "misc.h"


static struct DMstar_struct CDstar2[MAXDMSTAR];

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("COMPARE_CD - Compare CD against itself.\n");
    printf("Made in 2024 by Daniel Severin.\n");

    if (argc < 3) {
        printf("Usage: compare_ppm file1 file2\n");
        printf("    where files can be:\n");
        printf("        cd.txt = Current CD catalog at Vizier\n");
        printf("        1114.txt = NASA-ADC CD catalog (has some errors)\n");
        printf("        cd_vol1.txt = same as cd.txt but only 1st. Volume (Resultados XVI)\n");
        printf("        I88.txt = 1982 CD catalog version (has some errors)\n");
        printf("You should compare cd.txt against 1114.txt, or cd_vol1.txt against I88.txt\n");
        exit(-1);
    }

    /* leemos catalogo CD 2 */
    char buffer[64];
    snprintf(buffer, 64, "cat/%s", argv[2]);
    readDM(buffer);
    int stars = getDMStars();
    struct DMstar_struct *CDstar1 = getDMStruct();
    for (int i = 0; i < stars; i++) {
        memcpy(&CDstar2[i], &CDstar1[i], sizeof(DMstar_struct));
    }

    /* leemos catalogo CD 1 */
    snprintf(buffer, 64, "cat/%s", argv[1]);
    readDM(buffer);

    /* comparamos */
    int diffRA = 0;
    int diffDE = 0;
    int diffMag = 0;
    int sameStar = 0;
    for (int i2 = 0; i2 < stars; i2++) {
        if (CDstar2[i2].supplRef != ' ') continue;
        int declRef = CDstar2[i2].declRef;
        int numRef = CDstar2[i2].numRef;
        int i1 = getDMindex(true, declRef, numRef);
        if (fabs(CDstar1[i1].rah - CDstar2[i2].rah) > EPS ||
            fabs(CDstar1[i1].ramin - CDstar2[i2].ramin) > EPS ||
            fabs(CDstar1[i1].raseg - CDstar2[i2].raseg) > EPS) {
            printf("Difference in RA for CD %d°%d\n", declRef, numRef);
            diffRA++;
            continue;
        }
        if (fabs(CDstar1[i1].decldeg - CDstar2[i2].decldeg) > EPS ||
            fabs(CDstar1[i1].declmin - CDstar2[i2].declmin)) {
            printf("Difference in DE for CD %d°%d\n", declRef, numRef);
            diffDE++;
            continue;
        }
        if (CDstar1[i1].vmag > 19.0 || CDstar2[i2].vmag > 19.0) continue;
        if (fabs(CDstar1[i1].vmag - CDstar2[i2].vmag) > EPS) {
            printf("Difference in Mag for CD %d°%d\n", declRef, numRef);
            diffMag++;
            continue;
        }
        sameStar++;
    }
    printf("Errors in RA: %d,  DE: %d,  Mag: %d\n", diffRA, diffDE, diffMag);
    printf("Same Stars: %d\n", sameStar);
    return 0;
}