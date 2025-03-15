
/*
 * CROSS_GC - Compara registros de PPM y CD contra GC
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include "read_ppm.h"
#include "read_dm.h"
#include "read_gc.h"
#include "trig.h"
#include "misc.h"

#define CURATED true // true if curated CD catalog should be used
#define PRINT_WARNINGS false // true if print warnings about stars without CD star near them

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    char catName[20], cdName[20], ppmName[20];

    printf("CROSS_GC - Compare GC and PPM/CD catalogs.\n");
    printf("Made in 2025 by Daniel Severin.\n");

    /* leemos catalogo CD */
    readDM(CURATED ? "cat/cd_curated.txt" : "cat/cd.txt");
    struct DMstar_struct *CDstar = getDMStruct();
    int CDstars = getDMStars();

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, 1875.0);
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo GC */
	readGC();
	struct GCstar_struct *GCstar = getGCStruct();
	int GCstars = getGCStars();

	printf("\n***************************************\n");
    printf("Perform comparison between GC and PPM/CD...\n");

	FILE *crossCDStream = openCrossFile("results/cross_gc_cd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross_gc_ppm.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countCD = 0;
    int errors = 0;
    for (int gcIndex = 0; gcIndex < GCstars; gcIndex++) {
        double x = GCstar[gcIndex].x;
        double y = GCstar[gcIndex].y;
        double z = GCstar[gcIndex].z;
        float vmag = GCstar[gcIndex].vmag;
        int gcRef = GCstar[gcIndex].gcRef;

		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
        bool ppmFound = false;
		/* busca la PPM mas cercana y genera el cruzamiento */
		findPPMByCoordinates(x, y, z, &ppmIndex, &minDistance);
		if (minDistance < MAX_DIST_PPM) {
			float ppmVmag = PPMstar[ppmIndex].vmag;
			if (vmag > __FLT_EPSILON__ && ppmVmag > __FLT_EPSILON__ && !GCstar[gcIndex].dpl) {
				// Note: no fit is performed to convert scales of magnitudes
				float delta = fabs(vmag - ppmVmag);
				if (delta < MAX_MAGNITUDE) {
                    akkuDeltaError += delta * delta;
                    countDelta++;
                } else {
					printf("%d) Warning: GC has a near PPM %d star (vmag = %.1f) with a different magnitude (delta = %.1f).\n",
                        ++errors,
						PPMstar[ppmIndex].ppmRef,
						ppmVmag,
						delta);	
                    writeRegisterGC(gcIndex);
				}
			}

            akkuDistError += minDistance * minDistance;
            countDist++;
            ppmFound = true;

			snprintf(catName, 20, "GC %d", gcRef);
			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		} else {
            if (PRINT_WARNINGS) {
                printf("Warning: GC has no PPM star near it.\n");
                writeRegisterGC(gcIndex);
            }
		}
		
		int cdIndex = -1;
        minDistance = HUGE_NUMBER;
        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
		if (GCstar[gcIndex].Decld >= 22) {
		    for (int i = 0; i < CDstars; i++) {
			    double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
			    if (minDistance > dist) {
				    cdIndex = i;
				    minDistance = dist;
			    }
            }
			if (minDistance < MAX_DIST_CD) {
			    float cdVmag = CDstar[cdIndex].vmag;
			    if (vmag > __FLT_EPSILON__ && cdVmag < 29.9 && !GCstar[gcIndex].dpl) {
				    float delta = fabs(vmag - cdVmag);
				    if (delta >= MAX_MAGNITUDE) {
					    printf("%d) Warning: GC has a near CD with a different magnitude (delta = %.1f).\n",
                            ++errors,
                            delta);
                        writeRegisterGC(gcIndex);
					    writeRegister(cdIndex, false);
				    }
			    }

                countCD++;
                cdFound = true;

			    snprintf(catName, 20, "GC %d", gcRef);
			    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
			    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
			} else {
				if (PRINT_WARNINGS) {
					printf("Warning: GC has no CD star near it.\n");
                    writeRegisterGC(gcIndex);
				}
			}
		}

        if (!ppmFound && !cdFound) {
            printf("%d) Warning: GC is ALONE (no PPM or CD star near it).\n", ++errors);
            writeRegisterGC(gcIndex);
            if (GCstar[gcIndex].cum) {
                printf("  Possible cause: cumulus.\n");
            }
            if (GCstar[gcIndex].neb) {
                printf("  Possible cause: nebula.\n");
            }
            if (GCstar[gcIndex].RAs == 0) {
                printf("  Possible cause: lack of RA (s).\n");
            }
            if (GCstar[gcIndex].Decls == 0) {
                printf("  Possible cause: lack of Decl (s).\n");
            }
            if (GCstar[gcIndex].vmag > 9.9) {
                printf("  Possible cause: dim star.\n");
            }
            if (GCstar[gcIndex].Decld < 22) {
                printf("  Note: no CD coverage for stars below 22°.\n");
            }
            if (GCstar[gcIndex].Decld > 61) {
                printf("  Note: poor CD coverage for stars above 61°.\n");
            }
        }
	}
	fclose(crossPPMStream);
	fclose(crossCDStream);

    printf("Stars from GC identified with PPM = %d, and CD = %d\n", countDist, countCD);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)countDelta),
        countDelta);
    printf("Errors logged = %d\n", errors);
    return 0;
}
