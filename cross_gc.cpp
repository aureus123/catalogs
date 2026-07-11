
/*
 * CROSS_GC - Compara registros de PPM, CD y CPD contra GC
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "read_ppm.h"
#include "read_cpd.h"
#include "read_dm.h"
#include "read_gc.h"
#include "trig.h"
#include "misc.h"
#include "find_gsc.h"
#include "cross_utils.h"

#define CURATED true // true if curated CD catalog should be used
#define PRINT_WARNINGS false // true if print warnings about stars without CD star near them

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    char catName[20], cdName[20];

    printf("CROSS_GC - Compare GC and PPM/CD catalogs.\n");
    printf("Made in 2025 by Daniel Severin.\n");

    /* leemos catalogo CD */
    readDM(CURATED ? "cat/cd_curated.txt" : "cat/cd.txt");
    struct DMstar_struct *CDstar = getDMStruct();

    /** leemos catalogo CPD */
    readCPD(false, false);
    struct CPDstar_struct *CPDstar = getCPDStruct();

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    struct PPMstar_struct *PPMstar = preparePPM(1875.0, false);

    /* leemos catalogo GC */
	readGC();
	struct GCstar_struct *GCstar = getGCStruct();
	int GCstars = getGCStars();

	printf("\n***************************************\n");
    printf("Perform comparison between GC and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_gc_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_gc_cpd.csv");
    FILE *crossPPMStream, *crossSAOStream, *crossHDStream;
    openCrossSet("gc", &crossPPMStream, &crossSAOStream, &crossHDStream);
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/gc_unidentified.csv");
    FILE *catalogStream = openCatalogFile("likelihood/cat1875/gc.csv");

    CrossStats stats;
    for (int gcIndex = 0; gcIndex < GCstars; gcIndex++) {
        double ra = GCstar[gcIndex].RA1875;
        double decl = GCstar[gcIndex].Decl1875;
        int gcRef = GCstar[gcIndex].gcRef;
        snprintf(catName, 20, "GC %d", gcRef);

        /* revisa precesion en RA y Decl (usamos constantes de Struve) */
        double preRA = GCstar[gcIndex].preRA;
        if (fabs(preRA) > EPS) {
            double realPreRA = 3.072245 + 1.33695 * dsin(ra) * dtan(decl);
            double diff = fabs(preRA - realPreRA);
            if (diff > 0.0099) {
                printf("%d) Warning: GC %d reports %.3f on precession RA but it should be %.3f (diff=%.3f).\n",
                    ++stats.errors,
                    gcRef,
                    preRA,
                    realPreRA,
                    diff);
                writeRegisterGC(gcIndex);
            }
        }

        double preDecl = GCstar[gcIndex].preDecl;
        if (fabs(preDecl) > EPS) {
            double realPreDecl = 20.05425 * dcos(ra);
            double diff = fabs(preDecl - realPreDecl);
            if (diff > 0.099) {
                printf("%d) Warning: GC %d reports %.3f on precession DECL but it should be %.3f (diff=%.3f).\n",
                    ++stats.errors,
                    gcRef,
                    preDecl,
                    realPreDecl,
                    diff);
                writeRegisterGC(gcIndex);
            }
        }

        double x = GCstar[gcIndex].x;
        double y = GCstar[gcIndex].y;
        double z = GCstar[gcIndex].z;
        float gcVmag = GCstar[gcIndex].vmag;

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		findPPMByCoordinates(x, y, z, decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
		if (minDistance < MAX_DIST_PPM) {
			float ppmVmag = PPMstar[ppmIndex].vmag;
			if (gcVmag > __FLT_EPSILON__ && fabs(ppmVmag) > __FLT_EPSILON__ && !GCstar[gcIndex].dpl) {
                float vmag = compGCmagToVmag(gcVmag); // perform magnitude transformation
				float delta = fabs(vmag - ppmVmag);
				if (delta < MAX_MAGNITUDE) {
                    stats.akkuDeltaError += delta * delta;
                    stats.countDelta++;
                } else {
					printf("%d) Warning: GC has a near PPM %d star (vmag = %.1f) with a different magnitude (delta = %.1f).\n",
                        ++stats.errors,
						PPMstar[ppmIndex].ppmRef,
						ppmVmag,
						delta);
                    writeRegisterGC(gcIndex);
				}
			}

            stats.akkuDistError += minDistance * minDistance;
            stats.countDist++;
            ppmFound = true;

            writePPMCrossEntry(crossPPMStream, crossSAOStream, crossHDStream, catName, &PPMstar[ppmIndex], gcVmag, minDistance);
		} else {
            if (PRINT_WARNINGS) {
                printf("Warning: GC has no PPM star near it.\n");
                writeRegisterGC(gcIndex);
            }
		}

        /* si no encuentra una PPM cercana, prueba con GSC */
        bool gscFound = tryGSC(ppmFound, ra, decl, 1875.0, crossPPMStream, catName, gcVmag, &stats);

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        if (GCstar[gcIndex].Decld >= 18) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            findCPDByCoordinates(x, y, z, decl, &cpdIndex, &minDistance);
			if (minDistance < MAX_DIST_CPD) {
                stats.countCPD++;
                cpdFound = true;

			    snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
			    writeCrossEntry(crossCPDStream, catName, cdName, gcVmag, minDistance);
			} else {
				if (PRINT_WARNINGS) {
					printf("Warning: GC has no CPD star near it.\n");
                    writeRegisterGC(gcIndex);
				}
			}
        }

        bool cdFound = false;
		/* busca la CD mas cercana y genera el cruzamiento */
		if (GCstar[gcIndex].Decld >= 22) {
            int cdIndex = -1;
            minDistance = HUGE_NUMBER;
            findDMByCoordinates(x, y, z, decl, &cdIndex, &minDistance);
			if (minDistance < MAX_DIST_CD) {
			    float cdVmag = CDstar[cdIndex].vmag;
			    if (gcVmag > __FLT_EPSILON__ && cdVmag < 29.9 && !GCstar[gcIndex].dpl) {
				    float delta = fabs(gcVmag - cdVmag);
				    if (delta >= MAX_MAGNITUDE) {
					    printf("%d) Warning: GC has a near CD with a different magnitude (delta = %.1f).\n",
                            ++stats.errors,
                            delta);
                        writeRegisterGC(gcIndex);
					    writeRegister(cdIndex, false);
				    }
			    }

                stats.countCD++;
                cdFound = true;

			    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
			    writeCrossEntry(crossCDStream, catName, cdName, gcVmag, minDistance);
			} else {
				if (PRINT_WARNINGS) {
					printf("Warning: GC has no CD star near it.\n");
                    writeRegisterGC(gcIndex);
				}
			}
		}

        if (!ppmFound && !cdFound && !cpdFound) {
            if (gscFound) {
                writeUnidentified(unidentifiedStream, catName, x, y, z);
            } else {
                printf("%d) Warning: GC %d is ALONE (no PPM / CD / CPD / GSC star near it).\n", ++stats.errors, gcRef);
                writeRegisterGC(gcIndex);
                logCauses(catName, true,
                    GCstar[gcIndex].cum, GCstar[gcIndex].neb,
                    GCstar[gcIndex].RAs, decl, GCstar[gcIndex].Decls,
                    PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
            }
        }
        writeCatalogFile(catalogStream, catName, x, y, z, gcVmag);
	}
    fclose(unidentifiedStream);
    fclose(catalogStream);
    closeCrossSet(crossPPMStream, crossSAOStream, crossHDStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Stars from GC identified with PPM = %d, GSC-PPM = %d, CD = %d and CPD = %d\n",
        stats.countDist, stats.countGSC, stats.countCD, stats.countCPD);
    printRSMEDist(&stats);
    printRSMEMag(&stats);
    printf("Errors logged = %d\n", stats.errors);

    // Also, generate a file with all GC double stars
    double *gcX   = (double*) malloc(GCstars * sizeof(double));
    double *gcY   = (double*) malloc(GCstars * sizeof(double));
    double *gcZ   = (double*) malloc(GCstars * sizeof(double));
    double *gcMag = (double*) malloc(GCstars * sizeof(double));
    int    *gcRef = (int*)    malloc(GCstars * sizeof(int));
    for (int i = 0; i < GCstars; i++) {
        gcX[i]   = GCstar[i].x;
        gcY[i]   = GCstar[i].y;
        gcZ[i]   = GCstar[i].z;
        gcMag[i] = GCstar[i].vmag;
        gcRef[i] = GCstar[i].gcRef;
    }
    makeDoubles(GCstars, gcRef, gcX, gcY, gcZ, gcMag, "GC", "results/doubles/gc.csv");
    free(gcX);
    free(gcY);
    free(gcZ);
    free(gcMag);
    free(gcRef);
    return 0;
}
