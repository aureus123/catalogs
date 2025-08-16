
/*
 * CROSS_GC - Compara registros de PPM, CD y CPD contra GC
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include "read_ppm.h"
#include "read_cpd.h"
#include "read_dm.h"
#include "read_gc.h"
#include "trig.h"
#include "misc.h"
#include "find_gsc.h"

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

    /** leemos catalogo CPD */
    readCPD(false, false);
    struct CPDstar_struct *CPDstar = getCPDStruct();

    /* leemos catalogo PPM (pero no es necesario cruzarlo con DM) */
    readPPM(false, true, false, false, 1875.0);
    sortPPM();
	struct PPMstar_struct *PPMstar = getPPMStruct();
    int PPMstars = getPPMStars();

    /* leemos catalogo GC */
	readGC();
	struct GCstar_struct *GCstar = getGCStruct();
	int GCstars = getGCStars();

	printf("\n***************************************\n");
    printf("Perform comparison between GC and PPM/CD/CPD...\n");

	FILE *crossCDStream = openCrossFile("results/cross/cross_gc_cd.csv");
	FILE *crossCPDStream = openCrossFile("results/cross/cross_gc_cpd.csv");
	FILE *crossPPMStream = openCrossFile("results/cross/cross_gc_ppm.csv");
	FILE *unidentifiedStream = openUnidentifiedFile("results/cross/gc_unidentified.csv");

    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countCD = 0;
    int countCPD = 0;
    int errors = 0;
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
                    ++errors,
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
                    ++errors,
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
        float vmag = GCstar[gcIndex].vmag;

        bool ppmFound = false;
		int ppmIndex = -1;
		double minDistance = HUGE_NUMBER;
		/* busca la PPM mas cercana y genera el cruzamiento */
		findPPMByCoordinates(x, y, z, decl, &ppmIndex, &minDistance);
        double nearestPPMDistance = minDistance;
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

			snprintf(ppmName, 20, "PPM %d", PPMstar[ppmIndex].ppmRef);
			writeCrossEntry(crossPPMStream, catName, ppmName, minDistance);
		} else {
            if (PRINT_WARNINGS) {
                printf("Warning: GC has no PPM star near it.\n");
                writeRegisterGC(gcIndex);
            }
		}

        bool cpdFound = false;
        /* busca la CPD mas cercana y genera el cruzamiento */
        if (GCstar[gcIndex].Decld >= 18) {
            int cpdIndex = -1;
            minDistance = HUGE_NUMBER;
            findCPDByCoordinates(x, y, z, decl, &cpdIndex, &minDistance);
			if (minDistance < MAX_DIST_CPD) {
                countCPD++;
                cpdFound = true;

			    snprintf(cdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
			    writeCrossEntry(crossCPDStream, catName, cdName, minDistance);
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

			    snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
			    writeCrossEntry(crossCDStream, catName, cdName, minDistance);
			} else {
				if (PRINT_WARNINGS) {
					printf("Warning: GC has no CD star near it.\n");
                    writeRegisterGC(gcIndex);
				}
			}
		}

        if (!ppmFound && !cdFound && !cpdFound) {
            printf("%d) Warning: GC is ALONE (no PPM or CD or CPD star near it).\n", ++errors);
            writeRegisterGC(gcIndex);
            bool store = logCauses(catName, unidentifiedStream, x, y, z,
                GCstar[gcIndex].cum, GCstar[gcIndex].neb, GCstar[gcIndex].vmag,
                GCstar[gcIndex].RAs, decl, GCstar[gcIndex].Decls,
                PPMstar[ppmIndex].ppmRef, nearestPPMDistance);
            if (store) {
                if (!findGSCStar(ra, decl, 1875.0, MAX_DIST_GSC)) {
                    printf("  Warning: no nearby GSC star found.\n");
                }
            }
        }
	}
    fclose(unidentifiedStream);
	fclose(crossPPMStream);
	fclose(crossCPDStream);
	fclose(crossCDStream);

    printf("Stars from GC identified with PPM = %d, CD = %d and CPD = %d\n", countDist, countCD, countCPD);
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(akkuDistError / (double)countDist),
        countDist);
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(akkuDeltaError / (double)countDelta),
        countDelta);
    printf("Errors logged = %d\n", errors);
    return 0;
}
