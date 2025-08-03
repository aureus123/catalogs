/*
 * CROSS_TXT - Counts stars brighter than 8th magnitude in PPM catalog and matches with Coord.txt
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "trig.h"
#include "read_ppm.h"

#define MAX_CUSTOM_STARS 10000
#define THRESHOLD_PPM 15.0
#define THRESHOLD_MAG 1.5

/* Structure to store custom star data */
struct CustomStars {
    int index;
    int ppmIndex;
    double RA;
    double Decl;
    double imag;  /* instrumental magnitude */
};

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("CROSS_TXT - Cross custom catalog with PPM catalog.\n");
    printf("Made in 2025 by Daniel Severin.\n\n");

    /* leemos catalogo PPM */
    printf("Reading PPM catalog...\n");
    readPPM(false, true, false, false, 2000.0);
    sortPPM();
    int PPMstars = getPPMStars();
    struct PPMstar_struct *PPMstar = getPPMStruct();

    /* ahora leemos Coord.txt y Flux.txt y buscamos coincidencias con PPM */
    printf("\nReading Coord.txt and Flux.txt and finding PPM matches...\n");
    
    FILE *coordFile = fopen("Coord.txt", "rt");
    if (coordFile == NULL) {
        printf("Error: Cannot open Coord.txt file.\n");
        return 1;
    }

    FILE *fluxFile = fopen("Flux.txt", "rt");
    if (fluxFile == NULL) {
        printf("Error: Cannot open Flux.txt file.\n");
        fclose(coordFile);
        return 1;
    }

    char coordLine[256], fluxLine[256];
    int starIndex, fluxIndex;
    double RA, Decl, flux, background;
    int matches = 0;
    
    /* Array to store custom stars data */
    struct CustomStars customStars[MAX_CUSTOM_STARS];
    int customStarCount = 0;
    
    /* Skip header lines */
    if (fgets(coordLine, sizeof(coordLine), coordFile) == NULL) {
        printf("Error: Cannot read header from Coord.txt.\n");
        fclose(coordFile);
        fclose(fluxFile);
        return 1;
    }
    
    if (fgets(fluxLine, sizeof(fluxLine), fluxFile) == NULL) {
        printf("Error: Cannot read header from Flux.txt.\n");
        fclose(coordFile);
        fclose(fluxFile);
        return 1;
    }
    
    while (fgets(coordLine, sizeof(coordLine), coordFile) != NULL && 
           fgets(fluxLine, sizeof(fluxLine), fluxFile) != NULL) {
        if (sscanf(coordLine, "%d %lf %lf", &starIndex, &RA, &Decl) == 3 &&
            sscanf(fluxLine, "%d %lf %lf", &fluxIndex, &flux, &background) == 3) {
            
            /* Verify that both files have the same star index */
            if (starIndex != fluxIndex) {
                printf("Warning: Index mismatch at line %d (Coord: %d, Flux: %d)\n", 
                       starIndex, starIndex, fluxIndex);
                continue;
            }

            //printf("Star %d: RA = %f, Decl = %f, flux = %f, background = %f\n", 
            //       starIndex, RA, Decl, flux, background);
            
            double x, y, z;
            sph2rec(RA, Decl, &x, &y, &z);
            
            /* Find nearest PPM star within 15 arc seconds */
            int ppmIndex = -1;
            double minDistance = HUGE_NUMBER;
            findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);

            /* Omit brightest stars (or without Vmag) and those already matched */
            if (ppmIndex != -1
                    && minDistance < THRESHOLD_PPM
                    && !PPMstar[ppmIndex].discard
                    && PPMstar[ppmIndex].vmag > 3.0) {
                double imag = -2.5 * log10(flux);
                
                printf("Star %d with flux=%.2f is matched with PPM %d (dist = %.1f arcsec.)\n", 
                       starIndex, flux, PPMstar[ppmIndex].ppmRef, minDistance);
                
                /* Save match in custom stars structure */
                if (customStarCount < MAX_CUSTOM_STARS) {
                    customStars[customStarCount].index = starIndex;
                    customStars[customStarCount].ppmIndex = ppmIndex;
                    customStars[customStarCount].RA = RA;
                    customStars[customStarCount].Decl = Decl;
                    customStars[customStarCount].imag = imag;
                    customStarCount++;
                }
                
                PPMstar[ppmIndex].discard = true;
                matches++;
                continue;
            }
        }
    }
    
    fclose(coordFile);
    fclose(fluxFile);
    printf("\nTotal matches found: %d\n", matches);

    /* Perform linear regression fit between instrumental and PPM V magnitudes */
    if (customStarCount > 5) {
        printf("\nPerforming linear regression fit...\n");
        
        /* Calculate means for least squares fit */
        double sum_imag = 0.0, sum_vmag = 0.0, sum_imag_vmag = 0.0, sum_imag_sq = 0.0;
        
        for (int i = 0; i < customStarCount; i++) {
            int ppmIdx = customStars[i].ppmIndex;
            sum_imag += customStars[i].imag;
            sum_vmag += PPMstar[ppmIdx].vmag;
            sum_imag_vmag += customStars[i].imag * PPMstar[ppmIdx].vmag;
            sum_imag_sq += customStars[i].imag * customStars[i].imag;
        }
        
        double mean_imag = sum_imag / customStarCount;
        double mean_vmag = sum_vmag / customStarCount;
        
        /* Calculate slope and intercept */
        double denominator = sum_imag_sq - customStarCount * mean_imag * mean_imag;
        double factor = (sum_imag_vmag - customStarCount * mean_imag * mean_vmag) / denominator;
        double zeroPoint = mean_vmag - factor * mean_imag;
        
        printf("Linear fit: Vmag = %.3f + %.3f * imag\n", zeroPoint, factor);
        
        /* Calculate calibrated magnitudes, RMSE, and MAPE */
        double sum_sq_error = 0.0;
        double sum_abs_percentage_error = 0.0;
        printf("\nCalibrated magnitudes:\n");
         
        for (int i = 0; i < customStarCount; i++) {
            int ppmIndex = customStars[i].ppmIndex;
            double ppm_vmag = PPMstar[ppmIndex].vmag;
            double calibrated_vmag = zeroPoint + factor * customStars[i].imag;
            double error = fabs(calibrated_vmag - ppm_vmag);
            sum_sq_error += error * error;
            sum_abs_percentage_error += error / ppm_vmag * 100.0;

            if (error > THRESHOLD_MAG) {
                printf("Star %d has Vmag=%.1f while PPM %d has Vmag=%.1f, difference=%.1f\n", 
                        customStars[i].index, calibrated_vmag, PPMstar[ppmIndex].ppmRef, ppm_vmag, error);
            }
        }
        
        double rmse = sqrt(sum_sq_error / customStarCount);
        double mape = sum_abs_percentage_error / customStarCount;
        printf("\nRMSE: %.3f magnitudes\n", rmse);
        printf("MAPE: %.2f%%\n", mape);
    }

    return 0;
} 