/*
 * CROSS_TXT - Cross-identifies custom catalog against PPM catalog
 * Made in 2025 by Daniel E. Severin (partially created by Cursor AI)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "trig.h"
#include "read_ppm.h"

#define STRING_SIZE 14
#define MAX_CUSTOM_STARS 200
#define THRESHOLD_PPM 15.0
#define THRESHOLD_MAG 1.5
#define ZEROPOINT_IMAG 17.37
#define CROSS_OLD_CATALOGS true

/* Structure to store custom star data (histogram).
 * Note: instrumental magnitudes are given by indices,
 * e.g. customStars[112] gives consolidated info about stars of 11.2 imag. */
struct CustomStars {
    double ppmVmag;
    int count;
};

/*
 * readCrossFile - lee archivos de identificaciones cruzadas en formato CSV
 */
void readCrossFile(
        const char *ppm_file, struct PPMstar_struct *PPMstar, int PPMstars) {
    char buffer[1024];
    char targetRef[STRING_SIZE];
    int ppmRef, cdDeclRef, cdNumRef, cpdDeclRef, cpdNumRef;
    float vmag, minDistance;

    /* read cross file between PPM and target */
    printf("Reading cross file %s... ", ppm_file);
    FILE *stream = fopen(ppm_file, "rt");
    if (stream == NULL) {
        snprintf(buffer, 1024, "Cannot read %s", ppm_file);
        perror(buffer);
        exit(1);
    }
    bool first_line = true;
    while (fgets(buffer, 1023, stream) != NULL) {
        if (first_line) {
            // omit first line (header)
            first_line = false;
            continue;
        }
        sscanf(buffer, "%13[^,],PPM %d,%f,%f\n", targetRef, &ppmRef, &vmag, &minDistance);
        // printf("Cross: %s, PPM %d, dist = %.1f arcsec.\n", targetRef, ppmRef, minDistance);
        if (minDistance < __FLT_EPSILON__ || minDistance > THRESHOLD_PPM) {
            // omit identifications with zero distance (bug) or too far away
            continue;
        }
        for (int i = 0; i < PPMstars; i++) {
            if (PPMstar[i].ppmRef == ppmRef) {
                strncpy(PPMstar[i].dmString, targetRef, STRING_SIZE);
                break;
            }
        }
    }
    fclose(stream);
    printf("done!\n");
}

/* Join base directory and filename into out path */
static void joinPath(const char *base, const char *file, char *out, size_t outsz) {
    size_t len = strlen(base);
    if (len > 0 && (base[len - 1] == '/' || base[len - 1] == '\\')) {
        snprintf(out, outsz, "%s%s", base, file);
    } else {
        snprintf(out, outsz, "%s/%s", base, file);
    }
}

/*
 * main - comienzo de la aplicacion
 */
int main(int argc, char** argv)
{
    printf("CROSS_TXT - Cross custom catalog with PPM catalog.\n");
    printf("Made in 2025 by Daniel Severin.\n\n");

    if (argc < 2) {
        printf("Usage: %s --txt [baseDir]\n", argv[0]);
        printf("         It reads Coord.txt and Flux.txt from base folder and searches\n");
        printf("         for matches with PPM catalog. It also writes results to Cross.txt.\n");
        printf("Usage: %s --csv [csvFile]\n", argv[0]);
        printf("         It reads CSV cross-identification file and performs magnitude fits.\n");
        return 1;
    }

    const char *mode = argv[1];
    
    /* Array to store histogram of custom stars based on instrumental magnitudes. */
    struct CustomStars customStars[MAX_CUSTOM_STARS];
    for (int i = 0; i < MAX_CUSTOM_STARS; i++) {
        customStars[i].ppmVmag = 0.0;
        customStars[i].count = 0;
    }
    int customStarCount = 0;

    if (strcmp(mode, "--csv") == 0) {
        /* Mode is --csv, read CSV file and perform magnitude fits */
        if (argc < 3) {
            printf("Error: CSV file name is required for --csv mode.\n");
            return 1;
        }
        
        const char *csvFile = argv[2];
        printf("Reading CSV cross-identification file: %s\n", csvFile);
        printf("(original visual magnitudes are treated as instrumental magnitudes)\n");
        
        FILE *stream = fopen(csvFile, "rt");
        if (stream == NULL) {
            printf("Error: Cannot open %s file.\n", csvFile);
            return 1;
        }
        
        /* Load PPM catalog */
        printf("Reading PPM catalog...\n");
        readPPM(false, true, false, false, 2000.0);
        sortPPM();
        int PPMstars = getPPMStars();
        struct PPMstar_struct *PPMstar = getPPMStruct();
        
        /* Read CSV file and populate customStars histogram */
        char buffer[1024];
        char targetRef[STRING_SIZE];
        int ppmRef;
        float vmag, minDistance;
        bool first_line = true;
        int matches = 0;
        
        while (fgets(buffer, 1023, stream) != NULL) {
            if (first_line) {
                // omit first line (header)
                first_line = false;
                continue;
            }
            
            int scanned = sscanf(buffer, "%13[^,],PPM %d,%f,%f\n", targetRef, &ppmRef, &vmag, &minDistance);
            if (scanned != 4) {
                continue;
            }
            
            // Ignore if distance is too far of magnitude is absent
            if (minDistance < __FLT_EPSILON__ || minDistance > THRESHOLD_PPM || vmag < __FLT_EPSILON__) {
                continue;
            }
            
            // Find PPM star by ppmRef
            int ppmIndex = -1;
            for (int i = 0; i < PPMstars; i++) {
                if (PPMstar[i].ppmRef == ppmRef) {
                    ppmIndex = i;
                    break;
                }
            }
            
            // Ignore if PPM star not found or too bright
            if (ppmIndex == -1 || PPMstar[ppmIndex].vmag <= 1.0) {
                continue;
            }
            
            // Populate customStars histogram
            int index = (int)(10.0 * vmag + 0.5);
            if (index < 0 || index >= MAX_CUSTOM_STARS) {
                printf("Warning: Index out of range for vmag = %.2f (ppmRef = %d)\n", vmag, ppmRef);
                continue;
            }
            
            customStars[index].ppmVmag += PPMstar[ppmIndex].vmag;
            if (++customStars[index].count == 1) {
                customStarCount++;
            }
            matches++;
        }
        
        fclose(stream);

        // Discard those bins having less than 5 stars
        for (int i = 0; i < MAX_CUSTOM_STARS; i++) {
            if (customStars[i].count > 0 && customStars[i].count < 5) {
                customStars[i].ppmVmag = 0.0;
                customStars[i].count = 0;
                customStarCount--;
            }
        }
        printf("\nTotal matches processed: %d, custom stars registered: %d\n", matches, customStarCount);
        
    } else if (strcmp(mode, "--txt") == 0) {
        /* Mode is --txt, proceed with normal operation */
        const char *baseDir = (argc >= 3) ? argv[2] : ".";

        /* leemos Coord.txt y Flux.txt y buscamos coincidencias con PPM */
        printf("Reading Coord.txt and Flux.txt from base folder: %s\n", baseDir);

        char coordPath[256];
        char fluxPath[256];
        joinPath(baseDir, "Coord.txt", coordPath, sizeof(coordPath));
        joinPath(baseDir, "Flux.txt", fluxPath, sizeof(fluxPath));

        FILE *coordFile = fopen(coordPath, "rt");
        if (coordFile == NULL) {
            printf("Error: Cannot open %s file.\n", coordPath);
            return 1;
        }

        FILE *fluxFile = fopen(fluxPath, "rt");
        if (fluxFile == NULL) {
            printf("Error: Cannot open %s file.\n", fluxPath);
            fclose(coordFile);
            return 1;
        }

        /* Prepare Cross.txt output in the same base directory */
        char crossPath[256];
        joinPath(baseDir, "Cross.txt", crossPath, sizeof(crossPath));
        FILE *crossFile = fopen(crossPath, "wt");
        if (crossFile == NULL) {
            printf("Error: Cannot write %s file.\n", crossPath);
            fclose(coordFile);
            fclose(fluxFile);
            return 1;
        }
        fprintf(crossFile, "    Index   imag   Distance Identification\n");

        /* leemos catalogo PPM */
        printf("Reading PPM catalog...\n");
        readPPM(false, true, false, false, 2000.0);
        sortPPM();
        int PPMstars = getPPMStars();
        struct PPMstar_struct *PPMstar = getPPMStruct();

        /* leemos identificaciones cruzadas
            * Order de prioridad: UA, Lal, Lac, GC, ZC, OA, U, G */
        if (CROSS_OLD_CATALOGS) {
            readCrossFile("results/cross/cross_gilliss_ppm.csv", PPMstar, PPMstars);
            readCrossFile("results/cross/cross_usno_ppm.csv", PPMstar, PPMstars);
            //readCrossFile("results/cross/cross_oa_ppm.csv", PPMstar, PPMstars);
            //readCrossFile("results/cross/cross_zc_ppm.csv", PPMstar, PPMstars);
            readCrossFile("results/cross/cross_gc_ppm.csv", PPMstar, PPMstars);
            //readCrossFile("results/cross/cross_lalande_ppm.csv", PPMstar, PPMstars);
            //readCrossFile("results/cross/cross_lacaille_ppm.csv", PPMstar, PPMstars);
            readCrossFile("results/cross/cross_ua_ppm.csv", PPMstar, PPMstars);
        }

        char coordLine[256], fluxLine[256];
        int starIndex, fluxIndex;
        double RA, Decl, flux, background;
        int matches = 0;

        /* Skip header lines */
        if (fgets(coordLine, sizeof(coordLine), coordFile) == NULL) {
            printf("Error: Cannot read header from Coord.txt.\n");
            fclose(coordFile);
            fclose(fluxFile);
            fclose(crossFile);
            return 1;
        }

        if (fgets(fluxLine, sizeof(fluxLine), fluxFile) == NULL) {
            printf("Error: Cannot read header from Flux.txt.\n");
            fclose(coordFile);
            fclose(fluxFile);
            fclose(crossFile);
            return 1;
        }

        while (fgets(coordLine, sizeof(coordLine), coordFile) != NULL && 
                fgets(fluxLine, sizeof(fluxLine), fluxFile) != NULL) {
            
            int columnsCoordRead = sscanf(coordLine, "%d %lf %lf", &starIndex, &RA, &Decl);
            int columnsFluxRead = sscanf(fluxLine, "%d %lf %lf", &fluxIndex, &flux, &background);
            if (columnsCoordRead != 3 || columnsFluxRead != 3) {
                continue;
            }

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

            /* instrumental magnitude: zp - 2.5 log10(flux), for flux>=100 */        
            double imag = ZEROPOINT_IMAG - 2.5 * log10(fmax(flux, 100.0));

            /* Default write: unmatched */
            const char *identStr = "";
            char identBuf[128];

            /* Omit brightest stars (or without Vmag) and those already matched */
            bool ppmMatch = false;
            if (ppmIndex != -1
                    && minDistance < THRESHOLD_PPM
                    && !PPMstar[ppmIndex].discard) {
                snprintf(identBuf, sizeof(identBuf), "PPM %d / %s", PPMstar[ppmIndex].ppmRef, PPMstar[ppmIndex].dmString);
                identStr = identBuf;
                ppmMatch = true;
            }

            /* Write Cross.txt row, also save match in custom stars structure */
            if (ppmMatch) {
                fprintf(crossFile, "%8d  %6.2f  %8.1f  %s\n", starIndex, imag, minDistance, identStr);
                //printf("Star %d with flux=%.1f (imag=%.2f) is matched with PPM %d / %s (dist = %.1f arcsec.)\n", 
                //       starIndex, flux, imag, PPMstar[ppmIndex].ppmRef, PPMstar[ppmIndex].dmString, minDistance);
                
                /* Save match in custom stars structure, but only those no so bright or with photometric magnitude. */
                if (PPMstar[ppmIndex].vmag > 3.0) {
                    int index = (int)(10.0 * imag + 0.5);
                    if (index < 0 || index >= MAX_CUSTOM_STARS) {
                        printf("Warning: Index out of range at line %d (imag = %.2f)\n", starIndex, imag);
                        abort();
                    }
                    customStars[index].ppmVmag += PPMstar[ppmIndex].vmag;
                    if (++customStars[index].count == 1) {
                        customStarCount++;
                    }
                }
                
                PPMstar[ppmIndex].discard = true;
                matches++;
                continue;
            }
            
            // Write Cross.txt row, but without info about matching
            fprintf(crossFile, "%8d  %6.2f\n", starIndex, imag);
        }

        fclose(coordFile);
        fclose(fluxFile);
        fclose(crossFile);
        printf("\nTotal matches found: %d, custom stars registered: %d\n", matches, customStarCount);
    } else {
        printf("Error: Invalid mode '%s'. Use --csv or --txt\n", mode);
        return 1;
    }

    // Compute averages and display histogram
    printf("\nHistogram of custom stars:\n");
    printf("imag     ppmVmag     count\n");
    for (int i = 0; i < MAX_CUSTOM_STARS; i++) {
        if (customStars[i].count == 0) {
            continue;
        }
        double imag = i / 10.0;
        double avgPpmVmag = customStars[i].ppmVmag / customStars[i].count;
        printf("%-8.1f %-11.2f %d\n", imag, avgPpmVmag, customStars[i].count);
        customStars[i].ppmVmag = avgPpmVmag;
    }

    /* Perform constant regression fit between instrumental and PPM V magnitudes */
    if (customStarCount > 5) {
        printf("\nPerforming constant regression fit (assumes perfect linearity of instrumental magnitudes)...\n");
        
        /* Calculate constant fit: Vmag = zeroPoint + imag */
        double sum_vmag = 0.0, sum_imag = 0.0;
        double sum_weights = 0.0;
        
        for (int i = 0; i < MAX_CUSTOM_STARS; i++) {
            if (customStars[i].count == 0) {
                continue;
            }
            double ppm_vmag = customStars[i].ppmVmag;
            double imag = i / 10.0;
            double weight = log10(customStars[i].count);
            sum_vmag += ppm_vmag * weight;
            sum_imag += imag * weight;
            sum_weights += weight;
        }
        
        /* Calculate zero point */
        double zeroPoint = (sum_vmag - sum_imag) / sum_weights;
        
        printf(" 1) Constant fit: Vmag = %.3f + imag\n", zeroPoint);
        
        /* Calculate calibrated magnitudes, RMSE, and MAPE for constant fit */
        double sum_sq_error = 0.0;
        double sum_abs_percentage_error = 0.0;
        for (int i = 0; i < MAX_CUSTOM_STARS; i++) {
            if (customStars[i].count == 0) {
                continue;
            }
            double ppm_vmag = customStars[i].ppmVmag;
            double imag = i / 10.0;
            double calibrated_vmag = zeroPoint + imag;
            double error = fabs(calibrated_vmag - ppm_vmag);
            sum_sq_error += error * error;
            sum_abs_percentage_error += error / ppm_vmag * 100.0;
        }
        double rmse = sqrt(sum_sq_error / customStarCount);
        double mape = sum_abs_percentage_error / customStarCount;
        printf("\n    Constant fit RMSE: %.3f magnitudes\n", rmse);
        printf("    Constant fit MAPE: %.2f%%\n", mape);
        
        /* Now perform linear regression fit */
        printf("\nPerforming linear regression fit...\n");
        
        /* Calculate means for least squares fit */
        double sum_imag_lin = 0.0, sum_vmag_lin = 0.0;
        double sum_imag_vmag = 0.0, sum_imag_sq = 0.0;
        double sum_weights_lin = 0.0;
        
        for (int i = 0; i < MAX_CUSTOM_STARS; i++) {
            if (customStars[i].count == 0) {
                continue;
            }
            double imag = i / 10.0;
            double ppm_vmag = customStars[i].ppmVmag;
            double weight = log10(customStars[i].count);
            sum_imag_lin += imag * weight;
            sum_vmag_lin += ppm_vmag * weight;
            sum_imag_vmag += imag * ppm_vmag * weight;
            sum_imag_sq += imag * imag * weight;
            sum_weights_lin += weight;
        }
        
        double mean_imag = sum_imag_lin / sum_weights_lin;
        double mean_vmag = sum_vmag_lin / sum_weights_lin;
        
        /* Calculate slope and intercept */
        double denominator = sum_imag_sq - sum_weights_lin * mean_imag * mean_imag;
        double factor = (sum_imag_vmag - sum_weights_lin * mean_imag * mean_vmag) / denominator;
        double linearZeroPoint = mean_vmag - factor * mean_imag;
        
        printf(" 2) Linear fit: Vmag = %.3f + %.3f * imag\n", linearZeroPoint, factor);
        
        /* Calculate calibrated magnitudes, RMSE, and MAPE for linear fit */
        sum_sq_error = 0.0;
        sum_abs_percentage_error = 0.0;
        for (int i = 0; i < MAX_CUSTOM_STARS; i++) {
            if (customStars[i].count == 0) {
                continue;
            }
            double imag = i / 10.0;
            double ppm_vmag = customStars[i].ppmVmag;
            double calibrated_vmag = linearZeroPoint + factor * imag;
            double error = fabs(calibrated_vmag - ppm_vmag);
            sum_sq_error += error * error;
            sum_abs_percentage_error += error / ppm_vmag * 100.0;
        }
        rmse = sqrt(sum_sq_error / customStarCount);
        mape = sum_abs_percentage_error / customStarCount;
        printf("\n    Linear fit RMSE: %.3f magnitudes\n", rmse);
        printf("    Linear fit MAPE: %.2f%%\n", mape);
        
        /* Now perform quadratic regression fit (centered, solved stably) */
        printf("\nPerforming quadratic regression fit...\n");

        /* Center x to improve conditioning */
        double xmean = mean_imag;
        double S0 = 0.0;
        double Sx = 0.0, Sx2 = 0.0, Sx3 = 0.0, Sx4 = 0.0;
        double Sy = 0.0, Sxy = 0.0, Sx2y = 0.0;
        for (int i = 0; i < MAX_CUSTOM_STARS; i++) {
            if (customStars[i].count == 0) {
                continue;
            }
            double imag = i / 10.0;
            double x = imag - xmean;
            double y = customStars[i].ppmVmag;
            double weight = log10(customStars[i].count);
            double x2 = x * x;
            S0   += weight;
            Sx   += x * weight;
            Sx2  += x2 * weight;
            Sx3  += x2 * x * weight;
            Sx4  += x2 * x2 * weight;
            Sy   += y * weight;
            Sxy  += x * y * weight;
            Sx2y += x2 * y * weight;
        }

        /* Normal equations in centered basis: y = a + b*x + c*x^2 */
        double A[3][3] = {
            { S0,  Sx,  Sx2 },
            { Sx,  Sx2, Sx3 },
            { Sx2, Sx3, Sx4 }
        };
        double bvec[3] = { Sy, Sxy, Sx2y };
        double sol[3];
        if (!solve3x3(A, bvec, sol)) {
            printf("Error: Cannot perform quadratic regression (matrix ill-conditioned)\n");
        } else {
            /* Convert centered coefficients to original variable imag */
            double a_c = sol[0];
            double b_c = sol[1];
            double c_c = sol[2];
            double quadZeroPoint = a_c - b_c * xmean + c_c * xmean * xmean;
            double quadFactor    = b_c - 2.0 * c_c * xmean;
            double quad          = c_c;

            printf(" 3) Quadratic fit: Vmag = %.3f + %.3f * imag + %.3f * imag^2\n", quadZeroPoint, quadFactor, quad);

            /* Calculate calibrated magnitudes, RMSE, and MAPE for quadratic fit */
            sum_sq_error = 0.0;
            sum_abs_percentage_error = 0.0;
            for (int i = 0; i < MAX_CUSTOM_STARS; i++) {
                if (customStars[i].count == 0) {
                    continue;
                }
                double imag = i / 10.0;
                double trueV = customStars[i].ppmVmag;
                double yhat = quadZeroPoint + quadFactor * imag + quad * imag * imag;
                double err = fabs(yhat - trueV);
                sum_sq_error += err * err;
                sum_abs_percentage_error += err / trueV * 100.0;
            }

            rmse = sqrt(sum_sq_error / customStarCount);
            mape = sum_abs_percentage_error / customStarCount;
            printf("\n    Quadratic fit RMSE: %.3f magnitudes\n", rmse);
            printf("    Quadratic fit MAPE: %.2f%%\n", mape);
        }
    }

    return 0;
}
