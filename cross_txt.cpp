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
#define MAX_CUSTOM_STARS 10000
#define THRESHOLD_PPM 15.0
#define THRESHOLD_MAG 1.5
#define ZEROPOINT_IMAG 17.37
#define CROSS_OLD_CATALOGS true

/* Structure to store custom star data */
struct CustomStars {
    int index;
    int ppmIndex;
    double RA;
    double Decl;
    double imag;  /* instrumental magnitude */
};

/* Small 3x3 linear system solver with partial pivoting. Returns 1 on success, 0 on singular. */
static int solve3x3(double A[3][3], double b[3], double x[3]) {
    int p[3] = {0, 1, 2};
    // Partial pivoting for column 0
    int maxr = 0;
    if (fabs(A[1][0]) > fabs(A[maxr][0])) maxr = 1;
    if (fabs(A[2][0]) > fabs(A[maxr][0])) maxr = 2;
    if (maxr != 0) {
        for (int j = 0; j < 3; ++j) { double t = A[0][j]; A[0][j] = A[maxr][j]; A[maxr][j] = t; }
        double tb = b[0]; b[0] = b[maxr]; b[maxr] = tb;
    }
    if (fabs(A[0][0]) < 1e-15) return 0;
    double m10 = A[1][0] / A[0][0];
    double m20 = A[2][0] / A[0][0];
    for (int j = 1; j < 3; ++j) {
        A[1][j] -= m10 * A[0][j];
        A[2][j] -= m20 * A[0][j];
    }
    b[1] -= m10 * b[0];
    b[2] -= m20 * b[0];

    // Pivot for column 1
    maxr = 1;
    if (fabs(A[2][1]) > fabs(A[maxr][1])) maxr = 2;
    if (maxr != 1) {
        for (int j = 1; j < 3; ++j) { double t = A[1][j]; A[1][j] = A[maxr][j]; A[maxr][j] = t; }
        double tb = b[1]; b[1] = b[maxr]; b[maxr] = tb;
    }
    if (fabs(A[1][1]) < 1e-15) return 0;
    double m21 = A[2][1] / A[1][1];
    A[2][2] -= m21 * A[1][2];
    b[2] -= m21 * b[1];

    if (fabs(A[2][2]) < 1e-15) return 0;

    // Back substitution
    x[2] = b[2] / A[2][2];
    x[1] = (b[1] - A[1][2] * x[2]) / A[1][1];
    x[0] = (b[0] - A[0][1] * x[1] - A[0][2] * x[2]) / A[0][0];
    return 1;
}

/*
 * readCrossFile - lee archivos de identificaciones cruzadas en formato CSV
 */
void readCrossFile(
        const char *ppm_file, struct PPMstar_struct *PPMstar, int PPMstars) {
    char buffer[1024];
    char targetRef[STRING_SIZE];
    int ppmRef, cdDeclRef, cdNumRef, cpdDeclRef, cpdNumRef;
    float minDistance;

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
        sscanf(buffer, "%13[^,],PPM %d,%f\n", targetRef, &ppmRef, &minDistance);
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

    const char *baseDir = (argc >= 2) ? argv[1] : ".";

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

    /* Array to store custom stars data */
    struct CustomStars customStars[MAX_CUSTOM_STARS];
    int customStarCount = 0;
    
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
            if (PPMstar[ppmIndex].vmag > 3.0 &&
                    customStarCount < MAX_CUSTOM_STARS) {
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
        
        // Write Cross.txt row, but without info about matching
        fprintf(crossFile, "%8d  %6.2f\n", starIndex, imag);
    }
    
    fclose(coordFile);
    fclose(fluxFile);
    fclose(crossFile);
    printf("\nTotal matches found: %d\n", matches);

    /* Perform constant regression fit between instrumental and PPM V magnitudes */
    if (customStarCount > 5) {
        printf("\nPerforming constant regression fit (assumes perfect linearity of instrumental magnitudes)...\n");
        
        /* Calculate constant fit: Vmag = zeroPoint + imag */
        double sum_vmag = 0.0, sum_imag = 0.0;
        
        for (int i = 0; i < customStarCount; i++) {
            int ppmIndex = customStars[i].ppmIndex;
            sum_vmag += PPMstar[ppmIndex].vmag;
            sum_imag += customStars[i].imag;
        }
        
        /* Calculate zero point */
        double zeroPoint = (sum_vmag - sum_imag) / customStarCount;
        
        printf("Constant fit: Vmag = %.3f + imag\n", zeroPoint);
        
        /* Calculate calibrated magnitudes, RMSE, and MAPE for constant fit */
        double sum_sq_error = 0.0;
        double sum_abs_percentage_error = 0.0;
        printf("\nConstant fit worst calibrated magnitudes:\n");
        
        for (int i = 0; i < customStarCount; i++) {
            int ppmIndex = customStars[i].ppmIndex;
            double ppm_vmag = PPMstar[ppmIndex].vmag;
            double calibrated_vmag = zeroPoint + customStars[i].imag;
            double error = fabs(calibrated_vmag - ppm_vmag);
            sum_sq_error += error * error;
            sum_abs_percentage_error += error / ppm_vmag * 100.0;

            if (error > THRESHOLD_MAG) {
                printf("Star %d has Vmag=%.2f while PPM %d / %s has Vmag=%.2f, difference=%.2f\n", 
                        customStars[i].index, calibrated_vmag, PPMstar[ppmIndex].ppmRef, PPMstar[ppmIndex].dmString, ppm_vmag, error);
            }
        }
        
        double rmse = sqrt(sum_sq_error / customStarCount);
        double mape = sum_abs_percentage_error / customStarCount;
        printf("\nConstant fit RMSE: %.3f magnitudes\n", rmse);
        printf("Constant fit MAPE: %.2f%%\n", mape);
        
        /* Now perform linear regression fit */
        printf("\nPerforming linear regression fit...\n");
        
        /* Calculate means for least squares fit */
        double sum_imag_lin = 0.0, sum_vmag_lin = 0.0;
        double sum_imag_vmag = 0.0, sum_imag_sq = 0.0;
        
        for (int i = 0; i < customStarCount; i++) {
            int ppmIdx = customStars[i].ppmIndex;
            sum_imag_lin += customStars[i].imag;
            sum_vmag_lin += PPMstar[ppmIdx].vmag;
            sum_imag_vmag += customStars[i].imag * PPMstar[ppmIdx].vmag;
            sum_imag_sq += customStars[i].imag * customStars[i].imag;
        }
        
        double mean_imag = sum_imag_lin / customStarCount;
        double mean_vmag = sum_vmag_lin / customStarCount;
        
        /* Calculate slope and intercept */
        double denominator = sum_imag_sq - customStarCount * mean_imag * mean_imag;
        double factor = (sum_imag_vmag - customStarCount * mean_imag * mean_vmag) / denominator;
        double linearZeroPoint = mean_vmag - factor * mean_imag;
        
        printf("Linear fit: Vmag = %.3f + %.3f * imag\n", linearZeroPoint, factor);
        
        /* Calculate calibrated magnitudes, RMSE, and MAPE for linear fit */
        sum_sq_error = 0.0;
        sum_abs_percentage_error = 0.0;
        printf("\nLinear fit worst calibrated magnitudes:\n");
        
        for (int i = 0; i < customStarCount; i++) {
            int ppmIndex = customStars[i].ppmIndex;
            double ppm_vmag = PPMstar[ppmIndex].vmag;
            double calibrated_vmag = linearZeroPoint + factor * customStars[i].imag;
            double error = fabs(calibrated_vmag - ppm_vmag);
            sum_sq_error += error * error;
            sum_abs_percentage_error += error / ppm_vmag * 100.0;

            if (error > THRESHOLD_MAG) {
                printf("Star %d has Vmag=%.2f while PPM %d / %s has Vmag=%.2f, difference=%.2f\n", 
                        customStars[i].index, calibrated_vmag, PPMstar[ppmIndex].ppmRef, PPMstar[ppmIndex].dmString, ppm_vmag, error);
            }
        }
        
        rmse = sqrt(sum_sq_error / customStarCount);
        mape = sum_abs_percentage_error / customStarCount;
        printf("\nLinear fit RMSE: %.3f magnitudes\n", rmse);
        printf("Linear fit MAPE: %.2f%%\n", mape);
        
        /* Now perform quadratic regression fit (centered, solved stably) */
        printf("\nPerforming quadratic regression fit...\n");

        /* Center x to improve conditioning */
        double xmean = mean_imag;
        double S0 = (double)customStarCount;
        double Sx = 0.0, Sx2 = 0.0, Sx3 = 0.0, Sx4 = 0.0;
        double Sy = 0.0, Sxy = 0.0, Sx2y = 0.0;
        for (int i = 0; i < customStarCount; i++) {
            double x = customStars[i].imag - xmean;
            double y = PPMstar[customStars[i].ppmIndex].vmag;
            double x2 = x * x;
            Sx   += x;
            Sx2  += x2;
            Sx3  += x2 * x;
            Sx4  += x2 * x2;
            Sy   += y;
            Sxy  += x * y;
            Sx2y += x2 * y;
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

            printf("Quadratic fit: Vmag = %.3f + %.3f * imag + %.3f * imag^2\n", quadZeroPoint, quadFactor, quad);

            /* Calculate calibrated magnitudes, RMSE, and MAPE for quadratic fit */
            sum_sq_error = 0.0;
            sum_abs_percentage_error = 0.0;
            printf("\nQuadratic fit worst calibrated magnitudes:\n");

            for (int i = 0; i < customStarCount; i++) {
                int ppmIndex = customStars[i].ppmIndex;
                double trueV = PPMstar[ppmIndex].vmag;
                double x = customStars[i].imag;
                double yhat = quadZeroPoint + quadFactor * x + quad * x * x;
                double err = fabs(yhat - trueV);
                sum_sq_error += err * err;
                sum_abs_percentage_error += err / trueV * 100.0;

                if (err > THRESHOLD_MAG) {
                    printf("Star %d has Vmag=%.2f while PPM %d / %s has Vmag=%.2f, difference=%.2f\n",
                           customStars[i].index, yhat, PPMstar[ppmIndex].ppmRef, PPMstar[ppmIndex].dmString, trueV, err);
                }
            }

            rmse = sqrt(sum_sq_error / customStarCount);
            mape = sum_abs_percentage_error / customStarCount;
            printf("\nQuadratic fit RMSE: %.3f magnitudes\n", rmse);
            printf("Quadratic fit MAPE: %.2f%%\n", mape);
        }
    }

    return 0;
} 