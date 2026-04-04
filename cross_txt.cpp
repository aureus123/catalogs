/*
 * CROSS_TXT - Cross-identifies custom catalog against PPM catalog
 * Made in 2025 by Daniel E. Severin (partially created by Cursor AI and Claude AI)
 *
 * This tool cross-matches stars from a custom source (either extracted from
 * Coord.txt/Flux.txt files or from a pre-existing CSV cross-identification file)
 * against the PPM catalog, and performs photometric calibration by fitting
 * instrumental magnitudes to PPM V magnitudes.
 *
 * Three modes of operation:
 *
 * 1) --txt [baseDir]
 *    Reads Coord.txt (positions) and Flux.txt (fluxes) from the given base
 *    directory (defaults to current directory). For each star, computes its
 *    instrumental magnitude (imag = ZEROPOINT_IMAG - 2.5 * log10(flux)) and
 *    finds the nearest PPM match within THRESHOLD_PPM arcseconds. Writes all
 *    results to Cross.txt. Then builds a histogram of matched stars binned by
 *    instrumental magnitude (bins with fewer than 5 stars are discarded), and
 *    performs three regression fits (constant, linear, quadratic) of PPM V
 *    magnitude as a function of instrumental magnitude, reporting RMSE and MAPE
 *    for each. Finally writes Estim.txt with per-star estimated magnitudes
 *    from all three fits.
 *
 * 2) --csv csvFile
 *    Reads a pre-existing CSV cross-identification file (with columns:
 *    targetRef, "PPM" ppmRef, vmag, minDistance). The original visual
 *    magnitudes in the CSV are treated as instrumental magnitudes. Matches
 *    are looked up in the PPM catalog and the same histogram-based regression
 *    fits (constant, linear, quadratic) are performed as in --txt mode.
 *    Writes Estim.txt with per-star estimated magnitudes.
 *
 * 3) --var baseDir variablePPMid [controlPPMid] [sequencePPMid1] [sequencePPMid2]
 *              [controlVmag] [sequenceVmag1] [sequenceVmag2]
 *    Estimates the magnitude of a variable star. Reads Coord.txt/Flux.txt and
 *    performs cross-matching as in --txt mode (writing Cross.txt), but instead
 *    of histogram binning, it uses individual matched stars for calibration.
 *    The fitting strategy depends on the arguments provided:
 *
 *    a) Only variablePPMid given:
 *       Performs a linear fit (Vmag = a + b * imag) using an ensemble of all
 *       cross-matched PPM stars except the variable itself, with equal weights.
 *       Reports the fit, its RMSE and MAPE, and the estimated Vmag of the
 *       variable star.
 *
 *    b) controlPPMid also given:
 *       Same as (a), but additionally reports the estimated Vmag of the
 *       control star and the absolute error with respect to its PPM magnitude,
 *       as a quality check.
 *
 *    c) sequencePPMid1 also given:
 *       Instead of an ensemble fit, performs a constant fit (Vmag = zp + imag)
 *       derived from the single comparison star.
 *
 *    d) sequencePPMid2 also given:
 *       Instead of a constant fit, performs a linear fit (Vmag = a + b * imag)
 *       derived from the two comparison stars. Reports the variable and control
 *       star estimates, and indicates whether the variable star estimation is
 *       interpolated (its imag falls between the two sequence stars) or
 *       extrapolated.
 *
 *    Optionally, after all three PPM ids (cases b+c+d), exactly three Vmag
 *    values may be appended (controlVmag sequenceVmag1 sequenceVmag2). When
 *    present, these override the PPM catalog magnitudes for the control and
 *    sequence stars so that the fit is derived from user-supplied values
 *    instead.
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
#define CROSS_OLD_CATALOGS false

/* Structure to store custom star data (histogram).
 * Note: instrumental magnitudes are given by indices,
 * e.g. customStars[112] gives consolidated info about stars of 11.2 imag. */
struct CustomStars {
    double ppmVmag;
    int count;
};

/* Structure to store per-star data for Estim.txt output */
struct StarRecord {
    int index;
    double imag;
    double ppmVmag;  /* 0.0 if no match */
    bool hasMatch;
    char identification[128];
};

/*
 * readCrossFile - lee archivos de identificaciones cruzadas en formato CSV
 */
void readCrossFile(
        const char *ppm_file, struct PPMstar_struct *PPMstar, int PPMstars) {
    char buffer[1024], targetRef[STRING_SIZE], code[4];
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
        sscanf(buffer, "%13[^,],%3s %d,%f,%f\n", targetRef, code, &ppmRef, &vmag, &minDistance);
        if (strcmp(code, "PPM") != 0) continue;
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
        printf("Usage: %s --csv csvFile\n", argv[0]);
        printf("         It reads CSV cross-identification file and performs magnitude fits.\n");
        printf("Usage: %s --var baseDir variablePPMid [controlPPMid] [sequencePPMid1] [sequencePPMid2]\n", argv[0]);
        printf("              [controlVmag] [sequenceVmag1] [sequenceVmag2]\n");
        printf("         Estimates variable star magnitude from Coord.txt/Flux.txt data.\n");
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

    /* Per-star records for Estim.txt */
    struct StarRecord *starRecords = NULL;
    int starRecordCount = 0;
    int starRecordCapacity = 0;
    char estimPath[256];
    estimPath[0] = '\0';

    /* Fit parameters (set later if enough data) */
    bool fitComputed = false;
    double fit_constZP = 0.0;
    double fit_linZP = 0.0, fit_linSlope = 0.0;
    bool fit_quadOK = false;
    double fit_quadZP = 0.0, fit_quadSlope = 0.0, fit_quadCurve = 0.0;

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
        
        strcpy(estimPath, "Estim.txt");

        /* Load PPM catalog */
        printf("Reading PPM catalog...\n");
        readPPM(false, true, false, false, 2000.0);
        sortPPM();
        int PPMstars = getPPMStars();
        struct PPMstar_struct *PPMstar = getPPMStruct();
        
        /* Read CSV file and populate customStars histogram */
        char buffer[1024], targetRef[STRING_SIZE], code[4];
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
            
            sscanf(buffer, "%13[^,],%3s %d,%f,%f\n", targetRef, code, &ppmRef, &vmag, &minDistance);
            if (strcmp(code, "PPM") != 0) continue;
            
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

            /* Store star record for Estim.txt */
            if (starRecordCount >= starRecordCapacity) {
                starRecordCapacity = (starRecordCapacity == 0) ? 1024 : starRecordCapacity * 2;
                starRecords = (struct StarRecord *)realloc(starRecords, starRecordCapacity * sizeof(struct StarRecord));
            }
            {
                struct StarRecord *rec = &starRecords[starRecordCount++];
                rec->index = starRecordCount;
                rec->imag = vmag;
                rec->ppmVmag = PPMstar[ppmIndex].vmag;
                rec->hasMatch = true;
                snprintf(rec->identification, sizeof(rec->identification), "PPM %d", ppmRef);
            }
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
        joinPath(baseDir, "Estim.txt", estimPath, sizeof(estimPath));
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

            /* Store star record for Estim.txt */
            if (starRecordCount >= starRecordCapacity) {
                starRecordCapacity = (starRecordCapacity == 0) ? 1024 : starRecordCapacity * 2;
                starRecords = (struct StarRecord *)realloc(starRecords, starRecordCapacity * sizeof(struct StarRecord));
            }
            {
                struct StarRecord *rec = &starRecords[starRecordCount++];
                rec->index = starIndex;
                rec->imag = imag;
                if (ppmMatch) {
                    rec->ppmVmag = PPMstar[ppmIndex].vmag;
                    rec->hasMatch = true;
                    strncpy(rec->identification, identStr, sizeof(rec->identification) - 1);
                    rec->identification[sizeof(rec->identification) - 1] = '\0';
                } else {
                    rec->ppmVmag = 0.0;
                    rec->hasMatch = false;
                    rec->identification[0] = '\0';
                }
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
    } else if (strcmp(mode, "--var") == 0) {
        /* Mode is --var, estimate variable star magnitude */
        if (argc < 4) {
            printf("Error: baseDir and variablePPMid are required for --var mode.\n");
            return 1;
        }

        const char *baseDir = argv[2];
        int variablePPMid = atoi(argv[3]);
        int controlPPMid = (argc >= 5) ? atoi(argv[4]) : 0;
        int sequencePPMid1 = (argc >= 6) ? atoi(argv[5]) : 0;
        int sequencePPMid2 = (argc >= 7) ? atoi(argv[6]) : 0;
        bool hasControl = (argc >= 5);
        bool hasSeq1 = (argc >= 6);
        bool hasSeq2 = (argc >= 7);

        /* Custom vmags: if given, must be exactly 3 (positions 7, 8, 9 in argv) */
        if (argc > 7 && argc < 10) {
            printf("Error: Custom vmags must be provided as a group of exactly 3 "
                   "(controlVmag sequenceVmag1 sequenceVmag2).\n");
            return 1;
        }
        bool hasCustomVmags = (argc >= 10);
        double customControlVmag = hasCustomVmags ? atof(argv[7]) : 0.0;
        double customSeq1Vmag    = hasCustomVmags ? atof(argv[8]) : 0.0;
        double customSeq2Vmag    = hasCustomVmags ? atof(argv[9]) : 0.0;

        if (hasCustomVmags && !hasSeq2) {
            printf("Error: Custom vmags require all three PPM ids "
                   "(controlPPMid, sequencePPMid1, sequencePPMid2) to be given.\n");
            return 1;
        }

        printf("Variable star estimation mode.\n");
        printf("  Variable:   PPM %d\n", variablePPMid);
        if (hasControl) printf("  Control:    PPM %d\n", controlPPMid);

        /* Read Coord.txt and Flux.txt */
        char coordPath[256], fluxPath[256], crossPath[256];
        joinPath(baseDir, "Coord.txt", coordPath, sizeof(coordPath));
        joinPath(baseDir, "Flux.txt", fluxPath, sizeof(fluxPath));
        joinPath(baseDir, "Cross.txt", crossPath, sizeof(crossPath));

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
        FILE *crossFile = fopen(crossPath, "wt");
        if (crossFile == NULL) {
            printf("Error: Cannot write %s file.\n", crossPath);
            fclose(coordFile);
            fclose(fluxFile);
            return 1;
        }
        fprintf(crossFile, "    Index   imag   Distance Identification\n");

        /* Read PPM catalog */
        printf("Reading PPM catalog...\n");
        readPPM(false, true, false, false, 2000.0);
        sortPPM();
        int PPMstars = getPPMStars();
        struct PPMstar_struct *PPMstar = getPPMStruct();

        if (CROSS_OLD_CATALOGS) {
            readCrossFile("results/cross/cross_gilliss_ppm.csv", PPMstar, PPMstars);
            readCrossFile("results/cross/cross_usno_ppm.csv", PPMstar, PPMstars);
            readCrossFile("results/cross/cross_gc_ppm.csv", PPMstar, PPMstars);
            readCrossFile("results/cross/cross_ua_ppm.csv", PPMstar, PPMstars);
        }

        /* Per-match storage for variable star analysis */
        struct VarMatch { int starIndex; double imag; int ppmRef; double ppmVmag; };
        struct VarMatch *varMatches = NULL;
        int varMatchCount = 0, varMatchCapacity = 0;

        char coordLine[256], fluxLine[256];
        int starIndex, fluxIndex;
        double RA, Decl, flux, background;
        int matches = 0;

        /* Skip headers */
        if (fgets(coordLine, sizeof(coordLine), coordFile) == NULL) {
            printf("Error: Cannot read header from Coord.txt.\n");
            fclose(coordFile); fclose(fluxFile); fclose(crossFile);
            return 1;
        }
        if (fgets(fluxLine, sizeof(fluxLine), fluxFile) == NULL) {
            printf("Error: Cannot read header from Flux.txt.\n");
            fclose(coordFile); fclose(fluxFile); fclose(crossFile);
            return 1;
        }

        while (fgets(coordLine, sizeof(coordLine), coordFile) != NULL &&
               fgets(fluxLine, sizeof(fluxLine), fluxFile) != NULL) {
            int columnsCoordRead = sscanf(coordLine, "%d %lf %lf", &starIndex, &RA, &Decl);
            int columnsFluxRead = sscanf(fluxLine, "%d %lf %lf", &fluxIndex, &flux, &background);
            if (columnsCoordRead != 3 || columnsFluxRead != 3) continue;
            if (starIndex != fluxIndex) {
                printf("Warning: Index mismatch (Coord: %d, Flux: %d)\n", starIndex, fluxIndex);
                continue;
            }

            double x, y, z;
            sph2rec(RA, Decl, &x, &y, &z);

            int ppmIndex = -1;
            double minDistance = HUGE_NUMBER;
            findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);

            double imag = ZEROPOINT_IMAG - 2.5 * log10(fmax(flux, 100.0));

            bool ppmMatch = false;
            char identBuf[128];
            if (ppmIndex != -1 && minDistance < THRESHOLD_PPM && !PPMstar[ppmIndex].discard) {
                snprintf(identBuf, sizeof(identBuf), "PPM %d / %s",
                         PPMstar[ppmIndex].ppmRef, PPMstar[ppmIndex].dmString);
                ppmMatch = true;
            }

            if (ppmMatch) {
                fprintf(crossFile, "%8d  %6.2f  %8.1f  %s\n", starIndex, imag, minDistance, identBuf);

                /* Store match */
                if (varMatchCount >= varMatchCapacity) {
                    varMatchCapacity = (varMatchCapacity == 0) ? 1024 : varMatchCapacity * 2;
                    varMatches = (struct VarMatch *)realloc(varMatches,
                            varMatchCapacity * sizeof(struct VarMatch));
                }
                varMatches[varMatchCount].starIndex = starIndex;
                varMatches[varMatchCount].imag = imag;
                varMatches[varMatchCount].ppmRef = PPMstar[ppmIndex].ppmRef;
                varMatches[varMatchCount].ppmVmag = PPMstar[ppmIndex].vmag;
                varMatchCount++;

                PPMstar[ppmIndex].discard = true;
                matches++;
                continue;
            }

            fprintf(crossFile, "%8d  %6.2f\n", starIndex, imag);
        }

        fclose(coordFile);
        fclose(fluxFile);
        fclose(crossFile);
        printf("\nTotal matches found: %d\n", matches);

        /* Find variable star in matches */
        int varIdx = -1;
        for (int i = 0; i < varMatchCount; i++) {
            if (varMatches[i].ppmRef == variablePPMid) { varIdx = i; break; }
        }
        if (varIdx == -1) {
            printf("Error: Variable star PPM %d not found in cross-identification.\n", variablePPMid);
            free(varMatches);
            return 1;
        }

        /* Find control star */
        int ctrlIdx = -1;
        if (hasControl) {
            for (int i = 0; i < varMatchCount; i++) {
                if (varMatches[i].ppmRef == controlPPMid) { ctrlIdx = i; break; }
            }
            if (ctrlIdx == -1) {
                printf("Error: Control star PPM %d not found in cross-identification.\n", controlPPMid);
                free(varMatches);
                return 1;
            }
        }

        /* Find sequence stars */
        int seq1Idx = -1, seq2Idx = -1;
        if (hasSeq1) {
            for (int i = 0; i < varMatchCount; i++) {
                if (varMatches[i].ppmRef == sequencePPMid1) { seq1Idx = i; break; }
            }
            if (seq1Idx == -1) {
                printf("Error: Sequence star 1 (PPM %d) not found in cross-identification.\n", sequencePPMid1);
                free(varMatches);
                return 1;
            }
        }
        if (hasSeq2) {
            for (int i = 0; i < varMatchCount; i++) {
                if (varMatches[i].ppmRef == sequencePPMid2) { seq2Idx = i; break; }
            }
            if (seq2Idx == -1) {
                printf("Error: Sequence star 2 (PPM %d) not found in cross-identification.\n", sequencePPMid2);
                free(varMatches);
                return 1;
            }
        }

        /* Override PPM vmags with user-supplied values if provided */
        if (hasCustomVmags) {
            varMatches[ctrlIdx].ppmVmag  = customControlVmag;
            varMatches[seq1Idx].ppmVmag  = customSeq1Vmag;
            varMatches[seq2Idx].ppmVmag  = customSeq2Vmag;
            printf("  (Using user-supplied Vmag values instead of PPM catalog)\n");
        }

        if (hasSeq1)
            printf("  Sequence 1: PPM %d (Vmag = %.1f)\n", sequencePPMid1, varMatches[seq1Idx].ppmVmag);
        if (hasSeq2)
            printf("  Sequence 2: PPM %d (Vmag = %.1f)\n", sequencePPMid2, varMatches[seq2Idx].ppmVmag);

        /* Perform fitting */
        double fitIntercept = 0.0, fitSlope = 1.0;

        if (!hasSeq1) {
            /* Ensemble linear fit: all matched except variable, weight=1.
             * Exclude stars with no visual magnitude (vmag=0.0) or too bright
             * (vmag <= 3.0), same criteria as the --txt histogram. */
            double si = 0.0, sv = 0.0, siv = 0.0, sii = 0.0;
            int n = 0;
            for (int i = 0; i < varMatchCount; i++) {
                if (i == varIdx) continue;
                if (varMatches[i].ppmVmag <= 3.0) continue;
                double im = varMatches[i].imag;
                double vm = varMatches[i].ppmVmag;
                si += im;
                sv += vm;
                siv += im * vm;
                sii += im * im;
                n++;
            }
            if (n < 2) {
                printf("Error: Not enough stars for ensemble fit (need at least 2, got %d).\n", n);
                free(varMatches);
                return 1;
            }
            double mi = si / n, mv = sv / n;
            double denom = sii - n * mi * mi;
            fitSlope = (siv - n * mi * mv) / denom;
            fitIntercept = mv - fitSlope * mi;

            printf("\nEnsemble linear fit: Vmag = %.3f + %.3f * imag\n", fitIntercept, fitSlope);

            /* RMSE and MAPE */
            double sse = 0.0, sape = 0.0;
            for (int i = 0; i < varMatchCount; i++) {
                if (i == varIdx) continue;
                if (varMatches[i].ppmVmag <= 3.0) continue;
                double est = fitIntercept + fitSlope * varMatches[i].imag;
                double err = fabs(est - varMatches[i].ppmVmag);
                sse += err * err;
                sape += err / varMatches[i].ppmVmag * 100.0;
            }
            printf("    Ensemble RMSE: %.3f magnitudes\n", sqrt(sse / n));
            printf("    Ensemble MAPE: %.2f%%\n", sape / n);

        } else if (!hasSeq2) {
            /* Constant fit with one sequence star */
            fitSlope = 1.0;
            fitIntercept = varMatches[seq1Idx].ppmVmag - varMatches[seq1Idx].imag;

            printf("\nConstant fit: Vmag = %.3f + imag\n", fitIntercept);

        } else {
            /* Linear fit with two sequence stars */
            double s1i = varMatches[seq1Idx].imag, s1v = varMatches[seq1Idx].ppmVmag;
            double s2i = varMatches[seq2Idx].imag, s2v = varMatches[seq2Idx].ppmVmag;

            if (fabs(s2i - s1i) < 1e-9) {
                printf("Error: Sequence stars have identical instrumental magnitudes.\n");
                free(varMatches);
                return 1;
            }

            fitSlope = (s2v - s1v) / (s2i - s1i);
            fitIntercept = s1v - fitSlope * s1i;

            printf("\nLinear fit: Vmag = %.3f + %.3f * imag\n", fitIntercept, fitSlope);
        }

        /* Variable star estimation */
        double varEst = fitIntercept + fitSlope * varMatches[varIdx].imag;
        printf("\nVariable star PPM %d: imag = %.2f, estimated Vmag = %.2f",
               variablePPMid, varMatches[varIdx].imag, varEst);

        if (hasSeq2) {
            double minImag = fmin(varMatches[seq1Idx].imag, varMatches[seq2Idx].imag);
            double maxImag = fmax(varMatches[seq1Idx].imag, varMatches[seq2Idx].imag);
            if (varMatches[varIdx].imag >= minImag && varMatches[varIdx].imag <= maxImag) {
                printf(" (interpolated)");
            } else {
                printf(" (extrapolated)");
            }
        }
        printf("\n");

        /* Control star estimation */
        if (hasControl) {
            double ctrlEst = fitIntercept + fitSlope * varMatches[ctrlIdx].imag;
            double ctrlErr = fabs(ctrlEst - varMatches[ctrlIdx].ppmVmag);
            printf("Control star PPM %d: imag = %.2f, estimated Vmag = %.2f, PPM Vmag = %.1f, error = %.2f mag\n",
                   controlPPMid, varMatches[ctrlIdx].imag, ctrlEst, varMatches[ctrlIdx].ppmVmag, ctrlErr);
        }

        free(varMatches);
        return 0;

    } else {
        printf("Error: Invalid mode '%s'. Use --csv, --txt, or --var\n", mode);
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
        fit_constZP = zeroPoint;
        fitComputed = true;

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
        fit_linZP = linearZeroPoint;
        fit_linSlope = factor;

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
            fit_quadZP = quadZeroPoint;
            fit_quadSlope = quadFactor;
            fit_quadCurve = quad;
            fit_quadOK = true;

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

    /* Generate Estim.txt with estimated magnitudes */
    if (fitComputed && starRecordCount > 0 && estimPath[0] != '\0') {
        FILE *estimFile = fopen(estimPath, "wt");
        if (estimFile == NULL) {
            printf("Error: Cannot write %s file.\n", estimPath);
        } else {
            fprintf(estimFile, "   Index    Vmag  ConstEst    LinEst   QuadEst  Identification\n");
            for (int i = 0; i < starRecordCount; i++) {
                struct StarRecord *rec = &starRecords[i];
                double constEst = fit_constZP + rec->imag;
                double linEst = fit_linZP + fit_linSlope * rec->imag;

                fprintf(estimFile, "%8d", rec->index);

                if (rec->hasMatch && rec->ppmVmag > __FLT_EPSILON__) {
                    fprintf(estimFile, "  %6.1f", rec->ppmVmag);
                } else {
                    fprintf(estimFile, "        ");
                }

                fprintf(estimFile, "  %8.2f  %8.2f", constEst, linEst);

                if (fit_quadOK) {
                    double quadEst = fit_quadZP + fit_quadSlope * rec->imag + fit_quadCurve * rec->imag * rec->imag;
                    fprintf(estimFile, "  %8.2f", quadEst);
                }

                if (rec->hasMatch) {
                    fprintf(estimFile, "  %s", rec->identification);
                }

                fprintf(estimFile, "\n");
            }
            fclose(estimFile);
            printf("\nEstim.txt written to %s with %d stars.\n", estimPath, starRecordCount);
        }
    }

    free(starRecords);
    return 0;
}
