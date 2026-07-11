
/*
 * CROSS_UTILS - Codigo auxiliar compartido por cross_north y cross_south
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_ppm.h"
#include "read_cpd.h"
#include "read_dm.h"
#include "read_gc.h"
#include "trig.h"
#include "misc.h"
#include "find_gsc.h"
#include "cross_utils.h"

/*
 * preparePPM - lee PPM a la epoca dada y lo deja ordenado
 */
struct PPMstar_struct *preparePPM(double epoch, bool discardSouth) {
    readPPM(false, true, false, discardSouth, epoch);
    sortPPM();
    return getPPMStruct();
}

/*
 * readRAField - lee ascension recta y la devuelve en grados
 */
double readRAField(char *buffer, int col, int *RAh, int *RAm, int *RAs) {
    char cell[256];

    readFieldSanitized(buffer, cell, col, 2);
    *RAh = atoi(cell);
    double RA = (double) *RAh;
    readFieldSanitized(buffer, cell, col + 2, 2);
    *RAm = atoi(cell);
    RA += ((double) *RAm)/60.0;
    readFieldSanitized(buffer, cell, col + 4, 4);
    *RAs = atoi(cell);
    RA += (((double) *RAs)/100.0)/3600.0;
    return RA * 15.0; /* conversion horas a grados */
}

/*
 * readDeclField - lee declinacion (o NPD/SPD) y la devuelve en grados positivos
 */
double readDeclField(char *buffer, int colD, int widthD, int colM, int colS, int widthS,
        double divisor, int *Decld, int *Declm, int *Decls) {
    char cell[256];

    readFieldSanitized(buffer, cell, colD, widthD);
    *Decld = atoi(cell);
    double Decl = (double) *Decld;
    readFieldSanitized(buffer, cell, colM, 2);
    *Declm = atoi(cell);
    Decl += ((double) *Declm)/60.0;
    readFieldSanitized(buffer, cell, colS, widthS);
    *Decls = atoi(cell);
    Decl += (((double) *Decls)/divisor)/3600.0;
    return Decl;
}

/*
 * readMagIntHalf - lee magnitud entera con posible fraccion '5'
 */
float readMagIntHalf(char *buffer, int col, int width, int colHalf) {
    char cell[256];

    readField(buffer, cell, col, width);
    float vmag = 0.0;
    if (cell[0] != ' ') {
        vmag = (float) atoi(cell);
        readField(buffer, cell, colHalf, 1);
        if (cell[0] == '5') vmag += 0.5;
    }
    return vmag;
}

/*
 * formatCatLine - formatea una linea de registro
 */
void formatCatLine(char *catLine, int RAh, int RAm, int RAs, char sign, int Decld, int Declm, int Decls) {
    snprintf(catLine, 64, "%02dh %02dm %02ds%02d %c%02d°%02d'%02d\"%01d",
        RAh, RAm, RAs / 100, RAs % 100, sign, Decld, Declm, Decls / 10, Decls % 10);
}

/*
 * openCrossSet / closeCrossSet - terna de archivos de cruzamiento PPM/SAO/HD
 */
void openCrossSet(const char *tag, FILE **ppm, FILE **sao, FILE **hd) {
    char name[64];

    snprintf(name, 64, "results/cross/cross_%s_ppm.csv", tag);
    *ppm = openCrossFile(name);
    snprintf(name, 64, "results/cross/cross_%s_sao.csv", tag);
    *sao = openCrossFile(name);
    snprintf(name, 64, "results/cross/cross_%s_hd.csv", tag);
    *hd = openCrossFile(name);
}

void closeCrossSet(FILE *ppm, FILE *sao, FILE *hd) {
    fclose(ppm);
    fclose(sao);
    fclose(hd);
}

/*
 * crossWithPPM - busca la PPM mas cercana y genera el cruzamiento
 */
bool crossWithPPM(double x, double y, double z, double Decl, double vmag, double maxDist,
        const char *magWarnName, char *crossName, FILE *ppmStream, FILE *saoStream, FILE *hdStream,
        int *ppmIndexOut, double *nearestPPMDistance, struct CrossStats *stats) {
    struct PPMstar_struct *PPMstar = getPPMStruct();

    int ppmIndex = -1;
    double minDistance = HUGE_NUMBER;
    findPPMByCoordinates(x, y, z, Decl, &ppmIndex, &minDistance);
    *ppmIndexOut = ppmIndex;
    *nearestPPMDistance = minDistance;
    if (minDistance >= maxDist) return false;

    float vmagf = (float) vmag;
    if (magWarnName != NULL) {
        float ppmVmag = PPMstar[ppmIndex].vmag;
        if (vmagf > __FLT_EPSILON__ && fabs(ppmVmag) > __FLT_EPSILON__) {
            // Note: no fit is performed to convert scales of magnitudes
            float delta = fabs(vmagf - ppmVmag);
            if (delta < MAX_MAGNITUDE) {
                stats->akkuDeltaError += delta * delta;
                stats->countDelta++;
            } else {
                printf("%d) Warning: %s (vmag = %.1f) has a near PPM %d star (vmag = %.1f) with a different magnitude (delta = %.1f).\n",
                    ++stats->errors,
                    magWarnName,
                    vmagf,
                    PPMstar[ppmIndex].ppmRef,
                    ppmVmag,
                    delta);
            }
        }
    }

    stats->akkuDistError += minDistance * minDistance;
    stats->countDist++;
    if (crossName != NULL) {
        writePPMCrossEntry(ppmStream, saoStream, hdStream, crossName, &PPMstar[ppmIndex], vmag, minDistance);
    }
    return true;
}

/*
 * tryGSC - si no hubo PPM cercana, prueba con GSC
 */
bool tryGSC(bool ppmFound, double RA, double Decl, double epoch,
        FILE *crossPPMStream, char *catName, double vmag, struct CrossStats *stats) {
    if (ppmFound) return false;
    bool gscFound = findGSCStar(RA, Decl, epoch, MAX_DIST_GSC);
    if (gscFound && crossPPMStream != NULL) {
        char gscName[20];
        snprintf(gscName, 20, "GSC %s", getGSCId());
        writeCrossEntry(crossPPMStream, catName, gscName, vmag, getDist());
        stats->countGSC++;
    }
    return gscFound;
}

/*
 * crossWithCPD - busca la CPD mas cercana y genera el cruzamiento
 */
bool crossWithCPD(double x, double y, double z, double Decl1875, double vmag, double maxDist,
        FILE *stream, char *catName, struct CrossStats *stats, int *cpdIndexOut, double *minDistOut) {
    struct CPDstar_struct *CPDstar = getCPDStruct();

    int cpdIndex = -1;
    double minDistance = HUGE_NUMBER;
    findCPDByCoordinates(x, y, z, Decl1875, &cpdIndex, &minDistance);
    if (cpdIndexOut != NULL) *cpdIndexOut = cpdIndex;
    if (minDistOut != NULL) *minDistOut = minDistance;
    if (minDistance >= maxDist) return false;

    stats->countCPD++;
    if (catName != NULL) {
        char cpdName[20];
        snprintf(cpdName, 20, "CPD %d°%d", CPDstar[cpdIndex].declRef, CPDstar[cpdIndex].numRef);
        writeCrossEntry(stream, catName, cpdName, vmag, minDistance);
    }
    return true;
}

/*
 * crossWithCD - busca la CD mas cercana, compara magnitudes y genera el cruzamiento
 */
bool crossWithCD(double x, double y, double z, double Decl1875, double vmag,
        const char *magWarnName, FILE *stream, char *catName, struct CrossStats *stats,
        int *cdIndexOut, double *minDistOut) {
    struct DMstar_struct *CDstar = getDMStruct();

    int cdIndex = -1;
    double minDistance = HUGE_NUMBER;
    findDMByCoordinates(x, y, z, Decl1875, &cdIndex, &minDistance);
    if (cdIndexOut != NULL) *cdIndexOut = cdIndex;
    if (minDistOut != NULL) *minDistOut = minDistance;
    if (minDistance >= MAX_DIST_CD) return false;

    float vmagf = (float) vmag;
    float cdVmag = CDstar[cdIndex].vmag;
    if (magWarnName != NULL && vmagf > __FLT_EPSILON__ && cdVmag < 29.9) {
        float delta = fabs(vmagf - cdVmag);
        if (delta >= MAX_MAGNITUDE) {
            printf("%d) Warning: %s (vmag = %.1f) has a near CD with dif. magnitude (delta = %.1f). Check dpl.\n",
                ++stats->errors,
                magWarnName,
                vmagf,
                delta);
            writeRegister(cdIndex, false);
        }
    }

    stats->countCD++;
    if (catName != NULL) {
        char cdName[20];
        snprintf(cdName, 20, "CD %d°%d", CDstar[cdIndex].declRef, CDstar[cdIndex].numRef);
        writeCrossEntry(stream, catName, cdName, vmag, minDistance);
    }
    return true;
}

/*
 * storeStar - almacena una estrella para futuras identificaciones
 */
int storeStar(int *count, int max, const char *name, int *ref, double *X, double *Y, double *Z,
        double *mag, int refValue, double x, double y, double z, double magValue) {
    if (*count >= max) {
        printf("Error: too many %s stars.\n", name);
        exit(1);
    }
    int i = *count;
    ref[i] = refValue;
    X[i] = x;
    Y[i] = y;
    Z[i] = z;
    if (mag != NULL) mag[i] = magValue;
    (*count)++;
    return i;
}

/*
 * checkCrossRef - revisa una referencia cruzada contra un catalogo almacenado
 */
void checkCrossRef(const char *srcName, const char *catLine, const char *label,
        double x, double y, double z, int numRefCat, const struct StarList *list,
        bool breakAfterFirst, int gcIndexRegister, int *check, int *errors) {
    for (int i = 0; i < *list->count; i++) {
        if (list->ref[i] != numRefCat) continue;
        double dist = 3600.0 * calcAngularDistance(x, y, z, list->x[i], list->y[i], list->z[i]);
        if (dist > MAX_DIST_CROSS) {
            printf("%d) Warning: %s is FAR from %s %d (dist = %.1f arcsec).\n",
                ++(*errors),
                srcName,
                label,
                numRefCat,
                dist);
            if (catLine != NULL) printf("     Register %s: %s\n", srcName, catLine);
            if (gcIndexRegister >= 0) writeRegisterGC(gcIndexRegister);
        } else (*check)++;
        if (breakAfterFirst) break;
    }
}

/*
 * checkCrossRefWB - idem contra el catalogo WB (identificacion hora + numero)
 */
void checkCrossRefWB(const char *srcName, const char *catLine, bool starWord,
        double x, double y, double z, int RARef, int numRefCat,
        const struct StarList *list, const int *hourRef, int gcIndexRegister,
        int *check, int *errors) {
    for (int i = 0; i < *list->count; i++) {
        if (hourRef[i] != RARef) continue;
        if (list->ref[i] != numRefCat) continue;
        double dist = 3600.0 * calcAngularDistance(x, y, z, list->x[i], list->y[i], list->z[i]);
        if (dist > MAX_DIST_CROSS) {
            printf("%d) Warning: %s is FAR from WB %dh %d%s (dist = %.1f arcsec).\n",
                ++(*errors),
                srcName,
                RARef,
                numRefCat,
                starWord ? " star" : "",
                dist);
            if (catLine != NULL) printf("     Register %s: %s\n", srcName, catLine);
            if (gcIndexRegister >= 0) writeRegisterGC(gcIndexRegister);
        } else (*check)++;
    }
}

/*
 * checkBDRef - idem contra el catalogo Durchmusterung (BD/CD)
 */
void checkBDRef(const char *srcName, const char *catLine, bool printNotFound,
        bool signRef, int declRef, int numRefCat,
        double x, double y, double z, int *check, int *errors) {
    struct DMstar_struct *BDstar = getDMStruct();

    int index = getDMindex(signRef, declRef, numRefCat);
    if (index == -1) {
        if (printNotFound) {
            printf("DM not found for declRef = %d, numRef = %d.\n",
                declRef,
                numRefCat);
        }
        bye("Cannot find DM star!");
    }
    double dist = 3600.0 * calcAngularDistance(x, y, z, BDstar[index].x, BDstar[index].y, BDstar[index].z);
    if (dist > MAX_DIST_CROSS) {
        printf("%d) Warning: %s is FAR from BD star (dist = %.1f arcsec).\n",
            ++(*errors),
            srcName,
            dist);
        printf("     Register %s: %s\n", srcName, catLine);
        writeRegister(index, false);
    } else (*check)++;
}

/*
 * checkCrossRefGCcat - idem contra el catalogo GC ya leido con readGC()
 */
void checkCrossRefGCcat(const char *srcName, const char *catLine, int numRefCat,
        double x, double y, double z, int *check, int *errors) {
    struct GCstar_struct *GCstar = getGCStruct();
    int GCstars = getGCStars();

    for (int i = 0; i < GCstars; i++) {
        if (GCstar[i].gcRef != numRefCat) continue;
        double dist = 3600.0 * calcAngularDistance(x, y, z, GCstar[i].x, GCstar[i].y, GCstar[i].z);
        if (dist > MAX_DIST_CROSS) {
            printf("%d) Warning: %s is FAR from GC %d (dist = %.1f arcsec). Check if it is a 1/2 star in GC.\n",
                ++(*errors),
                srcName,
                numRefCat,
                dist);
            if (catLine != NULL) printf("     Register %s: %s\n", srcName, catLine);
        } else (*check)++;
    }
}

/*
 * checkYarnallRef - revisa una referencia a Yarnall (USNO 2a edicion) buscando
 * la estrella USNO mas cercana, ya que su numeracion difiere de la de
 * Yarnall-Frisby (USNO 3a edicion)
 */
void checkYarnallRef(const char *srcName, const char *catLine, int numRefCat,
        double x, double y, double z, bool requirePositiveIndex, int gcIndexRegister,
        const struct StarList *list, int *check, int *errors) {
    double minDistance = HUGE_NUMBER;
    int usnoIndex = -1;
    for (int i = 0; i < *list->count; i++) {
        double dist = 3600.0 * calcAngularDistance(x, y, z, list->x[i], list->y[i], list->z[i]);
        if (minDistance > dist) {
            usnoIndex = i;
            minDistance = dist;
        }
    }
    if (requirePositiveIndex && usnoIndex <= 0) return;
    if (minDistance > MAX_DIST_CROSS_YARNALL) {
        printf("%d) Warning: %s is FAR from Y %d / U %d (dist = %.1f arcsec).\n",
            ++(*errors),
            srcName,
            numRefCat,
            list->ref[usnoIndex],
            minDistance);
        if (catLine != NULL) printf("     Register %s: %s\n", srcName, catLine);
        if (gcIndexRegister >= 0) writeRegisterGC(gcIndexRegister);
    } else (*check)++;
}

/*
 * warnIfAloneNorth - advierte estrella sin PPM cercana ni GSC (estilo cross_north)
 */
void warnIfAloneNorth(bool ppmFound, double minDistance, double RA, double Decl, double epoch,
        const char *warnDesc, int *errors) {
    if (ppmFound || minDistance <= MAX_DIST_PPM_FAR) return;
    if (!findGSCStar(RA, Decl, epoch, MAX_DIST_GSC)) {
        printf("%d) Warning: %s is ALONE (nearest PPM star at %.1f arcsec).\n",
            ++(*errors),
            warnDesc,
            minDistance);
    }
}

/*
 * warnAlone - advierte estrella sola (no PPM / CD / CPD / GSC) y registra causas
 */
void warnAlone(int *errors, const char *warnDesc, const char *registerDesc, const char *catLine,
        char *catName, int RAs, double decl, int Decls, int ppmRef, double nearestPPMDistance) {
    printf("%d) Warning: %s is ALONE (no PPM / CD / CPD / GSC star near it).\n",
        ++(*errors),
        warnDesc);
    if (catLine != NULL) printf("     Register %s: %s\n", registerDesc, catLine);
    logCauses(catName, true,
        false, false, RAs, decl, Decls,
        ppmRef, nearestPPMDistance);
}

/*
 * warnAlonePPMGSC - advierte estrella sola (no PPM or GSC) y registra causas
 */
void warnAlonePPMGSC(int *errors, const char *warnDesc, char *catName,
        int RAs, double decl, int Decls, int ppmRef, double nearestPPMDistance) {
    printf("%d) Warning: %s is ALONE (no PPM or GSC star near it).\n",
        ++(*errors),
        warnDesc);
    logCauses(catName, false,
        false, false, RAs, decl, Decls,
        ppmRef, nearestPPMDistance);
}

/*
 * writeUnidentified - almacena una estrella no identificada (pero hallada en GSC)
 */
void writeUnidentified(FILE *stream, const char *catName, double x, double y, double z) {
    fprintf(stream, "%s,%.8f,%.8f,%.8f\n", catName, x, y, z);
}

/*
 * printRSMEDist / printRSMEMag - imprime los RSME acumulados
 */
void printRSMEDist(const struct CrossStats *stats) {
    printf("RSME of distance (arcsec) = %.2f  among a total of %d stars\n",
        sqrt(stats->akkuDistError / (double)stats->countDist),
        stats->countDist);
}

void printRSMEMag(const struct CrossStats *stats) {
    printf("RSME of visual magnitude = %.5f  among a total of %d stars\n",
        sqrt(stats->akkuDeltaError / (double)stats->countDelta),
        stats->countDelta);
}

/*
 * checkPrecessionRA / checkPrecessionDecl - revisa la precesion anual reportada
 * contra la calculada (usamos constantes de Struve)
 */
void checkPrecessionRA(double preRA, const char *srcName, const char *catLine,
        double RA, double Decl, double base, double factor, double tol, int *errors) {
    if (fabs(preRA) <= EPS) return;
    double realPreRA = base + factor * dsin(RA) * dtan(Decl);
    double diff = fabs(preRA - realPreRA);
    if (diff > tol) {
        printf("%d) Warning: %s reports %.3f on precession RA but it should be %.3f (diff=%.3f).\n",
            ++(*errors),
            srcName,
            preRA,
            realPreRA,
            diff);
        if (catLine != NULL) printf("     Register %s: %s\n", srcName, catLine);
    }
}

void checkPrecessionDecl(double preDecl, const char *srcName, const char *catLine,
        double RA, double coef, double tol, int *errors) {
    if (fabs(preDecl) <= EPS) return;
    double realPreDecl = coef * dcos(RA);
    double diff = fabs(preDecl - realPreDecl);
    if (diff > tol) {
        printf("%d) Warning: %s reports %.2f on precession DECL but it should be %.2f (diff=%.2f).\n",
            ++(*errors),
            srcName,
            preDecl,
            realPreDecl,
            diff);
        if (catLine != NULL) printf("     Register %s: %s\n", srcName, catLine);
    }
}

/*
 * parseGCScanLine - extrae numero GC (1a columna) y referencia (3a columna)
 * de una fila escaneada de GC; false si no hay referencia confiable
 */
bool parseGCScanLine(char *buffer, int *gcRef, char *ref) {
    /* 1a columna: número GC */
    *gcRef = atoi(buffer);
    if (*gcRef <= 0) return false;

    /* ubicamos la 3a columna (luego de la 2a coma) */
    char *p = strchr(buffer, ',');
    if (p == NULL) return false;
    p = strchr(p + 1, ',');
    if (p == NULL) return false;
    p++;
    while (*p == ' ') p++; /* saltamos espacios iniciales */

    /* copiamos la referencia, sin salto de línea ni espacios finales */
    int k = 0;
    for (int i = 0; p[i] != 0 && p[i] != '\n' && p[i] != '\r' && k < 63; i++) {
        ref[k++] = p[i];
    }
    while (k > 0 && ref[k - 1] == ' ') k--;
    ref[k] = 0;

    if (ref[0] == 0) return false; /* sin referencia */
    if (strchr(ref, '(') != NULL) return false; /* parentesis: identificación no confiable */
    return true;
}

/*
 * parseSOMLine - parsea una fila de los "Standards of Magnitude" de la UA
 */
bool parseSOMLine(char *buffer, char field[17][32], char *catName, char *catLine,
        double *RA, double *Decl) {
    /* separa la linea CSV en 17 campos */
    int nf = 0, j = 0;
    for (char *p = buffer; nf < 17; p++) {
        char c = *p;
        if (c == ',' || c == '\n' || c == '\r' || c == 0) {
            field[nf][j] = 0;
            nf++;
            j = 0;
            if (c != ',') break;
        } else if (j < 31) field[nf][j++] = c;
    }
    if (nf < 17) return false; /* linea incompleta */

    /* lee numeración como string, p.ej. "34a"/"34b" para las
       dobles (omite la cabecera) */
    if (atoi(field[0]) <= 0) return false;
    snprintf(catName, 20, "SOM %s", field[0]);

    /* lee ascension recta B1875.0 */
    int RAh = atoi(field[2]);
    int RAm = atoi(field[3]);
    int RAs = atoi(field[4]);
    double ra = (double) RAh;
    ra += ((double) RAm)/60.0;
    ra += ((double) RAs)/3600.0;
    *RA = ra * 15.0; /* conversion horas a grados */

    /* lee declinacion B1875.0 (el signo viene en dec_d, p.ej. "-0") */
    int Decld = atoi(field[5]);
    if (Decld < 0) Decld = -Decld;
    double Declm = atof(field[6]);
    double decl = (double) Decld;
    decl += Declm/60.0;
    if (field[5][0] == '-') decl = -decl;
    *Decl = decl;

    snprintf(catLine, 64, "%02dh %02dm %02ds %c%02d°%04.1f'",
        RAh, RAm, RAs, field[5][0] == '-' ? '-' : '+', Decld, Declm);
    return true;
}

/*
 * parseUALine - parsea una fila de cat/ua.txt
 */
bool parseUALine(char *buffer, char *serpens, struct UAstar_struct *ua) {
    char cell[256];

    /* no leemos estrellas sin coordenadas */
    readField(buffer, cell, 101, 1);
    if (cell[0] == ' ') return false;

    /* lee ascension recta B1875.0 */
    readFieldSanitized(buffer, cell, 100, 2);
    ua->RAh = atoi(cell);
    double RA = (double) ua->RAh;
    readFieldSanitized(buffer, cell, 103, 2);
    ua->RAm = atoi(cell);
    RA += ((double) ua->RAm)/60.0;
    readFieldSanitized(buffer, cell, 106, 2);
    ua->RAs = atoi(cell);
    RA += ((double) ua->RAs)/3600.0;
    ua->RA = RA * 15.0; /* conversion horas a grados */

    /* lee declinacion B1875.0 */
    readFieldSanitized(buffer, cell, 111, 2);
    ua->Decld = atoi(cell);
    double Decl = (double) ua->Decld;
    readFieldSanitized(buffer, cell, 114, 4);
    ua->Declm = atof(cell);
    Decl += ua->Declm/60.0;
    readField(buffer, cell, 110, 1);
    if (cell[0] == '-') Decl = -Decl;
    ua->Decl = Decl;

    snprintf(ua->catLine, 64, "%02dh %02dm %02ds %c%02d°%02.1f'",
        ua->RAh, ua->RAm, ua->RAs, cell[0], ua->Decld, ua->Declm);

    /* lee numeración de Gould y constelación, si existe */
    readFieldSanitized(buffer, cell, 1, 1);
    if (cell[0] == 'G') {
        ua->existsRef = true;
        readField(buffer, cell, 3, 3);
        ua->gouldRef = atoi(cell);
        readField(buffer, cell, 7, 3);
        copyWithoutSpaces(ua->cstRef, cell);
        if (ua->cstRef[0] == 'S' && ua->cstRef[1] == 'e' && ua->cstRef[2] == 'r') {
            /* Serpens tiene parte (a) y (b) */
            ua->cstRef[3] = *serpens;
            ua->cstRef[4] = 0;
            /* si es la ultima estrella de (a), actualiza a (b) */
            if (ua->gouldRef == 49) {
                *serpens = 'b';
            }
        }
        snprintf(ua->catgName, 20, "%dG %s", ua->gouldRef, ua->cstRef);
    } else {
        ua->existsRef = false;
        ua->gouldRef = 0;
        snprintf(ua->catgName, 11, "Annonymous");
    }
    return true;
}

/*
 * splitUARefs - separa las referencias a catalogos de la UA en subceldas
 */
int splitUARefs(char *buffer, char subcell[3][18]) {
    char cell[256];

    readField(buffer, cell, 82, 17);
    int count = 0;
    int j = 0;
    for (int i = 0; i < 17; i++) {
        char c = cell[i];
        if (j == 0 && c == ' ') continue; /* skip space at beggining if exists. */
        if (c == ',') {
            subcell[count++][j] = 0;
            j = 0;
            continue;
        }
        subcell[count][j++] = c;
    }
    subcell[count++][j] = 0;
    return count;
}
