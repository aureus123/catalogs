
/*
 * TRIG - Funciones trigonométricas y matematicas en gral.
 * Made in 2024 by Daniel E. Severin
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "trig.h"

/* Para uso de la libreria WCS: */
#define WCS_B1950 2 /* B1950(FK4) right ascension and declination */
extern "C" void wcsconp(int sys1, int sys2, double eq1, double eq2, double ep1, double ep2,
             double *dtheta, double *dphi, double *ptheta, double *pphi);

/* Constantes para makeDoubles */
static const double MAX_DIST_DOUBLE = 60.0;     // arcsec
static const double MIN_DIST_NODOUBLE = 300.0;  // arcsec
static const double MAG_MIN_DOUBLE = 3.0;
static const double MAG_MAX_DOUBLE = 8.0;


/*
 * Transforma coordenadas entre épocas (ep1 y ep2) en FK4
 * No se consideran movimientos propios.
 * 
 * eq1, eq2 - Epoca origen y destino
 * RA, Decl - Ascensión recta y declinación (ambos, en grados)
 */
void transform(double eq1, double eq2, double *RA, double *Decl) {
    double pmRA = 0.0;
    double pmDecl = 0.0;
    wcsconp(WCS_B1950, WCS_B1950, eq1, eq2, eq1, eq2, RA, Decl, &pmRA, &pmDecl);
}

/*
 * Funciones trigonométricas en grados sexagesimales
 */
double dcos(double angle)
{
    return cos(angle * PI / 180.0);
}

double dsin(double angle)
{
    return sin(angle * PI / 180.0);
}

double dtan(double angle)
{
    return tan(angle * PI / 180.0);
}

double datan2(double y, double x)
{
    if (x < 0.000001 && x > -0.000001) {
        /* caso x = 0 */
        if (y > 0.0) return 90.0;
        else return 270.0;
    }
    double angle = atan(y/x) * 180.0 / PI;
    if (x < 0.0) {
        if (y < 0.0) angle -= 180.0;
        else angle += 180.0;
    }
    if (angle < 0.0) angle += 360.0;
    return angle;
}

void sph2rec(double ra, double decl, double *x, double *y, double *z)
{
    *x = dcos(ra) * dcos(decl);
    *y = dsin(ra) * dcos(decl);
    *z = dsin(decl);
}

void rec2sph(double x, double y, double z, double *ra, double *decl)
{
    *ra = datan2(y, x);
    *decl = datan2(z, sqrt(x*x + y*y));
    if (*decl > 90.0) *decl -= 360.0;
}

/*
 * calcCosDist - calcula el producto escalar entre 2 coordenadas de la esfera unidad
 * puede ir de -1 (180 grados de distancia) hasta 1 (0 grados de distancia)
 */
double calcCosDistance(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return x1*x2 + y1*y2 + z1*z2;
}

/*
 * calcAngularDist - calcula la distancia angular entre 2 coordenadas de la esfera unidad
 * puede ir de 0 a 180 grados
 */
double calcAngularDistance(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double cosdist = calcCosDistance(x1, y1, z1, x2, y2, z2);
	if (cosdist < -1.0) cosdist = -1.0;
	if (cosdist > 1.0) cosdist = 1.0;
	return acos(cosdist) * 180.0 / PI;
}

/*
 * solve3x3 - small 3x3 linear system solver with partial pivoting.
 * Returns true on success, false on singular.
 */
bool solve3x3(double A[3][3], double b[3], double x[3]) {
    int p[3] = {0, 1, 2};
    // Partial pivoting for column 0
    int maxr = 0;
    if (fabs(A[1][0]) > fabs(A[maxr][0])) maxr = 1;
    if (fabs(A[2][0]) > fabs(A[maxr][0])) maxr = 2;
    if (maxr != 0) {
        for (int j = 0; j < 3; ++j) { double t = A[0][j]; A[0][j] = A[maxr][j]; A[maxr][j] = t; }
        double tb = b[0]; b[0] = b[maxr]; b[maxr] = tb;
    }
    if (fabs(A[0][0]) < 1e-15) return false;
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
    if (fabs(A[1][1]) < 1e-15) return false;
    double m21 = A[2][1] / A[1][1];
    A[2][2] -= m21 * A[1][2];
    b[2] -= m21 * b[1];

    if (fabs(A[2][2]) < 1e-15) return false;

    // Back substitution
    x[2] = b[2] / A[2][2];
    x[1] = (b[1] - A[1][2] * x[2]) / A[1][1];
    x[0] = (b[0] - A[0][1] * x[1] - A[0][2] * x[2]) / A[0][0];
    return true;
}

/*
 * compVmagToCDmag - dada una magnitud en Johnson V la convierte a la escala usada en CD
 */
double compVmagToCDmag(int decl_ref, double vmag)
{
    if (decl_ref <= -22 && decl_ref >= -31) {
        // Quadratic fit of 11659 stars, RSME = 0.3004
        return -0.157169 + 1.188316*vmag -0.022130*vmag*vmag;
    }
    if (decl_ref <= -32 && decl_ref >= -41) {
        // Quadratic fit of 11796 stars, RSME = 0.3013 
        return -1.517044 + 1.595675*vmag -0.050674*vmag*vmag;
    }
    if (decl_ref <= -42 && decl_ref >= -51) {
        // Quadratic fit of 9620 stars, RSME = 0.2714
        return -4.903298 + 2.435433*vmag -0.098302*vmag*vmag;
    }
    if (decl_ref <= -52 && decl_ref >= -61) {
        // Quadratic fit of 7672 stars, RSME = 0.3834
        return -1.347719 + 1.513742*vmag -0.040125*vmag*vmag;
    }
    if (decl_ref <= -62) {
        // Quadratic fit of 7372 stars, RSME = 0.3239
        return -4.814060 + 2.342925*vmag -0.087990*vmag*vmag;
    }

    // If we reach here, BD scale should be used
    // Quadratic fit of 1412 stars, RSME = 0.1155
    return 0.441147 + 1.699235*vmag -0.079943*vmag*vmag;

    // old CD fit for full sky:
    // return -0.046197294300226*vmag*vmag + 1.435555317168187*vmag - 0.557294231481489;
}

/**
 * compCDmagToVmag - dada una magnitud en la escala usada en CD la convierte a Johnson V
 * se obtuvo con la herramienta:
 * > ./cross_txt --csv results/cross/cross_cd_vol1_ppm.csv
 * > ./cross_txt --csv results/cross/cross_cd_vol2_ppm.csv
 * > ./cross_txt --csv results/cross/cross_cd_vol3_ppm.csv
 * > ./cross_txt --csv results/cross/cross_cd_vol4_ppm.csv
 * > ./cross_txt --csv results/cross/cross_cd_vol5_ppm.csv
 */
double compCDmagToVmag(int decl_ref, double cdVmag) {
    if (decl_ref <= -22 && decl_ref >= -31) {
        // Weighted quadratic fit of 15187 stars, RSME = 0.075, MAPE = 0.78%
        return -8.955 + 3.221 * cdVmag - 0.132 * cdVmag * cdVmag;
    }
    if (decl_ref <= -32 && decl_ref >= -41) {
        // Weighted quadratic fit of 15038 stars, RSME = 0.085, MAPE = 0.94%
        return -4.874 + 2.084 * cdVmag - 0.055 * cdVmag * cdVmag;
    }
    if (decl_ref <= -42 && decl_ref >= -51) {
        // Weighted linear fit of 25597 stars, RSME = 0.092, MAPE = 1.03%
        return -1.302 + 1.163 * cdVmag;
    }
    if (decl_ref <= -52 && decl_ref >= -61) {
        // Weighted quadratic fit of 12980 stars, RSME = 0.149, MAPE = 1.34%
        return -5.647 + 2.384 * cdVmag - 0.083 * cdVmag * cdVmag;
    }
    if (decl_ref <= -62) {
        // Weighted linear fit of 12218 stars, RSME = 0.052, MAPE = 0.55%
        return 0.154 + 0.981 * cdVmag;
    }

    abort();
    return 0;
}

/*
 * compGCmagToVmag - dada una magnitud en la escala usada en GC la convierte a Johnson V
*/
double compGCmagToVmag(double gcVmag) {
    // Quadratic fit of 28977 stars grouped in 48 magnitudes, ECM = 0.202
    return - 5.003 + 2.305 * gcVmag - 0.085 * gcVmag * gcVmag;
}

/*
 * makeDoubles - busca pares de estrellas dobles aisladas y los guarda en un CSV.
 *
 * Una estrella i forma una doble con j si:
 *   - mag[i] y mag[j] estan en [3, 8]
 *   - dist(i, j) < MAX_DIST_DOUBLE
 *   - dentro de MIN_DIST_NODOUBLE alrededor de i no hay otras estrellas que j
 *     (y simetricamente, alrededor de j no hay otras que i)
 * Cada par se emite una sola vez, con la estrella mas brillante primero.
 * Se reporta la ascension recta (en horas) y declinacion (en grados enteros)
 * de la estrella mas brillante, habitualmente en la época 1875, junto con los
 * nombres, magnitudes y distancia entre ambas estrellas.
 *
 * Algoritmo O(n^2): se hace una pasada para contar, por cada estrella, cuantos
 * vecinos tiene dentro de MIN_DIST_NODOUBLE y cual es uno de ellos. Luego se
 * recorren las estrellas y se emiten las que tienen exactamente un vecino y
 * cuyo vecino tambien tiene exactamente uno.
 */
void makeDoubles(int n, int *ref, double *X, double *Y, double *Z,
                 double *mag, const char *desig, const char *filename)
{
    int *nearCount = (int*) calloc(n, sizeof(int));
    int *nearIdx = (int*) malloc(n * sizeof(int));
    double *nearDist = (double*) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        nearIdx[i] = -1;
        nearDist[i] = HUGE_NUMBER;
    }

    /* Pasada O(n^2): contamos vecinos dentro de MIN_DIST_NODOUBLE para cada estrella. */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d = 3600.0 * calcAngularDistance(X[i], Y[i], Z[i],
                                                    X[j], Y[j], Z[j]);
            if (d < MIN_DIST_NODOUBLE) {
                if (nearCount[i] == 0) {
                    nearIdx[i] = j;
                    nearDist[i] = d;
                }
                nearCount[i]++;
                if (nearCount[j] == 0) {
                    nearIdx[j] = i;
                    nearDist[j] = d;
                }
                nearCount[j]++;
            }
        }
    }

    FILE *stream = fopen(filename, "wt");
    if (stream == NULL) {
        perror("Cannot write doubles file");
        exit(1);
    }
    fprintf(stream, "RAh,Decl,name1,mag1,name2,mag2,dist\n");

    int count = 0;
    for (int i = 0; i < n; i++) {
        if (nearCount[i] != 1) continue;
        int j = nearIdx[i];
        if (i > j) continue;                  /* emit each pair only once */
        if (nearCount[j] != 1) continue;      /* simetria del aislamiento */
        if (nearDist[i] >= MAX_DIST_DOUBLE) continue;
        if (mag[i] < MAG_MIN_DOUBLE || mag[i] > MAG_MAX_DOUBLE) continue;
        if (mag[j] < MAG_MIN_DOUBLE || mag[j] > MAG_MAX_DOUBLE) continue;

        int a = i, b = j;
        if (mag[b] < mag[a]) { a = j; b = i; }   /* brighter first */

        double ra, decl;
        rec2sph(X[a], Y[a], Z[a], &ra, &decl);
        if (ra < 0.0) ra += 360.0;
        int rah = (int) round(ra / 15.0);
        if (rah == 24) rah = 0;
        int decl_int = (int) round(decl);

        fprintf(stream, "%d,%d,%s %d,%.1f,%s %d,%.1f,%.1f\n",
                rah, decl_int,
                desig, ref[a], mag[a],
                desig, ref[b], mag[b],
                nearDist[i]);
        count++;
    }

    fclose(stream);
    free(nearCount);
    free(nearIdx);
    free(nearDist);

    printf("Doubles found in %s catalog: %d (saved to %s)\n", desig, count, filename);
}
