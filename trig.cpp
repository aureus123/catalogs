
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
 * compVmagToCDmag - dada una magnitud en Johnson V la convierte a la escala usada en CD
 */
double compVmagToCDmag(int decl_ref, double vmag)
{
    if (decl_ref <= -22 && decl_ref >= -31) {
        // Quadratic fit of 11659 stars, ECM = 0.3004
        return -0.157169 + 1.188316*vmag -0.022130*vmag*vmag;
    }
    if (decl_ref <= -32 && decl_ref >= -41) {
        // Quadratic fit of 11796 stars, ECM = 0.3013 
        return -1.517044 + 1.595675*vmag -0.050674*vmag*vmag;
    }
    if (decl_ref <= -42 && decl_ref >= -51) {
        // Quadratic fit of 9620 stars, ECM = 0.2714
        return -4.903298 + 2.435433*vmag -0.098302*vmag*vmag;
    }
    if (decl_ref <= -52 && decl_ref >= -61) {
        // Quadratic fit of 7672 stars, ECM = 0.3834
        return -1.347719 + 1.513742*vmag -0.040125*vmag*vmag;
    }
    if (decl_ref <= -62) {
        // Quadratic fit of 7372 stars, ECM = 0.3239
        return -4.814060 + 2.342925*vmag -0.087990*vmag*vmag;
    }

    // If we reach here, BD scale should be used
    // Quadratic fit of 1412 stars, ECM = 0.1155
    return 0.441147 + 1.699235*vmag -0.079943*vmag*vmag;

    // old CD fit:
    // return -0.046197294300226*vmag*vmag + 1.435555317168187*vmag - 0.557294231481489;
}
