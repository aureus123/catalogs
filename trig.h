
/*
 * TRIG - Header
 */

#define PI 3.1415926535897932384
#define COSDISTTOL  0.999998476913288 /* cos(6 arcmin) */
#define COSDISTDPLTOL 0.999999924785823 /* cos(80 arcsec) */
#define HUGE_NUMBER 9999999999
#define EPS 1E-8

void transform(double eq1, double eq2, double *RA, double *Decl);
double dcos(double angle);
double dsin(double angle);
double dtan(double angle);
double datan2(double y, double x);
void sph2rec(double ra, double decl, double *x, double *y, double *z);
void rec2sph(double x, double y, double z, double *ra, double *decl);
double calcCosDistance(double x1, double y1, double z1, double x2, double y2, double z2);
double calcAngularDistance(double x1, double y1, double z1, double x2, double y2, double z2);
bool solve3x3(double A[3][3], double b[3], double x[3]);
double compVmagToCDmag(int decl_ref, double vmag);
