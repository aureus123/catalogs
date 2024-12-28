
/*
 * READ_PPM - Lee catálogo PPM
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "read_ppm.h"
#include "misc.h"
#include "trig.h"

/* Para uso de la libreria WCS: */
#define WCS_J2000 1 /* J2000(FK5) right ascension and declination */
#define WCS_B1950 2 /* B1950(FK4) right ascension and declination */
extern "C" void wcsconp(int sys1, int sys2, double eq1, double eq2, double ep1, double ep2,
             double *dtheta, double *dphi, double *ptheta, double *pphi);
//SYNOPSIS:   wcsconp(sys1, sys2, eq1, eq2, ep1, ep2, dtheta, dphi, ptheta, pphi)
//  int     sys1;   /* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC) */
//  int     sys2;   /* Output coordinate system (J2000, B1950, ECLIPTIC, GALACTIC) */
//  double  eq1;    /* Input equinox (default of sys1 if 0.0) */
//  double  eq2;    /* Output equinox (default of sys2 if 0.0) */
//  double  ep1;    /* Input Besselian epoch in years (for proper motion) */
//  double  ep2;    /* Output Besselian epoch in years (for proper motion) */
//  double  *dtheta; /* Right ascension in degrees: Input in sys1, returned in sys2 */
//  double  *dphi;  /* Declination in degrees: Input in sys1, returned in sys2 */
//  double  *ptheta; /* Right ascension proper motion in RA degrees/year: Input in sys1, returned in sys2 */
//  double  *pphi;  /* Declination proper motion in Dec degrees/year: Input in sys1, returned in sys2 */

static struct PPMstar_struct PPMstar[MAXPPMSTAR];
static int PPMstars;

/*
 * getPPMstars - devuelve la cantidad de estrellas de PPM leidas
 */
int getPPMStars()
{
    return PPMstars;
}

/*
 * getPPMstruct - devuelve la estructura PPM
 */
struct PPMstar_struct *getPPMStruct()
{
    return &PPMstar[0];
}

/*
 * revise - revisa si mas de una estrella PPM se condice con una de DM
 */
bool revise(int ppmIndex) {
    int dmIndex = PPMstar[ppmIndex].dmIndex;
    bool warning = false;
    for (int i = 0; i < PPMstars; i++) {
      if (i == ppmIndex) continue;
      if (PPMstar[i].dmIndex == dmIndex) {
        if (PPMstar[i].discard) {
          printf("     Note: Also associated to discarded PPM %d (dist = %.1f arcsec).\n", PPMstar[i].ppmRef, PPMstar[i].dist);
        } else {
          printf("     Warning: Also associated to PPM %d (dist = %.1f arcsec).\n", PPMstar[i].ppmRef, PPMstar[i].dist);
          warning = true;
        }
      }
    }
    return warning;
}

/*
 * read_ppm - lee base de datos PPM
 *
   Bytes  Format  Units   Label    Explanations
   2-  7   I6     ---     PPM     *[181732/378910]+ Designation of the star
  10- 18   A9     ---     DM      *Durchmusterung, BD or CD (10-12 decl., 13-17 number, 18 component [AB])
  20- 23   F4.1   mag     Mag     *Magnitude, Visual if Flag5 is 'V'
  25- 26   A2     ---     Sp      *Spectral type
  28- 29   I2     h       RAh      Right Ascension for the Equinox=J2000.0 and
                                   Epoch=J2000.0, on the system of FK5
  31- 32   I2     min     RAm      Right Ascension J2000 (minutes)
  34- 39   F6.3   s       RAs      Right Ascension J2000 (seconds)
      42   A1     ---     DE-      Declination J2000 (sign)
  43- 44   I2     deg     DEd      Declination for the equinox and epoch
                                   J2000.0, on the system of FK5
  46- 47   I2     arcmin  DEm      Declination J2000 (minutes)
  49- 53   F5.2   arcsec  DEs      Declination J2000 (seconds)
  56- 62   F7.4   s/yr    pmRA    *Proper motion in RA, J2000
  64- 69   F6.3 arcsec/yr pmDE    *Proper motion in DE, J2000
  71- 72   I2     ---     Npos    *Number of individual positions used
  74- 75   I2     10mas   e_RA    *Mean error of RA
  77- 78   I2     10mas   e_DE    *Mean error of DE
  80- 83   F4.1   mas/yr  e_pmRA  *Mean error of pmRA
  85- 88   F4.1   mas/yr  e_pmDE  *Mean error of pmDE
  90- 94   F5.2   yr    EpRA-1900 *Weighted mean epoch, RA and pmRA
  96-100   F5.2   yr    EpDE-1900 *Weighted mean epoch, DE and pmDE
 102-107   I6     ---     SAO      [1/258997]? SAO Designation
 109-114   I6     ---     HD       [1/359083]? Henry Draper Designation
 117-125   A9     ---     CPD      Cape Photographic Durchmusterung Designation
     127   A1     ---     Flag1   *'P' - problem, 'C' - comment
     128   A1     ---     Flag2   *'D' - double star
     129   A1     ---     Flag3   *    - not used
     130   A1     ---     Flag4   *'F' - member of FK5
     131   A1     ---     Flag5   *'R' - remark, 'V' - V Mag (CPC-2)
--------------------------------------------------------------------------------
Note on PPM:
  Designation of the star in PPM South starting with No. 181732,
  see chapter 3 of the author's description file "desc.txt".
Note on DM:
  Designation of the star in the Bonner Durchmusterung (for zones from
  -02 to -22 degrees) or in the Cordoba Durchmusterung (zones -23 to
  -89 degrees).
Note on Mag and Sp:
  See chapter 3 of description file "desc.txt".
Note on pmRA:
  Proper motion in right ascension for epoch and equinox J2000.0,
  on the system of FK5 given in seconds of time per Julian year.
Note on pmDE:
  Proper motion in Declination for epoch and equinox J2000.0,
  on the system of FK5 given in seconds of arc per Julian year.
Note on Npos:
  Number of individual published positions used for the
  derivation of the position and proper motion given.
Note on e_RA:
  Mean error of right ascension at the mean epoch of right ascension EpRA,
  multiplied by the cosine of declination, given in units of 0.01
  seconds of arc.
Note on e_DE:
  Mean error of declination at the mean epoch of declination EpDE, given
  in units of 0.01 seconds of arc.
Note on e_pmRA:
  Mean error of the proper motion in right ascension, multiplied by
  the cosine of declination, given in units of 0.001 seconds of arc
  per Julian year.
Note on e_pmDE:
  mean error of the proper motion in declination given in units of
  0.001 seconds of arc per Julian year
Note on EpRA-1900:
  Weighted mean epoch of the measured positions used for the
  derivation of RA and pmRA, given in years since 1900.0.
Note on EpDE-1900:
  Weighted mean epoch of the measured positions used for the
  derivation of DE and pmDE, given in years since 1900.0.
Note on Flag1:
  P - problem case, preferably not to be used as astrometric reference star.
  C - a critical comment is given in the List of Critical Comments.
      Not to be used as astrometric reference star.
Note on Flag2:
  D - double star, preferably not to be used as astrometric reference star.
Note on Flag3:
  Not used, void for consistency with PPM (north).
Note on Flag4:
  F - member of FK5, mostly bright stars, original FK5 data are given.
Note on Flag5:
  R - a remark is given in the List of Remarks on Individual Stars.
  V - the magnitude is a photographic V magnitude copied from CPC-2.
 *
 * si useDurch = true, lee la identificacion cruzada con Durchmusterung y
 *   usa flag allSky, sino lee todo el catálogo PPM sin la identificación a BD/CD
 * si allSky = true, de -22 al polo sur; si es false, hasta -31 inclusive
 */
void readPPM(bool useDurch, bool allSky, double targetYear)
{
    FILE *stream;
    char buffer[1024];
    char cell[256];

    int DMstars = getDMStars();
    struct DMstar_struct *DMstar = getDMStruct();

    stream = fopen("cat/ppm.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read ppm.txt");
        exit(1);
    }

    PPMstars = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
      bool zoneSign;
      int declRefAbs, numRef, declRef;
      if (useDurch) {
        readField(buffer, cell, 10, 1);
        zoneSign = (cell[0] == '-');
        readField(buffer, cell, 11, 2);
        declRefAbs = atoi(cell);

        if (isCD()) {
          /* para CD: de -23 hasta el polo sur, excepto que allSky=false en cuyo es el 1er. volumen */
          if (!zoneSign) continue;
          if (declRefAbs < 23) continue;
          if (!allSky && declRefAbs > 31) continue;
        } else {
          /* para BD: entre -01 y +19 */
          if (declRefAbs > 19 && !zoneSign) continue;
          if (declRefAbs > 1 && zoneSign) continue;
        }
        declRef = zoneSign ? -declRefAbs : declRefAbs;

        readField(buffer, cell, 13, 5);
        numRef = atoi(cell);
        if (numRef == 0) continue;
      }

      /* lee identificacion PPM */
      readField(buffer, cell, 2, 6);
      int ppmRef = atoi(cell);

      /* lee ascension recta J2000 */
      readField(buffer, cell, 28, 2);
      double RA = atof(cell);
      readField(buffer, cell, 31, 2);
      RA += atof(cell)/60.0;
      readField(buffer, cell, 34, 6);
      RA += atof(cell)/3600.0;
      RA *= 15.0; /* conversion horas a grados */

      /* lee declinacion J2000 */
      readField(buffer, cell, 42, 1);
      char sign = cell[0];
      readField(buffer, cell, 43, 2);
      double Decl = atof(cell);
      readField(buffer, cell, 46, 2);
      Decl += atof(cell)/60.0;
      readField(buffer, cell, 49, 5);
      Decl += atof(cell)/3600.0;
      if (sign == '-') Decl = -Decl; /* incorpora signo negativo en caso de ser necesario */

      /* lee mov. propio en asc. recta */
      readField(buffer, cell, 56, 7);
      double pmRA = atof(cell);
      pmRA /= 240; /* conversion de s/yr a grados/yr (juliano) */
      
      /* lee mov. propio en declinacion */
      readField(buffer, cell, 64, 6);
      double pmDecl = atof(cell);
      pmDecl /= 3600; /* conversion de arcsec/yr a grados/yr (juliano) */

      /* convierte coordenadas al Siglo XIX: como hay que introducir la epoca besseliana de J2000, usamos la formula de SOFA:
        * B = 1900 + (2451545 - 2415020.31352) / 365.242198781 = 2000.001278 */
      double RAtarget = RA;
      double Decltarget = Decl;
      wcsconp(WCS_J2000, WCS_B1950, 0.0, targetYear, 2000.001278, targetYear, &RAtarget, &Decltarget, &pmRA, &pmDecl);

      /* calcula coordenadas rectangulares (segun el catálogo destino) */
      double x, y, z;
      sph2rec(RAtarget, Decltarget, &x, &y, &z);

      /* lee magnitud (si Flag5 != 'V' la magnitud es fotografica o hay una remark --> poner 0.0) */
      double vmag = 0.0;
      readField(buffer, cell, 131, 1);
      if (cell[0] == 'V') {
          readField(buffer, cell, 20, 4);
          vmag = atof(cell);
      }

      /* lee si la estrella es "problematica" */
      char problem = 0;
      readField(buffer, cell, 127, 1);
      if (cell[0] == 'P' || cell[0] == 'C') problem = 1;
      readField(buffer, cell, 128, 1);
      if (cell[0] == 'D') problem = 1;

      int dmIndex = -1;
      double minDistance = HUGE_NUMBER;
      if (useDurch) {
        /* lee el numero DM y busca la estrella asociada en DM
        * en caso de haber más de una, escoger la de menor distancia */
        char dm[20];
        snprintf(dm, 20, "DM %s%d°%d", zoneSign ? "-" : "+", declRefAbs, numRef);
        int i = getDMindex(zoneSign, declRef, numRef);
        while (i != -1) {
          double dist = 3600.0 * calcAngularDistance(x, y, z, DMstar[i].x, DMstar[i].y, DMstar[i].z);
          if (minDistance > dist) {
            dmIndex = i;
            minDistance = dist;
          }
          i = DMstar[i].next;
        }
        if (dmIndex == -1) {
          printf("Star %s not found (corresponding to PPM %d). Discarding PPM star.\n", dm, ppmRef);
          continue;
        }

        /* asocia la DM a la PPM mas cercana */
        if (DMstar[dmIndex].catIndex == -1) {
          DMstar[dmIndex].catIndex = PPMstars;
        } else {
          /* ya hay otra PPM con misma DM asociada */
          int previousPPMIndex = DMstar[dmIndex].catIndex;
          if (minDistance < PPMstar[previousPPMIndex].dist) {
            printf("PPM %d is removed because PPM %d is nearer to same %s\n", PPMstar[previousPPMIndex].ppmRef, ppmRef, dm);
            PPMstar[previousPPMIndex].discard = true;
            DMstar[dmIndex].catIndex = PPMstars;
          } else {
            printf("PPM %d is removed because PPM %d is nearer to same %s\n", ppmRef, PPMstar[previousPPMIndex].ppmRef, dm);
            continue;
          }
        }
      }

      /* la almacena en memoria */
      if (PPMstars == MAXPPMSTAR) bye("Maximum amount reached!\n");
      PPMstar[PPMstars].ppmRef = ppmRef;
      PPMstar[PPMstars].discard = false;
      PPMstar[PPMstars].vmag = vmag;
      PPMstar[PPMstars].problem = problem;
      PPMstar[PPMstars].dmIndex = dmIndex;
      PPMstar[PPMstars].x = x;
      PPMstar[PPMstars].y = y;
      PPMstar[PPMstars].z = z;
      PPMstar[PPMstars].dist = minDistance;

      /* proxima estrella */
      PPMstars++;
    }
    printf("Stars read from PPM: %d\n", PPMstars);
    fclose(stream);
}

/* 
 * findPPMByCoordinates - busca la estrella PPM más cercana
 * Aquí (x, y, z) son las coord rectangulares en el año target.
 * minDistanceOutput debe ser una cota de la distancia a buscar.
 * El resultado se almacena en (ppmIndexOutput, minDistanceOutput).
*/
void findPPMByCoordinates(double x, double y, double z, int *ppmIndexOutput, double *minDistanceOutput) {
  int ppmIndex = -1;
  double minDistance = *minDistanceOutput;
	for (int i = 0; i < PPMstars; i++) {
        double dist = 3600.0 * calcAngularDistance(x, y, z, PPMstar[i].x, PPMstar[i].y, PPMstar[i].z);
        if (minDistance > dist) {
        	ppmIndex = i;
        	minDistance = dist;
        }
    }
    *ppmIndexOutput = ppmIndex;
    *minDistanceOutput = minDistance;
}
