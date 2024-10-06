
/*
 * READ_CPD - Lee catalogo CPD
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "read_cpd.h"
#include "misc.h"
#include "trig.h"

#define MAX_DECL 81
#define MAX_NUM 12500

static struct CPDstar_struct CPDstar[MAXCPDSTAR];
static int CPDstars;

/* Mapa para acceder rápido según declinación y num:
   el índice devuelto es la estrella DM decl num, pero de haber
   más de una (por suppl) se utiliza la lista enlazada */
static int mapCPDindex[MAX_DECL][MAX_NUM];

/*
 * getCPDStars - devuelve la cantidad de estrellas de CPD leidas
 */
int getCPDStars()
{
    return CPDstars;
}

/*
 * getCPDindex - devuelve el índice a partir de la declinación y num
 */
int getCPDindex(int declRef, int numRef) {
    declRef = abs(declRef);
    if (declRef >= MAX_DECL) bye("Declination error!");
    if (numRef >= MAX_NUM) bye("Number error!");
    return mapCPDindex[declRef][numRef];
}

/*
 * setCPDindex - asigna el índice a partir de la declinación y num
 */
void setCPDindex(int declRef, int numRef, int index) {
    declRef = abs(declRef);
    if (declRef >= MAX_DECL) {
        printf("Num: %d %d\n", declRef, numRef);
        bye("Declination error!");
    }
    if (numRef >= MAX_NUM) bye("Number error!");
    mapCPDindex[declRef][numRef] = index;
}

/*
 * getCPDStruct - devuelve la estructura CPD
 */
struct CPDstar_struct *getCPDStruct()
{
    return &CPDstar[0];
}

/*
 * read_cpd - lee base de datos del Cape Photographic Durchmusterung cpd.txt
 * 
   Bytes Format   Units    Label       Explanations
   1-  2  A2      ---      ---         [CP] The catalog prefix
   3-  5  I3      deg      zone        [-18/-89] The declination zone
   6- 10  I5      ---      num         The number of the star within the zone
      11  A1      ---      suppl      *[a-c D] star in corrigenda
  12- 15  F4.1    mag      mag        *Estimated photographic magnitude
  16- 17  I2      h        RAh         Hours of right ascension, 1875
  18- 19  I2      min      RAm         Minutes of right ascension, 1875
  20- 23  F4.1    s        RAs        *Seconds of right ascension, 1875
      24  A1      ---      DE-         [-] Sign of declination
  25- 26  I2      deg      DEd         Degree of declination, 1875
  27- 32  F6.3    arcmin   DEm        *Minutes of declination, 1875
 *
 * Nota: se descartan objetos nebulares y suplementarias
 * (no hay suplementarias antes de la declinación -60)
 * 
 * true = utiliza catálogo 4011
 * false = utiliza catálogo 4005
 */
void readCPD(bool catalog)
{
    FILE *stream;
    char buffer[1024];
    char cell[256];

    int CDstars = getDMStars();
    struct DMstar_struct *CDstar = getDMStruct();

    stream = fopen("cat/cpd.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read cpd.txt");
        exit(1);
    }

    for (int decl = 0; decl < MAX_DECL; decl++) {
        for (int num = 0; num < MAX_NUM; num++) {
            mapCPDindex[decl][num] = -1;
        }
    }
 
    CPDstars = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
        /* lee la zona de declinacion */
        readField(buffer, cell, 3, 3);
        int declRef = atoi(cell);
        if (declRef > -22 || declRef <= -MAX_DECL) continue;

        /* lee numeracion y caracter suplementario */
        readField(buffer, cell, 6, 5);
        int numRef = atoi(cell);
        char supplRef = buffer[11-1];
        if (supplRef == 'D') continue;
        if (supplRef != ' ') {
            printf("Ommitting star CPD %d°%d%c\n", declRef, numRef, supplRef);
            continue;
        }

        /* lee magnitud fotografica */
        readField(buffer, cell, 12, 4);
        double pmag = atof(cell);
        if (pmag > 11.4) {
            /* vmag no es una magnitud, si no un codigo */
            if ((pmag > 19.9 && pmag < 20.1) || (pmag > 29.9 && pmag < 30.1)) continue;
            printf("Unknown code: %f, CPD %d°%d\n", pmag, declRef, numRef);
            exit(1);
        }

        /* lee ascension recta B1875.0 */
        readField(buffer, cell, 16, 2);
        float rah = atof(cell);
        double RA = rah;
        readField(buffer, cell, 18, 2);
        float ramin = atof(cell);
        RA += ramin/60.0;
        readField(buffer, cell, 20, 4);
        float raseg = atof(cell);
        RA += raseg/3600.0;
        RA *= 15.0; /* conversion horas a grados */

        /* lee declinacion B1875.0 */
        readField(buffer, cell, 25, 2);
        float decldeg = atof(cell);
        double Decl = decldeg;
        readField(buffer, cell, 27, 6);
        float declmin = atof(cell);
        Decl += declmin/60.0;
        Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RA, Decl, &x, &y, &z);

        /* la almacena en memoria */
		if (CPDstars == MAXCPDSTAR) bye("Maximum amount reached!\n");
        CPDstar[CPDstars].discard = true; /* todas comienzan descartadas a menos que se cruce */
        CPDstar[CPDstars].declRef = declRef;
        CPDstar[CPDstars].numRef = numRef;
        CPDstar[CPDstars].RA1875 = RA;
        CPDstar[CPDstars].Decl1875 = Decl;
        CPDstar[CPDstars].pmag = pmag;
        CPDstar[CPDstars].dmIndex = -1; /* a ser rellenado en la siguiente fase */
        CPDstar[CPDstars].x = x;
        CPDstar[CPDstars].y = y;
        CPDstar[CPDstars].z = z;
        CPDstar[CPDstars].dist = 0.0; /* a ser rellenado en la siguiente fase */
        CPDstar[CPDstars].declRef2 = -1; /* idem */
        CPDstar[CPDstars].numRef2 = -1; /* idem */

        setCPDindex(declRef, numRef, CPDstars);

        CPDstars++;
    }
    printf("Stars read from Cape Photographic Durchmusterung: %d\n", CPDstars);
    fclose(stream);

    if (catalog) {
        /* siguiente fase: leer identificación cruzada del catálogo 4011 (Bonnet) */
        stream = fopen("cat/4011.txt", "rt");
        if (stream == NULL) {
            perror("Cannot read 4011.txt");
            exit(1);
        }
    } else {
        /* siguiente fase: leer identificación cruzada del catálogo 4005 */
        stream = fopen("cat/4005.txt", "rt");
        if (stream == NULL) {
            perror("Cannot read 4005.txt");
            exit(1);
        }
    }

    int crossed = 0;
    int numRefCDprevious = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
        int declRefCP, numRefCP, declRefCD, numRefCD;
        if (catalog) { /* CATALOGO 4011 */
            /* lee la zona de declinacion de CPD y numero */
            readField(buffer, cell, 11, 3);
            declRefCP = atoi(cell);
            if (declRefCP > -22) continue;
            readField(buffer, cell, 15, 5);
            numRefCP = atoi(cell);
            if (numRefCP == 0) continue;
            readField(buffer, cell, 14, 1);
            if (cell[0] != ' ') {
                printf("Something different from space in col 14 for CPD %d°%d\n", declRefCP, numRefCP);
                exit(1);
            }

            /* lee la estrella de CD asociada, si existe */
            readField(buffer, cell, 23, 3);
            if (cell[1] == '.') {
                //printf("CPD %d°%d with no associated CD star. Discarding it.\n", declRefCP, numRefCP);
                continue;
            }
            declRefCD = atoi(cell);
            readField(buffer, cell, 27, 5);
            numRefCD = atoi(cell);
            if (numRefCD == 0) continue;
            readField(buffer, cell, 26, 1);
            if (cell[0] != ' ') {
                printf("Something different from space in col 26 for CD %d°%d\n", declRefCD, numRefCD);
                exit(1);
            }
        } else { /* CATALOGO 4005 */
            /* necesitamos que una sea CPD y la otra sea CD */
            readField(buffer, cell, 12, 1);
            bool isSourceCP = cell[0] == '4';
            bool isSourceCD = cell[0] == '2';
            readField(buffer, cell, 25, 1);
            bool isTargetCP = cell[0] == '4';
            bool isTargetCD = cell[0] == '2';
            int declRefSource, numRefSource, declRefTarget, numRefTarget;
            if (isSourceCP && isTargetCD) {
                declRefSource = 2;
                numRefSource = 6;
                declRefTarget = 15;
                numRefTarget = 19;
            } else if (isSourceCD && isTargetCP) {
                declRefSource = 15;
                numRefSource = 19;
                declRefTarget = 2;
                numRefTarget = 6;
            } else continue;

            /* lee la zona de declinacion de CPD y numero */
            readField(buffer, cell, declRefSource, 3);
            declRefCP = atoi(cell);
            readField(buffer, cell, numRefSource, 5);
            numRefCP = atoi(cell);
            if (numRefCP == 0) bye("Error in CPD num!");

            /* lee la estrella de CD asociada, si existe */
            readField(buffer, cell, declRefTarget, 3);
            declRefCD = atoi(cell);
            readField(buffer, cell, numRefTarget, 5);
            numRefCD = atoi(cell);
            if (numRefCD == 0) bye("Error in CD num!");
        }

        /* buscar la estrella en el catalogo CPD */
        int catIndex = getCPDindex(declRefCP, numRefCP);
        if (catIndex == -1) continue;

        /* obtiene posición en la esfera unidad */
        double x = CPDstar[catIndex].x;
        double y = CPDstar[catIndex].y;
        double z = CPDstar[catIndex].z;

        /* Busca la estrella asociada en DM; en caso de haber más
         * de una, escoger la de menor distancia */
        int dmIndex = -1;
        double minDistance = 9999999999;
        int i = getDMindex(true, declRefCD, numRefCD);
        while (i != -1) {
          double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
          if (minDistance > dist) {
            dmIndex = i;
            minDistance = dist;
          }
          i = CDstar[i].next;
        }
        if (dmIndex == -1) {
          printf("CD %d°%d not found (corresponding to CPD %d°%d). Discarding CPD star.\n",
            declRefCD, numRefCD, declRefCP, numRefCP);
          continue;
        }

        if (catalog) {
          /* chequeo para ver si hay "continuidad" en los números */  
          if (declRefCD == declRefCP) {
            if (abs(numRefCD - numRefCDprevious) > 50) {
              if (declRefCD == declRefCP && numRefCP > 2) {
                printf("Warning: CPD %d°%d asoociated to CD %d, but previous one was CD %d.\n",
                    declRefCP, numRefCP, numRefCD, numRefCDprevious);
              } 
            }
            numRefCDprevious = numRefCD;
          }
        }

        /* se almacena la lectura de la identificacion */
        CPDstar[catIndex].discard = false;
        CPDstar[catIndex].dmIndex = dmIndex;
        CPDstar[catIndex].dist = minDistance;
        CDstar[dmIndex].catIndex = catIndex;
        crossed++;
    }
    printf("Number of CPD stars cross-identified with CD stars: %d\n", crossed);
    fclose(stream);

    /*  leer identificación cruzada del catálogo 4019 (Rappaport) */
    stream = fopen("cat/4019.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read 4019.txt");
        exit(1);
    }
    while (fgets(buffer, 1023, stream) != NULL) {
        /* lee la zona de declinacion de CPD y numero */
        readField(buffer, cell, 12, 3);
        int declRefCP = atoi(cell);
        if (declRefCP <= -MAX_DECL) continue;
        readField(buffer, cell, 15, 5);
        int numRefCP = atoi(cell);
        int index = getCPDindex(declRefCP, numRefCP);
        if (index == -1) continue;
        if (CPDstar[index].discard) continue;

        /* guarda la zona de declinacion de CD y numero */
        readField(buffer, cell, 1, 3);
        int declRefCD = atoi(cell);
        readField(buffer, cell, 4, 5);
        int numRefCD = atoi(cell);

        int dmIndex = CPDstar[index].dmIndex;
        if (declRefCD != CDstar[dmIndex].declRef || numRefCD != CDstar[dmIndex].numRef) {
            /* las diferencias son almacenadas */
            printf("Warning: cross-identifications differ for CPD %d°%d: CD %d°%d vs %d°%d\n",
                declRefCP, numRefCP, declRefCD, numRefCD, CDstar[dmIndex].declRef, CDstar[dmIndex].numRef);
            CPDstar[index].declRef2 = declRefCD;
            CPDstar[index].numRef2 = numRefCD;
        }
    }
    fclose(stream);
}
