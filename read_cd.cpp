
/*
 * READ_CD - Lee catálogo CD y lo ubica en memoria
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include "read_dm.h"
#include "misc.h"
#include "trig.h"

#define MAX_DECL 90
#define MAX_NUM 20000

static struct DMstar_struct CDstar[MAXDMSTAR];
static int CDstars;
static int CDstarsTomo16, firstIndexTomo16;
static int CDstarsTomo17, firstIndexTomo17;
static int CDstarsTomo18, firstIndexTomo18;
static int CDstarsTomo21a, firstIndexTomo21a;
static int CDstarsTomo21b, firstIndexTomo21b;

/* Mapa para acceder rápido según declinación y num:
   el índice devuelto es la estrella DM decl num, pero de haber
   más de una (por suppl) se utiliza la lista enlazada */
static int mapCDindex[MAX_DECL][MAX_NUM];

/*
 * getDMStars - devuelve la cantidad de estrellas de CD leidas
 */
int getDMStars()
{
    return CDstars;
}

/*
 * getDMindex - devuelve el índice a partir del signo, declinación y num
 */
int getDMindex(bool signRef, int declRef, int numRef) {
    if (!signRef) bye("Sign error!");
    declRef = abs(declRef);
    if (declRef >= MAX_DECL) bye("Declination error!");
    if (numRef >= MAX_NUM) bye("Number error!");
    return mapCDindex[declRef][numRef];
}

/*
 * setDMindex - asigna el índice a partir del signo, declinación y num
 */
void setDMindex(bool signRef, int declRef, int numRef, int index) {
    if (!signRef) bye("Sign error!");
    declRef = abs(declRef);
    if (declRef >= MAX_DECL) bye("Declination error!");
    if (numRef >= MAX_NUM) bye("Number error!");
    mapCDindex[declRef][numRef] = index;
}

/*
 * isCD - dice que es el catálogo CD
 */
bool isCD()
{
    return true;
}

/*
 * getDMStruct - devuelve la estructura CD
 */
struct DMstar_struct *getDMStruct()
{
    return &CDstar[0];
}

/*
 * writeRegister - escribe en pantalla un registro de CD en formato ONA
 */
void writeRegister(int cdIndex)
{
    struct DMstar_struct *s = &CDstar[cdIndex];

    int declRef = s->declRef;
    int numRef = s->numRef;
    int page = 0;
    if (declRef <= -22 && declRef >= -31) {
        page = 605 * (cdIndex - firstIndexTomo16) / CDstarsTomo16;
    }
    if (declRef <= -32 && declRef >= -41) {
        page = 539 * (cdIndex - firstIndexTomo17) / CDstarsTomo17;
    }
    if (declRef <= -42 && declRef >= -51) {
        page = 503 * (cdIndex - firstIndexTomo18) / CDstarsTomo18;
    }
    if (declRef <= -52 && declRef >= -61) {
        page = 305 * (cdIndex - firstIndexTomo21a) / CDstarsTomo21a;
    }
    if (declRef <= -62) {
        // TO DO: obtener número de páginas
        page = 400 * (cdIndex - firstIndexTomo21b) / CDstarsTomo21b;
    }
    printf("     Register CD %d°%d (en %.0fh):  %.1f | %.0fm%.1fs | %.1f'     (pag. %d)\n",
                declRef,
                numRef,
                s->rah,
                s->vmag,
                s->ramin,
                s->raseg,
                s->declmin,
                page);
}

/*
 * read_cd - lee base de datos del Durchmusterung
   Bytes Format   Units    Label       Explanations
    1-  2  A2      ---      ---         [CD] The catalog prefix
    3-  5  I3      deg      zone        [-22/-89] The declination zone
    6- 10  I5      ---      num         The number of the star within the zone
       11  A1      ---      suppl      *[a-c D] star in corrigenda
   12- 15  F4.1    mag      mag        *Estimated visual magnitude
   16- 17  I2      h        RAh         Hours of right ascension, 1875
   18- 19  I2      min      RAm         Minutes of right ascension, 1875
   20- 23  F4.1    s        RAs        *Seconds of right ascension, 1875
       24  A1      ---      DE-         [-] Sign of declination
   25- 26  I2      deg      DEd         Degree of declination, 1875
   27- 30  F4.1    arcmin   DEm         Minutes of declination, 1875
 *
 * Nota: se descartan estrellas variables y objetos nebulares
 */
void readDM(const char *filename)
{
    FILE *stream;
    char buffer[1024];
    char cell[256];

    stream = fopen(filename, "rt");
    if (stream == NULL) {
        snprintf(buffer, 1024, "Cannot read %s", filename);
        perror(buffer);
        exit(1);
    }

    for (int decl = 0; decl < MAX_DECL; decl++) {
        for (int num = 0; num < MAX_NUM; num++) {
            mapCDindex[decl][num] = -1;
        }
    }
 
    CDstars = 0;
    CDstarsTomo16 = 0;
    CDstarsTomo17 = 0;
    CDstarsTomo18 = 0;
    CDstarsTomo21a = 0;
    CDstarsTomo21b = 0;
    firstIndexTomo16 = -1;
    firstIndexTomo17 = -1;
    firstIndexTomo18 = -1;
    firstIndexTomo21a = -1;
    firstIndexTomo21b = -1;
    while (fgets(buffer, 1023, stream) != NULL) {
        /* lee la zona de declinacion */
        readField(buffer, cell, 3, 3);
        int declRef = atoi(cell);

        /* lee numeracion y caracter suplementario */
        readField(buffer, cell, 6, 5);
        int numRef = atoi(cell);
        readField(buffer, cell, 11, 1);
        char supplRef = cell[0];
        if (supplRef == 'D') continue;

        /* lee magnitud visual */
        readField(buffer, cell, 12, 4);
        double vmag = atof(cell);
        if (vmag > 12.1) {
            /* vmag no es una magnitud, si no un codigo:
             * 20.0 = neb
             * 30.0 = var
             * 40.0 = M4 
             * 50.0 = 47 Tuc */
            if ((vmag > 19.9 && vmag < 20.1) || (vmag > 39.9 && vmag < 40.1) || (vmag > 49.9 && vmag < 50.1)) continue;
            if (vmag > 29.9 && vmag < 30.1) {
                /* estrella variable */
            } else {
                printf("Unknown code: %f for CD %d°%d. Setting as variable.\n", vmag, declRef, numRef);
                vmag = 30.0;
                continue;
            }
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
        readField(buffer, cell, 27, 4);
        float declmin = atof(cell);
        Decl += declmin/60.0;
        Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RA, Decl, &x, &y, &z);

        /* la almacena en memoria */
		if (CDstars == MAXDMSTAR) bye("Maximum amount reached!\n");
        CDstar[CDstars].signRef = true;
        CDstar[CDstars].declRef = declRef;
        CDstar[CDstars].numRef = numRef;
        CDstar[CDstars].supplRef = supplRef;
        CDstar[CDstars].rah = rah;
        CDstar[CDstars].ramin = ramin;
        CDstar[CDstars].raseg = raseg;
        CDstar[CDstars].decldeg = decldeg;
        CDstar[CDstars].declmin = declmin;
        CDstar[CDstars].RA1875 = RA;
        CDstar[CDstars].Decl1875 = Decl;
        CDstar[CDstars].vmag = vmag;
        CDstar[CDstars].catIndex = -1; // se rellenará luego
        CDstar[CDstars].x = x;
        CDstar[CDstars].y = y;
        CDstar[CDstars].z = z;
        CDstar[CDstars].next = -1;

        /* si ya hay otra de misma identificación, la enlaza */
        int index = getDMindex(true, declRef, numRef);
        if (index == -1) {
            setDMindex(true, declRef, numRef, CDstars);
        } else {
            while (CDstar[index].next != -1) index = CDstar[index].next;
            CDstar[index].next = CDstars;
        } 

        /* crea metadata para conocer el número de página en el Tomo */
        if (declRef <= -22 && declRef >= -31) {
            if (firstIndexTomo16 == -1) firstIndexTomo16 = CDstars;
            CDstarsTomo16++;
        }
        if (declRef <= -32 && declRef >= -41) {
            if (firstIndexTomo17 == -1) firstIndexTomo17 = CDstars;
            CDstarsTomo17++;
        }
        if (declRef <= -42 && declRef >= -51) {
            if (firstIndexTomo18 == -1) firstIndexTomo18 = CDstars;
            CDstarsTomo18++;
        }
        if (declRef <= -52 && declRef >= -61) {
            if (firstIndexTomo21a == -1) firstIndexTomo21a = CDstars;
            CDstarsTomo21a++;
        }
        if (declRef <= -62) {
            if (firstIndexTomo21b == -1) firstIndexTomo21b = CDstars;
            CDstarsTomo21b++;
        }

        /* proxima estrella */
        CDstars++;
    }
    printf("Stars read from Cordoba Durchmusterung: %d\n", CDstars);
    printf("   Tomo XVI: %d, Tomo XVII: %d, Tomo XVIII: %d, Tomo XXIa: %d, Tomo XXIb: %d\n",
                CDstarsTomo16, CDstarsTomo17, CDstarsTomo18, CDstarsTomo21a, CDstarsTomo21b);
    fclose(stream);
}
