
/*
 * READ_BD - Lee catálogo BD y lo ubica en memoria
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include "read_dm.h"
#include "misc.h"
#include "trig.h"

#define MAX_DECL 20
#define MAX_NUM 7000

static struct DMstar_struct BDstar[MAXDMSTAR];
static int BDstars;
static int reindexByIndex[MAXDMSTAR];

/* Mapa para acceder rápido según signo, declinación y num:
   primer argumento: 0 = positive, 1 = negative
   el índice devuelto es la estrella DM decl num, pero de haber
   más de una (por suppl) se utiliza la lista enlazada */
static int mapBDindex[2][MAX_DECL][MAX_NUM];

/*
 * getBDStars - devuelve la cantidad de estrellas de BD leidas
 */
int getDMStars()
{
    return BDstars;
}

/*
 * getDMindex - devuelve el índice a partir del signo, declinación y num
 */
int getDMindex(bool signRef, int declRef, int numRef) {
    declRef = abs(declRef);
    if (declRef >= MAX_DECL) bye("Declination error!");
    if (numRef >= MAX_NUM) bye("Number error!");
    return mapBDindex[signRef ? 1 : 0][declRef][numRef];
}

/*
 * setDMindex - asigna el índice a partir del signo, declinación y num
 */
void setDMindex(bool signRef, int declRef, int numRef, int index) {
    declRef = abs(declRef);
    if (declRef >= MAX_DECL) bye("Declination error!");
    if (numRef >= MAX_NUM) bye("Number error!");
    mapBDindex[signRef ? 1 : 0][declRef][numRef] = index;
}

/*
 * isCD - dice que es el catálogo BD
 */
bool isCD()
{
    return false;
}

/*
 * getDMStruct - devuelve la estructura CD
 */
struct DMstar_struct *getDMStruct()
{
    return &BDstar[0];
}

/*
 * writeRegister - escribe en pantalla un registro de BD en su formato
 */
void writeRegister(int bdIndex, bool neighbors)
{
    struct DMstar_struct *s = &BDstar[bdIndex];

    int declRef = s->declRef;
    int numRef = s->numRef;
    int page = 1 + 378 * reindexByIndex[bdIndex] / BDstars;
    printf("     Register BD %d°%d (en %.0fh):  %.1f | %.0fm%.1fs | %.1f'     (pag. %d)\n",
                declRef,
                numRef,
                s->rah,
                s->vmag,
                s->ramin,
                s->raseg,
                s->declmin,
                page);
    if (neighbors) {
        printf("   Neighbors:\n");
        writeRegister(bdIndex - 1, false);
        writeRegister(bdIndex + 1, false);
    }
}

/*
 * read_bd - lee base de datos del Durchmusterung
   Bytes      Format   Units    Label       Explanations

    1- 2        A2      ---      ---        [BD] The catalog prefix
       3        A1      ---      zonesign   [+-]The sign of the declination zone
    4- 5        I2      deg      zone       The declination zone
    6-10        I5      ---      num        The number of the star within the zone
      11        A1      ---      suppl     *[a-c D*?M] Note
   12-15      F4.1      mag      mag       *Estimated visual magnitude
   16-17        I2      h        RAh        Right Ascension 1855 (hours)
   18-19        I2      min      RAm        Right Ascension 1855 (minutes)
   20-23      F4.1      s        RAs        Right Ascension, 1855 (seconds)
   24-24        A1      ---      DE-        [+-]Sign of declination
   25-26        I2      deg      DEd        Declination 1855 (degrees)
   27-30      F4.1      arcmin   DEm        Declination 1855 (minutes)
 *
 * Nota: se descartan estrellas variables y objetos nebulares
 * Solo se lee primer tomo (declinaciones -1 a 19)
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

    for (int sign = 0; sign < 2; sign++) {
        for (int decl = 0; decl < MAX_DECL; decl++) {
            for (int num = 0; num < MAX_NUM; num++) {
                mapBDindex[sign][decl][num] = -1;
            }
        }
    }
 
    BDstars = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
        /* lee la zona de declinacion (es necesario también conocer el signo para
         * diferenciar la declinación +00 de la -00). */
        readField(buffer, cell, 3, 1);
        bool zoneSign = (cell[0] == '-');
        readField(buffer, cell, 4, 2);
        int declRef = atoi(cell);
        if (zoneSign) declRef = -declRef;
        if (declRef >= 20) continue;

        /* lee numeracion y caracter suplementario */
        readField(buffer, cell, 6, 5);
        int numRef = atoi(cell);
        readField(buffer, cell, 11, 1);
        char supplRef = cell[0];
        if (supplRef == 'D') continue;
        if (supplRef == '*') {
            printf("Star already corrected (BD %d°%d)\n", declRef, numRef);
            continue;
        }
        
        /* lee magnitud visual */
        readField(buffer, cell, 12, 4);
        double vmag = atof(cell);
        if (vmag > 12.1) {
            /* vmag no es una magnitud, si no un codigo:
             * 20.0 = neb
             * 30.0 = var
             * 40.0 = nova or nova?
             * 50.0 = cluster */
            if ((vmag > 19.9 && vmag < 20.1) || (vmag > 39.9 && vmag < 40.1) || (vmag > 49.9 && vmag < 50.1)) continue;
            if (vmag > 29.9 && vmag < 30.1) {
                /* estrella variable */
            } else {
                printf("Unknown code: %f, BD %d°%d\n", vmag, declRef, numRef);
                exit(1);
            }
        }

        /* lee ascension recta B1855.0 */
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

        /* lee declinacion B1855.0 */
        readField(buffer, cell, 24, 1);
        bool sign = (cell[0] == '-');
        readField(buffer, cell, 25, 2);
        float decldeg = atof(cell);
        double Decl = decldeg;
        readField(buffer, cell, 27, 4);
        float declmin = atof(cell);
        Decl += declmin/60.0;
        if (sign) Decl = -Decl; /* incorpora signo negativo solo si es necesario */

        /* calcula coordenadas rectangulares */
        double x, y, z;
        sph2rec(RA, Decl, &x, &y, &z);

        /* la almacena en memoria */
		if (BDstars == MAXDMSTAR) bye("Maximum amount reached!\n");
        BDstar[BDstars].signRef = zoneSign;
        BDstar[BDstars].declRef = declRef;
        BDstar[BDstars].numRef = numRef;
        BDstar[BDstars].supplRef = supplRef;
        BDstar[BDstars].rah = rah;
        BDstar[BDstars].ramin = ramin;
        BDstar[BDstars].raseg = raseg;
        BDstar[BDstars].decldeg = decldeg;
        BDstar[BDstars].declmin = declmin;
        BDstar[BDstars].RA1855 = RA;
        BDstar[BDstars].Decl1855 = Decl;
        BDstar[BDstars].vmag = vmag;
        BDstar[BDstars].catIndex = -1; // se rellenará luego
        BDstar[BDstars].x = x;
        BDstar[BDstars].y = y;
        BDstar[BDstars].z = z;
        BDstar[BDstars].next = -1;

        /* si ya hay otra de misma identificación, la enlaza */
        int index = getDMindex(zoneSign, declRef, numRef);
        if (index == -1) {
            setDMindex(zoneSign, declRef, numRef, BDstars);
        } else {
            while (BDstar[index].next != -1) index = BDstar[index].next;
            BDstar[index].next = BDstars;
        }

        /* proxima estrella */
        BDstars++;
    }
    printf("Stars read from Bonner Durchmusterung: %d\n", BDstars);

    /* Ahora calculamos el indice pero con declinaciones ascendentes, aquí decl = -2 es declinación -01, decl = 1 es -00
       y para decl >= 0 ya es la declinación positiva. */
    int reindex = 0;
    for (int decl = -2; decl <= 19; decl++) {
        for (int i = 0; i < BDstars; i++) {
            if (decl == -2) {
                if (BDstar[i].declRef != -1 || !BDstar[i].signRef) continue;
            } else {
                if (decl == -1) {
                    if (BDstar[i].declRef != 0 || !BDstar[i].signRef) continue;
                } else {
                    if (BDstar[i].declRef != decl || BDstar[i].signRef) continue;
                }
            }
            reindexByIndex[i] = reindex;
            reindex++;
        }
    }
    if (reindex != BDstars) {
        printf("Some star is missing or duplicated :(\n");
        exit(1);
    }
    fclose(stream);
}


/* 
 * findDMByCoordinates - busca la estrella BD más cercana
 * Aquí (x, y, z) son las coord rectangulares en 1855.
 * minDistanceOutput debe ser una cota de la distancia a buscar.
 * El resultado se almacena en (bdIndexOutput, minDistanceOutput).
*/
void findDMByCoordinates(double x, double y, double z, int *bdIndexOutput, double *minDistanceOutput) {
    int bdIndex = -1;
    double minDistance = *minDistanceOutput;
	for (int i = 0; i < BDstars; i++) {
        double dist = 3600.0 * calcAngularDistance(x, y, z, BDstar[i].x, BDstar[i].y, BDstar[i].z);
        if (minDistance > dist) {
        	bdIndex = i;
        	minDistance = dist;
        }
    }
    *bdIndexOutput = bdIndex;
    *minDistanceOutput = minDistance;
}
