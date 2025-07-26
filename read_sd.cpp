/*
 * READ_SD - Lee catalogo SD
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_dm.h"
#include "read_sd.h"
#include "trig.h"
#include "misc.h"

static struct SDstar_struct SDstar[MAXSDSTAR];
static int SDstars = 0;

/*
 * getSDStars - devuelve la cantidad de estrellas de SD leidas
 */
int getSDStars()
{
    return SDstars;
}

/*
 * getSDStruct - devuelve la estructura SD
 */
struct SDstar_struct *getSDStruct()
{
    return &SDstar[0];
}

/*
 * writeRegisterSD - escribe en pantalla un registro de SD en su formato
 * (solo usar en caso de leer únicamente la declinación -22, de otro modo
 * la página reportada no es correcta)
 */
void writeRegisterSD(int sdIndex)
{
    struct SDstar_struct *s = &SDstar[sdIndex];

    int numRef = s->numRef;
    int page = 438 + 22 * sdIndex / SDstars;
    printf("     Register SD %d°%d (en %.0fh):  %.1f | %.0fm%.1fs | %.1f'     (pag. %d)\n",
                s->declRef,
                numRef,
                s->rah,
                s->vmag,
                s->ramin,
                s->raseg,
                s->declmin,
                page);
}

/*
 * readSD - lee estrellas del Southern Durchmusterung
 * (mismo formato que BD, ver read_bd.cpp).
 * onlyDecl22: si es true, solo lee decl -22 y hace identificacion cruzada con CD.
 *             si es false, lee todas las declinaciones disponibles.
 */
void readSD(bool onlyDecl22) {
    FILE *stream;
    char buffer[1024];
    char cell[256];

    int CDstars = 0;
    struct DMstar_struct *CDstar;
    
    if (onlyDecl22) {
        CDstars = getDMStars();
        CDstar = getDMStruct();
    }

    stream = fopen("cat/sd.txt", "rt");
    if (stream == NULL) {
        perror("Cannot read sd.txt");
        exit(1);
    }

    SDstars = 0;
    while (fgets(buffer, 1023, stream) != NULL) {
        /* lee la zona de declinacion */
        readField(buffer, cell, 3, 3);
        int declRef = atoi(cell);
        if (onlyDecl22 && declRef != -22) continue;

        /* lee numeracion y caracter suplementario */
        readField(buffer, cell, 6, 5);
        int numRef = atoi(cell);
        char supplRef = buffer[11-1];
        if (supplRef == 'D') continue;
        if (supplRef != ' ') {
            printf("Ommitting star SD %d°%d%c\n", declRef, numRef, supplRef);
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
                printf("Unknown code: %f, SD %d°%d\n", vmag, declRef, numRef);
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
        readField(buffer, cell, 25, 2);
        float decldeg = atof(cell);
        double Decl = decldeg;
        readField(buffer, cell, 27, 6);
        float declmin = atof(cell);
        Decl += declmin/60.0;
        Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */

        /* calcula coordenadas rectangulares, pero de 1875.0 */
        double RA1875 = RA;
        double Decl1875 = Decl;
        transform(1855.0, 1875.0, &RA1875, &Decl1875);
        double x, y, z;
        sph2rec(RA1875, Decl1875, &x, &y, &z);

        /* la almacena en memoria */
		if (SDstars == MAXSDSTAR) bye("Maximum amount reached!\n");
        SDstar[SDstars].discard = true; /* todas comienzan descartadas a menos que se cruce */
        SDstar[SDstars].declRef = declRef;
        SDstar[SDstars].numRef = numRef;
        SDstar[SDstars].rah = rah;
        SDstar[SDstars].ramin = ramin;
        SDstar[SDstars].raseg = raseg;
        SDstar[SDstars].decldeg = decldeg;
        SDstar[SDstars].declmin = declmin;
        SDstar[SDstars].RA1855 = RA;
        SDstar[SDstars].Decl1855 = Decl;
        SDstar[SDstars].RA1875 = RA1875;
        SDstar[SDstars].Decl1875 = Decl1875;
        SDstar[SDstars].vmag = vmag;
        SDstar[SDstars].cdIndex = -1; /* a ser rellenado en la siguiente fase */
        SDstar[SDstars].x = x;
        SDstar[SDstars].y = y;
        SDstar[SDstars].z = z;
        SDstar[SDstars].dist = 0.0; /* a ser rellenado en la siguiente fase */
        SDstars++;
    }
    printf("Stars read from Southern Durchmusterung: %d\n", SDstars);
    fclose(stream);

    if (onlyDecl22) {
        /* siguiente fase: leer identificación cruzada del catálogo 4005 */
        stream = fopen("cat/4005.txt", "rt");
        if (stream == NULL) {
            perror("Cannot read 4005.txt");
            exit(1);
        }

        int crossed = 0;
        while (fgets(buffer, 1023, stream) != NULL) {
            /* necesitamos que una sea SD y la otra sea CD */
            readField(buffer, cell, 12, 1);
            bool isSourceSD = cell[0] == ' ';
            bool isSourceCD = cell[0] == '2';
            readField(buffer, cell, 25, 1);
            bool isTargetSD = cell[0] == ' ';
            bool isTargetCD = cell[0] == '2';
            int declRefSource, numRefSource, declRefTarget, numRefTarget;
            if (isSourceSD && isTargetCD) {
                declRefSource = 2;
                numRefSource = 6;
                declRefTarget = 15;
                numRefTarget = 19;
            } else if (isSourceCD && isTargetSD) {
                declRefSource = 15;
                numRefSource = 19;
                declRefTarget = 2;
                numRefTarget = 6;
            } else continue;

            /* lee la zona de declinacion de SD y numero */
            readField(buffer, cell, declRefSource, 3);
            int declRefSD = atoi(cell);
            if (declRefSD != -22) continue;
            readField(buffer, cell, numRefSource, 5);
            int numRefSD = atoi(cell);
            if (numRefSD == 0) bye("Error in SD num!");

            /* lee la estrella de CD asociada, si existe */
            readField(buffer, cell, declRefTarget, 3);
            int declRefCD = atoi(cell);
            readField(buffer, cell, numRefTarget, 5);
            int numRefCD = atoi(cell);
            if (numRefCD == 0) bye("Error in CD num!");

            /* buscar la estrella en el catalogo SD */
            int catIndex = -1;
            for (int i = 0; i < SDstars; i++) {
                if (numRefSD == SDstar[i].numRef) {
                    catIndex = i;
                    break;
                }
            }
            if (catIndex == -1) {
                printf("Warning: Star SD %d°%d not found\n", declRefSD, numRefSD);
                continue;
            }

            /* obtiene posición en la esfera unidad */
            double x = SDstar[catIndex].x;
            double y = SDstar[catIndex].y;
            double z = SDstar[catIndex].z;

            /* Busca la estrella asociada en CD; en caso de haber más
            * de una, escoger la de menor distancia */
            int cdIndex = -1;
            double minDistance = 9999999999;
            int i = getDMindex(true, declRefCD, numRefCD);
            while (i != -1) {
                double dist = 3600.0 * calcAngularDistance(x, y, z, CDstar[i].x, CDstar[i].y, CDstar[i].z);
                if (minDistance > dist) {
                    cdIndex = i;
                    minDistance = dist;
                }
                i = CDstar[i].next;
            }
            if (cdIndex == -1) {
            printf("CD %d°%d not found (corresponding to SD %d°%d). Discarding SD star.\n",
                declRefCD, numRefCD, declRefSD, numRefSD);
            continue;
            }

            /* se almacena la lectura de la identificacion */
            SDstar[catIndex].discard = false;
            SDstar[catIndex].cdIndex = cdIndex;
            SDstar[catIndex].dist = minDistance;
            CDstar[cdIndex].catIndex = catIndex;
            crossed++;
        }
        printf("Number of SD stars cross-identified with CD stars: %d\n", crossed);
        fclose(stream);
    }
}

/* 
 * findSDByCoordinates - busca la estrella SD más cercana
 * Aquí (x, y, z) son las coord rectangulares en 1875.
 * (aquí "decl" aun no está implementado y no se utiliza)
 * minDistanceOutput debe ser una cota de la distancia a buscar.
 * El resultado se almacena en (sdIndexOutput, minDistanceOutput).
*/
void findSDByCoordinates(double x, double y, double z, double decl, int *sdIndexOutput, double *minDistanceOutput) {
    int sdIndex = -1;
    double minDistance = *minDistanceOutput;
    for (int i = 0; i < SDstars; i++) {
        double dist = 3600.0 * calcAngularDistance(x, y, z, SDstar[i].x, SDstar[i].y, SDstar[i].z);
        if (minDistance > dist) {
            sdIndex = i;
            minDistance = dist;
        }
    }
    *sdIndexOutput = sdIndex;
    *minDistanceOutput = minDistance;
} 