
/*
 * MISC - Funciones misceláneas
 * Made in 2024 by Daniel E. Severin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "misc.h"

/*
 * bye - muestra un error y aborta
 */
void bye(const char *string)
{
	printf("%s", string);
	exit(1);
}

/*
 * readField - lee un campo de una longitud determinada en "bytes" (donde "initial" comienza a partir de 1)
 */
void readField(char *buffer, char *cell, int initial, int bytes)
{
    int length = strlen(buffer);
    if (initial > length) {
        /* solo lee "ceros" (en caso de que la entrada sea mas corta que lo que se desea leer */
        for (int i = 0; i < bytes; i++) cell[i] = 0;
    }
    else {
        /* lee el campo */
        for (int i = 0; i < bytes; i++) cell[i] = buffer[i+initial-1];
    }
    cell[bytes] = 0;
}

/*
 * readFieldSanitized - lee un campo de una longitud determinada en "bytes"
 * En caso de encontrar espacios o "&", los reemplaza por ceros.
 */
void readFieldSanitized(char *buffer, char *cell, int initial, int bytes) {
    readField(buffer, cell, initial, bytes);
    for (int i = 0; i < bytes; i++) {
        if (cell[i] == ' ' || cell[i] == '&') cell[i] = '0';
    }
}

/*
 * openCrossFile - abre un archivo de identificación cruzada
 */
FILE *openCrossFile(const char *name)
{
    FILE *stream = fopen(name, "wt");
    if (stream == NULL) {
        perror("Cannot write in cross file");
        exit(1);
    }
    fprintf(stream, "index1,index2,dist\n");
    return stream;
}

/*
 * writeCrossEntry - escribe una entrada en un archivo de identificación cruzada
 */
void writeCrossEntry(FILE *stream, char *index1, char *index2, double dist)
{
    fprintf(stream, "%s,%s,%.2f\n", index1, index2, dist);
}

/*
 * openPositionFile - abre un archivo de posiciones
 */
FILE *openPositionFile(const char *name)
{
    FILE *stream = fopen(name, "wt");
    if (stream == NULL) {
        perror("Cannot write in position file");
        exit(1);
    }
    fprintf(stream, "index,decl,num,ref,dist10\n");
    return stream;
}

/*
 * writePositionEntry - escribe una entrada en un archivo de posiciones
 */
void writePositionEntry(FILE *stream, int index, int decl, int num, const char *ref, double dist)
{
    fprintf(stream, "%d,%d,%d,%s,%.0f\n", index, decl, num, ref, 10.0 * dist);
}

/*
 * openMagnitudeFile - abre un archivo de magnitudes
 */
FILE *openMagnitudeFile(const char *name)
{
    FILE *stream = fopen(name, "wt");
    if (stream == NULL) {
        perror("Cannot write in magnitude file");
        exit(1);
    }
    fprintf(stream, "index,decl,num,ref,delta10\n");
    return stream;
}

/*
 * writeMagnitudeEntry - escribe una entrada en un archivo de magnitudes
 */
void writeMagnitudeEntry(FILE *stream, int index, int decl, int num, const char *ref, double delta)
{
    fprintf(stream, "%d,%d,%d,%s,%.0f\n", index, decl, num, ref, 10.0 * delta);
}

/*
 * openUnidentifiedFile - abre un archivo para estrellas no identificadas
 */
FILE *openUnidentifiedFile(const char *name)
{
    FILE *stream = fopen(name, "wt");
    if (stream == NULL) {
        perror("Cannot write in unidentified file");
        exit(1);
    }
    fprintf(stream, "name,x,y,z\n");
    return stream;
}

/*
 * logCauses - escribe posibles causas de falta de identificacion
 * También, en caso que stream != null, almacena la estrella en un archivo, junto
 * con sus coordenadas rectangulares en 1875.0, siempre que se den ciertas causas
 * razonables para identificarla con una estrella débil.
 */
void logCauses(char *name, FILE *stream, double x, double y, double z,
        bool cumulus, bool nebula, double vmag,
        int RAs, double Decl, int Decls,
        int ppmRef, double nearestPPMDistance) {
    bool store = true;
    if (cumulus) {
        printf("  Possible cause: cumulus.\n");
        store = false;
    }
    if (nebula) {
        printf("  Possible cause: nebula.\n");
        store = false;
    }
    if (RAs == 0) {
        printf("  Possible cause: lack of RA (s).\n");
        store = false;
    }
    if (Decls == 0) {
        printf("  Possible cause: lack of Decl (s).\n");
        store = false;
    }
    if (vmag >= 8.0) {
        printf("  Possible cause: dim star.\n");
    } else store = false;
    if (vmag < 0.1) {
        printf("  Possible cause: no magnitude (cumulus?).\n");
        store = false;
    }
    if (ppmRef != -1) {
        printf("  Note: Nearest PPM %d at %.1f arcsec.\n", ppmRef, nearestPPMDistance);
    }
    if (Decl > -18.0) {
        printf("  Note: no CD/CPD coverage for stars below 18°.\n");
    }
    if (Decl < -61.0) {
        printf("  Note: poor CD coverage for stars above 61°.\n");
    }
    if (store && stream != nullptr) {
        fprintf(stream, "%s,%.12f,%.12f,%.12f\n", name, x, y, z);
    }
}

/*
 * copy - Copia sin espacios
 */
void copyWithoutSpaces(char *dest, char *src) {
	bool previousIsSpace = true;
	int srcPtr = 0;
	int destPtr = 0;
	while (src[srcPtr] != 0) {
		if (src[srcPtr] != ' ') {
			dest[destPtr++] = src[srcPtr];
			previousIsSpace = false;
		} else {
			if (!previousIsSpace) {
				// omit further consecutive spaces
				previousIsSpace = true;
				dest[destPtr++] = src[srcPtr];
			}
		}
		srcPtr++;
	}
	dest[destPtr] = 0;
}

