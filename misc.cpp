
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
 * En caso de encontrar espacios, los reemplaza por ceros.
 */
void readFieldSanitized(char *buffer, char *cell, int initial, int bytes) {
    readField(buffer, cell, initial, bytes);
    for (int i = 0; i < bytes; i++) {
        if (cell[i] == ' ') cell[i] = '0';
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
