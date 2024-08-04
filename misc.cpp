
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
