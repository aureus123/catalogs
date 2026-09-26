
/*
 * READ_GC - Lee catálogo general argentino
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "trig.h"
#include "misc.h"
#include "read_gc.h"

struct GCstar_struct GCstar[MAXGCSTAR];

int GCstars;

/*
 * getGCStars - devuelve la cantidad de estrellas de GC leidas
 */
int getGCStars()
{
    return GCstars;
}

/*
 * getGCStruct - devuelve la estructura GC
 */
struct GCstar_struct *getGCStruct()
{
    return &GCstar[0];
}

/*
 * writeRegister - escribe en pantalla un registro de GC
 */
void writeRegisterGC(int index) {
	printf("     Register GC %d: mag = %.1f, RA = %02dh%02dm%02ds%02d, DE = %02d°%02d'%02d''%01d (pag %d) %s%s%s\n",
		GCstar[index].gcRef,
		GCstar[index].vmag,
		GCstar[index].RAh,
		GCstar[index].RAm,
		GCstar[index].RAs / 100,
		GCstar[index].RAs % 100,
		GCstar[index].Decld,
		GCstar[index].Declm,
		GCstar[index].Decls / 10,
		GCstar[index].Decls % 10,
		GCstar[index].page,
		GCstar[index].dpl ? "dpl " : "",
		GCstar[index].cum ? "cum " : "",
		GCstar[index].neb ? "neb " : ""
	);
}

/*
 * Lee coordenadas rectangulares de una estrella GC
 * Si devuelve true, la estrella fue encontrada y se cargaron las coordenadas e índice
 */
bool getGCStarData(int gcRef, int *index, double *x, double *y, double *z)
{
    bool found = false;
    for (int i = 0; i < GCstars; i++) {
        if (GCstar[i].gcRef != gcRef) continue;
        *x = GCstar[i].x;
        *y = GCstar[i].y;
        *z = GCstar[i].z;
        *index = i;
        found = true;
        break;
    }
	return found;
}

/*
 * hasBlanks - devuelve true si el campo contiene algun espacio (dato ausente o incompleto)
 */
static bool hasBlanks(char *buffer, int initial, int bytes)
{
	char cell[256];
	readField(buffer, cell, initial, bytes);
	for (int i = 0; i < bytes; i++) {
		if (cell[i] == ' ' || cell[i] == 0) return true;
	}
	return false;
}

/*
 * readRA - lee ascension recta B1875.0 (en grados)
 * Devuelve true si los segundos estan completos (sin espacios)
 */
static bool readRA(char *buffer, int *RAh, int *RAm, int *RAs, double *RA)
{
	char cell[256];
	readFieldSanitized(buffer, cell, 16, 2);
	*RAh = atoi(cell);
	*RA = (double) *RAh;
	readFieldSanitized(buffer, cell, 18, 2);
	*RAm = atoi(cell);
	*RA += ((double) *RAm)/60.0;
	readFieldSanitized(buffer, cell, 20, 4);
	*RAs = atoi(cell);
	*RA += (((double) *RAs)/100.0)/3600.0;
	*RA *= 15.0; /* conversion horas a grados */
	return !hasBlanks(buffer, 20, 4);
}

/*
 * readDecl - lee declinacion B1875.0 (en grados)
 * Devuelve true si los segundos estan completos (sin espacios)
 */
static bool readDecl(char *buffer, int *Decld, int *Declm, int *Decls, double *Decl)
{
	char cell[256];
	readFieldSanitized(buffer, cell, 39, 2);
	*Decld = atoi(cell);
	*Decl = (double) *Decld;
	readFieldSanitized(buffer, cell, 41, 2);
	*Declm = atoi(cell);
	*Decl += ((double) *Declm)/60.0;
	readFieldSanitized(buffer, cell, 43, 3);
	*Decls = atoi(cell);
	*Decl += (((double) *Decls)/10.0)/3600.0;
	*Decl = -*Decl; /* incorpora signo negativo (en nuestro caso, siempre) */
	return !hasBlanks(buffer, 43, 3);
}

/*
 * Lee estrellas del Primer Catalogo Argentino, coordenadas 1875.0
 * supuestamente todas estas estrellas deberian estar incluidas en el catálogo CD
 * (excepto las que están fuera de la faja, y algunas de CD marcadas como "dobles")
 */
void readGC()
{
    FILE *stream;
    char buffer[1024], cell[256];
	double vmag;
    int page = 1;
    int entry = 0;
	int cumulus = 0;
	int nebulae = 0;
	int variables = 0;
	bool pendingRA = false;   /* la ultima estrella almacenada tiene segundos de RA ausentes/incompletos */
	bool pendingDecl = false; /* idem para Decl */
	int fixedRA = 0;
	int fixedDecl = 0;
	int unresolved = 0;
    GCstars = 0;
	
	// Lee Catálogo General Argentino
	stream = fopen("cat/gc.txt", "rt");
	if (stream == NULL) {
		perror("Cannot read gc.txt");
		exit(1);
	}

	vmag = 0.0;
	while (fgets(buffer, 1023, stream) != NULL) {
		entry++;
		if ((entry-53) % 70 == 0) {
			page++;
		}

		/* omite cualquier observacion que no sea la primera, salvo que la primera tenga
		   segundos ausentes/incompletos en RA o Decl: en tal caso se toman de la siguiente
		   observacion completa (las siguientes observaciones no traen magnitud ni precesiones,
		   por lo que solo se reemplaza la coordenada incompleta) */
		/* descarta las estrellas "1/2" (byte 6 = 1), distintas de la estrella con igual numeracion */
		if (buffer[6-1] == '1') continue;

		/* numero de observacion (byte 7) */
		readField(buffer, cell, 7, 1);
		if (atoi(cell) != 1) {
			if (!pendingRA && !pendingDecl) continue;
			readField(buffer, cell, 1, 5);
			struct GCstar_struct *st = &GCstar[GCstars - 1];
			if (atoi(cell) != st->gcRef) continue;
			int h, m, s;
			double val;
			/* solo se toman los segundos; horas/grados y minutos se mantienen de la primera observacion */
			if (pendingRA && readRA(buffer, &h, &m, &s, &val)) {
				st->RAs = s;
				st->RA1875 = 15.0 * (st->RAh + st->RAm/60.0 + (s/100.0)/3600.0);
				pendingRA = false;
				fixedRA++;
			}
			if (pendingDecl && readDecl(buffer, &h, &m, &s, &val)) {
				st->Decls = s;
				st->Decl1875 = -(st->Decld + st->Declm/60.0 + (s/10.0)/3600.0);
				pendingDecl = false;
				fixedDecl++;
			}
			sph2rec(st->RA1875, st->Decl1875, &st->x, &st->y, &st->z);
			continue;
		}
		if (pendingRA || pendingDecl) {
			printf("Warning: GC %d has incomplete seconds with no later complete observation\n", GCstar[GCstars - 1].gcRef);
			unresolved++;
		}

		/* lee numeracion */
		readField(buffer, cell, 1, 5);
		int gcRef = atoi(cell);

		/* ver si es cumulo, nebulosa o variable */
		char type = buffer[11-1];
		bool cum = false;
		bool neb = false;
		if (type == 'C') {
			cum = true;
			cumulus++;
		}
		if (type == 'N') {
			neb = true;
			nebulae++;
		}
		if (type == 'V') {
			vmag = 0.0;
			variables++;
		}
		else {
			/* lee magnitud (excepto si son espacios, en cuyo caso la magnitud y variabilidad es de la entrada anterior) */
			readField(buffer, cell, 8, 3);
			if (cell[0] != ' ') {
				if (cell[2] == ' ') cell[2] = '0';
				vmag = atof(cell)/10.0;
			}
		}

		/* lee epoca en que fue hecha la observacion */
		// readField(buffer, cell, 12, 4);
		// double epoch = (atof(cell)/100.0) + 1800.0;

		/* lee ascension recta y declinacion B1875.0 */
		int RAh, RAm, RAs, Decld, Declm, Decls;
		double RA, Decl;
		bool completeRA = readRA(buffer, &RAh, &RAm, &RAs, &RA);
		bool completeDecl = readDecl(buffer, &Decld, &Declm, &Decls, &Decl);

		/* lee precesiones y chequea, si es requerido */
        readFieldSanitized(buffer, cell, 24, 7);
        double preRA = atof(cell) / 1000.0;
        readFieldSanitized(buffer, cell, 46, 6);
        double preDecl = atof(cell) / 1000.0;

		/* calcula coordenadas rectangulares) */
		double x, y, z;
		sph2rec(RA, Decl, &x, &y, &z);

		if (GCstars == MAXGCSTAR) {
			printf("Max amount reached!\n");
			exit(1);
		}

		/* almacena la estrella */
		GCstar[GCstars].gcRef = gcRef;
		GCstar[GCstars].RAh = RAh;
		GCstar[GCstars].RAm = RAm;
		GCstar[GCstars].RAs = RAs;
		GCstar[GCstars].Decld = Decld;
		GCstar[GCstars].Declm = Declm;
		GCstar[GCstars].Decls = Decls;
		GCstar[GCstars].RA1875 = RA;
		GCstar[GCstars].Decl1875 = Decl;
		GCstar[GCstars].preRA = preRA;
		GCstar[GCstars].preDecl = preDecl;
		GCstar[GCstars].x = x;
		GCstar[GCstars].y = y;
		GCstar[GCstars].z = z;
		GCstar[GCstars].vmag = vmag;
		GCstar[GCstars].page = page;
		GCstar[GCstars].dpl = false;
		GCstar[GCstars].cum = cum;
		GCstar[GCstars].neb = neb;

		/* proxima estrella */
		GCstars++;
		pendingRA = !completeRA;
		pendingDecl = !completeDecl;
		//printf("Pos %d: id=%d RA=%.4f Decl=%.4f (%.2f) Vmag=%.1f\n", GCstars, gcRef, RA, Decl, epoch, vmag);
	}
	if (pendingRA || pendingDecl) {
		printf("Warning: GC %d has incomplete seconds with no later complete observation\n", GCstar[GCstars - 1].gcRef);
		unresolved++;
	}
	printf("Stars read from Catalogo General Argentino: %d\n", GCstars);
	printf("   Incomplete seconds fixed from later observations: RA %d, Decl %d (unresolved %d)\n",
		fixedRA,
		fixedDecl,
		unresolved);

	/* Ahora vamos a identificar las dobles */
	for (int i = 0; i < GCstars - 1; i++) {
		for (int j = i + 1; j < GCstars; j++) {
			double dist = 3600.0 * calcAngularDistance(GCstar[i].x, GCstar[i].y, GCstar[i].z, GCstar[j].x, GCstar[j].y, GCstar[j].z);
			if (dist < DPL_DISTANCE) {
				GCstar[i].dpl = true;
				GCstar[j].dpl = true;
			}
		}
	}
	int countDpl = 0;
	for (int i = 0; i < GCstars; i++) {
		if (GCstar[i].dpl) countDpl++;
	}
	printf("   Doubles: %d, Cumulus: %d, Nebulae: %d, Variables: %d\n",
		countDpl,
		cumulus,
		nebulae,
		variables);
}
