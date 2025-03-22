
/*
 * READ_GC - Lee catálogo general argentino
 * Made in 2025 by Daniel E. Severin
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
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

		/* omite cualquier observacion que no sea la primera */
		readField(buffer, cell, 6, 2);
		if (atoi(cell) != 1) continue;

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

		/* lee ascension recta B1875.0 */
		readFieldSanitized(buffer, cell, 16, 2);
		int RAh = atoi(cell);
		double RA = (double) RAh;
		readFieldSanitized(buffer, cell, 18, 2);
		int RAm = atoi(cell);
		RA += ((double) RAm)/60.0;
		readFieldSanitized(buffer, cell, 20, 4);
		int RAs = atoi(cell);
		RA += (((double) RAs)/100.0)/3600.0;
		RA *= 15.0; /* conversion horas a grados */

		/* lee declinacion B1875.0 */
		readFieldSanitized(buffer, cell, 39, 2);
		int Decld = atoi(cell);
		double Decl = (double) Decld;
		readFieldSanitized(buffer, cell, 41, 2);
		int Declm = atoi(cell);
		Decl += ((double) Declm)/60.0;
		readFieldSanitized(buffer, cell, 43, 3);
		int Decls = atoi(cell);
		Decl += (((double) Decls)/10.0)/3600.0;
		Decl = -Decl; /* incorpora signo negativo (en nuestro caso, siempre) */

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
		//printf("Pos %d: id=%d RA=%.4f Decl=%.4f (%.2f) Vmag=%.1f\n", GCstars, gcRef, RA, Decl, epoch, vmag);
	}
	printf("Stars read from Catalogo General Argentino: %d\n", GCstars);

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

/*
 * logCauses - escribe posibles causas de falta de identificacion
 */
void logCauses(bool cumulus, bool nebula, double vmag, int RAs, double Decl, int Decls, int ppmIndex, double nearestPPMDistance) {
    if (cumulus) {
        printf("  Possible cause: cumulus.\n");
    }
    if (nebula) {
        printf("  Possible cause: nebula.\n");
    }
    if (RAs == 0) {
        printf("  Possible cause: lack of RA (s).\n");
    }
    if (Decls == 0) {
        printf("  Possible cause: lack of Decl (s).\n");
    }
    if (vmag > 9.9) {
        printf("  Possible cause: dim star.\n");
    }
    if (vmag < 0.1) {
        printf("  Possible cause: no magnitude (cumulus?).\n");
    }
	if (ppmIndex != -1) {
		printf("  Note: Nearest PPM %d at %.1f arcsec.\n", ppmIndex, nearestPPMDistance);
	}    
    if (Decl > -18.0) {
        printf("  Note: no CD/CPD coverage for stars below 18°.\n");
    }
    if (Decl < -61.0) {
        printf("  Note: poor CD coverage for stars above 61°.\n");
    }
}
