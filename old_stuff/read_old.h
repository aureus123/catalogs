
/*
 * READ_OLD - Header
 */

#define MAXGCSTAR 30000
#define MAX_MAGNITUDE 1.0
#define MAX_DISTANCE 360.0   // 6 minutos de arco
#define MAX_DIST_PPM 30.0    // 1/2 minuto de arco
#define MAX_DIST_CD 90.0    // 3/2 minuto de arco

struct GCstar_struct {
	bool discard; /* true si debe ser descartada (por doble) */
    int gcRef; /* identificador con numero */
    int obs; /* cantidad de observaciones */
	int RAh, RAm, RAs, Decld, Declm, Decls; /* detalle primera observacion */
	double x, y, z; /* primera obs. en coordenadas rectangulares */
    double RA1875[10], Decl1875[10], epoch[10]; /* coordenadas y epoca de la observacion */
    double vmag; /* magnitud visual */
    int page; /* pagina donde se encuentra */
	/* CD Cross-index */
    int cdIndex; /* indice (no identificador) a la estrella CD más cercana */
	double dist; /* distancia a la estrella CD */
	int cdIndexWithinMag; /* indice a la CD más cercana, pero dentro del rango de magnitud */
	double distWithinMag;
	/* Yarnall Cross-index */
	int yarnallRef; /* identificador a catalogo de la USNO */
	char yarnallCat[29]; /* referencia a otros catalogos */
	double distYarnall; /* distancia a USNO */
	double vmagYarnall; /* magnitud USNO */
	/* Weiss Cross-index */
	int weissRef, oeltzenRef; /* identificador a catalogo de Weiss y OA */
	double distWeiss; /* distancia a Weiss */
	double distOeltzen; /* distancia a OA */
	double vmagWeiss; /* magnitud Weiss */
	/* Stone Cross-index */
	int stoneRef, lacailleRef; /* identificador a catalogo de Stone y Lacaille */
	double distStone; /* distancia a Stone */
	double vmagStone; /* magnitud Stone */
	/* Gillis Cross-index */
	int giRef; /* identificador a catalogo de Gillis */
	char giCat[19]; /* referencia a otros catalogos */
	double distGi; /* distance a Gillis */
	double vmagGi; /* magnitud Gillis */
};

int getGCStars();
struct GCstar_struct *getGCStruct();
void writeRegisterGC(int index);
void readGC(bool mode, int fictRAh, int fictRAm, int fictRAs, int fictDecld, int fictDeclm, int fictDecls);
