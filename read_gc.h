
/*
 * READ_GC - Header
 */

#define MAXGCSTAR 33000
#define MAX_MAGNITUDE 2.0
#define MAX_DIST_CD 90.0     // 1.5 minutos de arco
#define MAX_DIST_CPD 20.0    // 20 segundos de arco
#define MAX_DIST_PPM 15.0    // 15 segundos de arco
#define DPL_DISTANCE 30.0    // 30 segundos de arco para considerar que una estrella es doble

struct GCstar_struct {
    int gcRef; /* identificador con numero */
	int RAh, RAm, RAs, Decld, Declm, Decls; /* detalle primera observacion */
	double RA1875, Decl1875; /* coordenadas ecuatoriales (1875) */
	double x, y, z; /* primera obs. en coordenadas rectangulares (1875) */
    double vmag; /* magnitud visual */
    int page; /* pagina donde se encuentra */
	bool dpl; /* true si hay otra estrella GC muy cerca */
	bool cum; /* true si es cumulo */
	bool neb; /* true si es nebulosa */
};

int getGCStars();
struct GCstar_struct *getGCStruct();
void writeRegisterGC(int index);
void readGC();
void logCauses(bool cumulus, bool nebula, double vmag, int RAs, double Decl, int Decls);
