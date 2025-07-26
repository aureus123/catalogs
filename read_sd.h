/*
 * READ_SD - Header
 */

#define MAXSDSTAR 150000

struct SDstar_struct {
    bool discard; /* true if should not be considered */
    int declRef; /* declinacion */
    int numRef; /* identificador con numero */
    float rah, ramin, raseg, decldeg, declmin; /* Ascension recta y declinacion, en partes */
    double RA1855, Decl1855, RA1875, Decl1875, vmag; /* coordenadas y magnitud visual */
    int cdIndex; /* Indice a CD */
    double x, y, z; /* coordenadas rectangulares en circulo unidad */
    double dist; /* distancia angular a su CD asociada (en arcsec) */
};

int getSDStars();
struct SDstar_struct *getSDStruct();
void writeRegisterSD(int sdIndex);
void readSD(bool onlyDecl22);
void findSDByCoordinates(double x, double y, double z, double decl, int *sdIndexOutput, double *minDistanceOutput); 