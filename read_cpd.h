
/*
 * READ_CPD - Header
 */

#define MAXCPDSTAR 454812

struct CPDstar_struct {
    bool discard; /* true if should not be considered */
    int declRef, numRef; /* identificador con declinacion y numero */
    double RA1875, Decl1875, pmag; /* coordenadas y magnitud fotografica */
    int dmIndex; /* Indice a DM (CD) */
    double x, y, z; /* coordenadas rectangulares en circulo unidad */
    double dist; /* distancia angular a su CD asociada (en arcsec) */
    int declRef2, numRef2; /* identificador a CD de cat. 4019 */
};

int getCPDStars();
struct CPDstar_struct *getCPDStruct();
void readCPD(bool cross, bool catalog);
void findCPDByCoordinates(double x, double y, double z, double decl, int *cpdIndexOutput, double *minDistanceOutput);
