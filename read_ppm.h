
/*
 * READ_PPM - Header
 */

#define MAXPPMSTAR 468861

struct PPMstar_struct {
    int ppmRef; /* identificador PPM */
    bool discard; /* true if should not be considered (e.g. other PPM star nearer to the same CD record) */
    double polarDist; /* distancia polar desde el sur (para sort) */
    double vmag; /* Magnitud Johnson V. Si vmag = 0.0 la magnitud original era fotografica y no la consideramos */
    char problem; /* = 1 solo si la estrella esta marcada como doble o "problematica" */
    int dmIndex; /* Indice a BD o CD */
    char dmString[14]; /* String con la identificación */
    double x, y, z; /* coordenadas rectangulares en circulo unidad */
    double dist; /* distancia angular a su CD asociada (en arcsec) */
};

int getPPMStars();
struct PPMstar_struct *getPPMStruct();
bool revise(int ppmIndex);
void readPPM(bool useDurch, bool allSky, bool discard_north, bool discard_south, double targetYear);
void sortPPM();
void findPPMByCoordinates(double x, double y, double z, double decl, int *ppmIndex, double *minDistance);
