
/*
 * READ_PPM - Header
 */

#define MAXPPMSTAR 468586

struct PPMstar_struct {
    int ppmRef; /* identificador PPM */
    bool discard; /* true if should not be considered (e.g. other PPM star nearer to the same CD record) */
    double RA1855, Decl1855; /* posiciones para comparar contra BD */
    double RA1875, Decl1875; /* posiciones para comparar contra CD */
    double vmag; /* Magnitud Johnson V. Si vmag = 0.0 la magnitud original era fotografica y no la consideramos */
    char problem; /* = 1 solo si la estrella esta marcada como doble o "problematica" */
    int dmIndex; /* Indice a BD o CD */
    double x, y, z; /* coordenadas rectangulares en circulo unidad */
    double dist; /* distancia angular a su CD asociada (en arcsec) */
};

int getPPMStars();
struct PPMstar_struct *getPPMStruct();
bool revise(int ppmIndex);
void readPPM(bool allSky);
