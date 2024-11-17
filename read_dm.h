
/*
 * READ_DM - Header (CD y BD)
 */

#define MAXDMSTAR 613959

/* Nota: una estrella se distingue con los 4 identificadores:
     signRef, declRef, numRef y supplRef */
struct DMstar_struct {
    bool signRef; /* true if sign of declRef is negative */
    int declRef, numRef; /* identificador con declinacion y numero */
    char supplRef; /* identificador suplementario (a,b,c) */
    float rah, ramin, raseg, decldeg, declmin; /* Ascension recta y declinacion, en partes */
    double RA1855, Decl1855; /* coordenadas (BD) */
    double RA1875, Decl1875; /* coordenadas (CD) */
    double vmag; /* magnitud visual */
    int catIndex; /* Indice al indice del otro catálogo en caso de existir, o -1 sino */
    double x, y, z; /* coordenadas rectangulares en circulo unidad */
    int next; /* lista a la siguiente estrella de misma decl y num, o -1 = última */
};

int getDMStars();
int getDMindex(bool signRef, int declRef, int numRef);
bool isCD();
struct DMstar_struct *getDMStruct();
void writeRegister(int dmIndex, bool neighbors);
void readDM(const char *filename);
void findByCoordinates(double x, double y, double z, int *index, double *minDistance);
