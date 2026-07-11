
/*
 * CROSS_UTILS - Header
 * Codigo auxiliar compartido por cross_north y cross_south
 * (incluir despues de read_ppm.h, read_dm.h, read_gc.h, trig.h y misc.h)
 */

#define MAX_DIST_CROSS 120.0
#define MAX_DIST_CROSS_YARNALL 60.0
#define MAX_DIST_PPM_FAR 120.0

/* estadisticas acumuladas durante el cruzamiento de un catalogo */
struct CrossStats {
    int countDist = 0;
    double akkuDistError = 0.0;
    int countDelta = 0;
    double akkuDeltaError = 0.0;
    int countGSC = 0;
    int countCD = 0;
    int countCPD = 0;
    int errors = 0;
};

/* referencia a un catalogo almacenado en memoria (arreglos paralelos) */
struct StarList {
    int *count;
    int *ref;
    double *x, *y, *z;
};

/* lee PPM a la epoca dada y lo deja ordenado; devuelve la estructura */
struct PPMstar_struct *preparePPM(double epoch, bool discardSouth);

/* lee ascension recta en horas (col, 2), minutos (col+2, 2) y segundos
   (col+4, 4, en centesimas); devuelve grados */
double readRAField(char *buffer, int col, int *RAh, int *RAm, int *RAs);

/* lee declinacion (o NPD/SPD) en grados (colD, widthD), minutos (colM, 2) y
   segundos (colS, widthS, a dividir por divisor); devuelve grados positivos */
double readDeclField(char *buffer, int colD, int widthD, int colM, int colS, int widthS,
    double divisor, int *Decld, int *Declm, int *Decls);

/* lee magnitud entera (col, width) con posible '5' fraccional en colHalf */
float readMagIntHalf(char *buffer, int col, int width, int colHalf);

/* formatea una linea de registro "hh mm ss.ss sdd°mm'ss"s */
void formatCatLine(char *catLine, int RAh, int RAm, int RAs, char sign, int Decld, int Declm, int Decls);

/* abre/cierra la terna de archivos results/cross/cross_<tag>_{ppm,sao,hd}.csv */
void openCrossSet(const char *tag, FILE **ppm, FILE **sao, FILE **hd);
void closeCrossSet(FILE *ppm, FILE *sao, FILE *hd);

/* busca la PPM mas cercana; si esta a menos de maxDist acumula estadisticas,
   compara magnitudes (si magWarnName != NULL) y escribe el cruzamiento
   (si crossName != NULL); siempre devuelve el indice y distancia mas cercanos */
bool crossWithPPM(double x, double y, double z, double Decl, double vmag, double maxDist,
    const char *magWarnName, char *crossName, FILE *ppmStream, FILE *saoStream, FILE *hdStream,
    int *ppmIndexOut, double *nearestPPMDistance, struct CrossStats *stats);

/* si no hubo PPM cercana, prueba con GSC; si la halla y crossPPMStream != NULL,
   escribe el cruzamiento y cuenta en stats->countGSC */
bool tryGSC(bool ppmFound, double RA, double Decl, double epoch,
    FILE *crossPPMStream, char *catName, double vmag, struct CrossStats *stats);

/* busca la CPD mas cercana y genera el cruzamiento (si catName != NULL) */
bool crossWithCPD(double x, double y, double z, double Decl1875, double vmag, double maxDist,
    FILE *stream, char *catName, struct CrossStats *stats, int *cpdIndexOut, double *minDistOut);

/* busca la CD mas cercana, compara magnitudes (si magWarnName != NULL)
   y genera el cruzamiento (si catName != NULL) */
bool crossWithCD(double x, double y, double z, double Decl1875, double vmag,
    const char *magWarnName, FILE *stream, char *catName, struct CrossStats *stats,
    int *cdIndexOut, double *minDistOut);

/* almacena una estrella para futuras identificaciones; devuelve su indice */
int storeStar(int *count, int max, const char *name, int *ref, double *X, double *Y, double *Z,
    double *mag, int refValue, double x, double y, double z, double magValue);

/* revisa una referencia cruzada contra un catalogo almacenado: si la distancia
   supera MAX_DIST_CROSS advierte (con linea de registro si catLine != NULL y
   writeRegisterGC si gcIndexRegister >= 0), sino incrementa *check */
void checkCrossRef(const char *srcName, const char *catLine, const char *label,
    double x, double y, double z, int numRefCat, const struct StarList *list,
    bool breakAfterFirst, int gcIndexRegister, int *check, int *errors);

/* idem contra el catalogo WB, cuya identificacion es hora + numero */
void checkCrossRefWB(const char *srcName, const char *catLine, bool starWord,
    double x, double y, double z, int RARef, int numRefCat,
    const struct StarList *list, const int *hourRef, int gcIndexRegister,
    int *check, int *errors);

/* idem contra el catalogo Durchmusterung (BD/CD) por su identificador */
void checkBDRef(const char *srcName, const char *catLine, bool printNotFound,
    bool signRef, int declRef, int numRefCat,
    double x, double y, double z, int *check, int *errors);

/* idem contra el catalogo GC ya leido con readGC() */
void checkCrossRefGCcat(const char *srcName, const char *catLine, int numRefCat,
    double x, double y, double z, int *check, int *errors);

/* revisa una referencia a Yarnall (USNO 2a ed.): como su numeracion difiere de
   la de Yarnall-Frisby (USNO 3a ed.), busca la estrella USNO mas cercana */
void checkYarnallRef(const char *srcName, const char *catLine, int numRefCat,
    double x, double y, double z, bool requirePositiveIndex, int gcIndexRegister,
    const struct StarList *list, int *check, int *errors);

/* advierte estrella sola (sin PPM cercana ni GSC), estilo cross_north */
void warnIfAloneNorth(bool ppmFound, double minDistance, double RA, double Decl, double epoch,
    const char *warnDesc, int *errors);

/* advierte estrella sola (no PPM / CD / CPD / GSC) y registra causas */
void warnAlone(int *errors, const char *warnDesc, const char *registerDesc, const char *catLine,
    char *catName, int RAs, double decl, int Decls, int ppmRef, double nearestPPMDistance);

/* advierte estrella sola (no PPM or GSC) y registra causas */
void warnAlonePPMGSC(int *errors, const char *warnDesc, char *catName,
    int RAs, double decl, int Decls, int ppmRef, double nearestPPMDistance);

/* almacena una estrella no identificada (pero hallada en GSC) */
void writeUnidentified(FILE *stream, const char *catName, double x, double y, double z);

/* imprime los RSME acumulados */
void printRSMEDist(const struct CrossStats *stats);
void printRSMEMag(const struct CrossStats *stats);

/* revisa la precesion anual reportada contra la calculada (ctes. de Struve) */
void checkPrecessionRA(double preRA, const char *srcName, const char *catLine,
    double RA, double Decl, double base, double factor, double tol, int *errors);
void checkPrecessionDecl(double preDecl, const char *srcName, const char *catLine,
    double RA, double coef, double tol, int *errors);

/* extrae numero GC y referencia (3a columna) de una fila escaneada de GC;
   false si no hay referencia confiable */
bool parseGCScanLine(char *buffer, int *gcRef, char *ref);

/* parsea una fila de cat/UA_standards.csv:
   number,name,ra_h,ra_m,ra_s,dec_d,dec_m,mag,BD_num,BD_mag,
   Lal_num,Lal_mag,WB_num,WB_mag,ARGEL,Albany,HEIS
   false si la linea es incompleta o cabecera */
bool parseSOMLine(char *buffer, char field[17][32], char *catName, char *catLine,
    double *RA, double *Decl);

/* datos comunes de una fila de cat/ua.txt */
struct UAstar_struct {
    double RA, Decl;
    int RAh, RAm, RAs, Decld;
    double Declm;
    bool existsRef;
    int gouldRef;
    char cstRef[5];
    char catgName[20];
    char catLine[64];
};

/* parsea una fila de cat/ua.txt (coordenadas, numeracion de Gould y
   constelacion); false si la estrella no tiene coordenadas */
bool parseUALine(char *buffer, char *serpens, struct UAstar_struct *ua);

/* separa las referencias a catalogos (columnas 82-98) en subceldas */
int splitUARefs(char *buffer, char subcell[3][18]);
