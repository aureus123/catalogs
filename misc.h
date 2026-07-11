
/*
 * MISC - Header
 */

void bye(const char *string);
void readField(char *buffer, char *cell, int initial, int bytes);
void readFieldSanitized(char *buffer, char *cell, int initial, int bytes);
FILE *openCrossFile(const char *name);
void writeCrossEntry(FILE *stream, char *index1, char *index2, double mag, double dist);
FILE *openUnidentifiedFile(const char *name);
FILE *openCatalogFile(const char *name);
void writeCatalogFile(FILE *stream, const char *name, double x, double y, double z, double mag);
void logCauses(char *name, bool durchCoverage,
    bool cumulus, bool nebula,
    int RAs, double Decl, int Decls,
    int ppmRef, double nearestPPMDistance);
void copyWithoutSpaces(char *dest, char *src);
