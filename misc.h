
/*
 * MISC - Header
 */

void bye(const char *string);
void readField(char *buffer, char *cell, int initial, int bytes);
void readFieldSanitized(char *buffer, char *cell, int initial, int bytes);
FILE *openCrossFile(const char *name);
void writeCrossEntry(FILE *stream, char *index1, char *index2, double dist);
FILE *openPositionFile(const char *name);
void writePositionEntry(FILE *stream, int index, int decl, int num, const char *ref, double dist);
FILE *openMagnitudeFile(const char *name);
void writeMagnitudeEntry(FILE *stream, int index, int decl, int num, const char *ref, double delta);
