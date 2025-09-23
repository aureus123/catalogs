
/*
 * FIND_GSC - Header
 */

#include <stdbool.h>

#define MAX_DIST_GSC 20.0   // 20 segundos de arco

const char* getGSCId();
double getDist();
bool findGSCStar(double RA, double Decl, double epoch, double minDistanceOutput);
