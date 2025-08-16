/*
 * FIND_GSC - Interface to Guide Star Catalogue (GSC) tool
 * Made in 2025 by Daniel E. Severin (partially created by Cursor AI)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "find_gsc.h"
#include "trig.h"
#include "misc.h"

// #define FAKE_GSC    // uncomment this line to avoid using "gsc" tool

/*
 * findGSCStar - Search for stars in the Guide Star Catalogue near given coordinates
 */
bool findGSCStar(double RA, double Decl, double epoch, double minDistanceOutput) {
#ifdef FAKE_GSC
    return true;
#else    
    char command[512];
    char result_line[256];
    FILE *pipe;
    
    // Convert coordinates from Besselian epoch to J2000.0
    double RAtarget = RA;
    double Decltarget = Decl;
    transform(epoch, 2000.0, &RAtarget, &Decltarget);
    
    // Convert threshold from arcseconds to arcminutes
    double radiusArcmin = minDistanceOutput / 60.0;
    
    // Ensure coordinates are in valid ranges
    if (RAtarget < 0.0) RAtarget += 360.0;
    if (RAtarget >= 360.0) RAtarget -= 360.0;
    
    // Build GSC command with proper environment variable and path
    const char* base_path = getenv("HOME");
    if (base_path == nullptr) bye("Error: HOME environment variable is not set.\n");
    
    snprintf(command, sizeof(command), 
                "GSCDAT=%s/catalogs/gsc %s/catalogs/gsc/bin/gsc -c \"%.8f %+.8f\" -r %.2f -p 2 2>/dev/null",
                base_path, base_path, RAtarget, Decltarget, radiusArcmin);
    
    // Execute GSC tool
    pipe = popen(command, "r");
    if (pipe == nullptr) {
        printf("Warning: Could not execute GSC tool\n");
        return false;
    }
    
    // Check if we got any results
    bool found = false;
    
    while (fgets(result_line, sizeof(result_line), pipe) != nullptr) {
        // Skip empty lines and comments
        if (strlen(result_line) < 10 || result_line[0] == '#' || result_line[0] == '!') {
            continue;
        }
        
        // Check if this looks like a GSC star result (starts with digits)
        if (isdigit(result_line[0])) {
            found = true;
            
            // Extract GSC ID (first field)
            char gscId[16];
            float gscRA, gscDecl;
            sscanf(result_line, "%15s %f %f", gscId, &gscRA, &gscDecl);
            double gx, gy, gz;
            sph2rec((double) gscRA, (double)gscDecl, &gx, &gy, &gz);
            double tx, ty, tz;
            sph2rec(RAtarget, Decltarget, &tx, &ty, &tz);    
            double dist = 3600.0 *calcAngularDistance(gx, gy, gz, tx, ty, tz);
            
            printf("  Note: found GSC star %s near the target star at %.2f arcsec.\n", gscId, dist);
        }

        if (found) { break; }
    }
    pclose(pipe);
    return found;
#endif
}
