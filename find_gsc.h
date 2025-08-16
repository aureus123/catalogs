/*
 * FIND_GSC - Interface to Guide Star Catalogue (GSC) tool
 * Made in 2025 by Daniel E. Severin
 * 
 * This module provides functionality to search the Guide Star Catalogue
 * for stars near given coordinates using the external GSC tool.
 */

#include <stdbool.h>

/*
 * findGSCStar - Search for stars in the Guide Star Catalogue near given coordinates
 * 
 * Parameters:
 *   RA              - Right ascension in decimal degrees
 *   Decl            - Declination in decimal degrees  
 *   epoch           - Epoch of coordinates in Besselian years (e.g., 1875.0, 1850.0)
 *   minDistanceOutput - Search threshold in arcseconds
 * 
 * Returns:
 *   true if GSC stars found within threshold, false otherwise
 * 
 * Notes:
 *   - Coordinates are automatically precessed to J2000.0 before GSC search
 *   - The GSC tool expects Julian epoch coordinates and arcminute radius
 *   - Uses GSCDAT environment variable to locate GSC data
 *   - Calls external gsc tool: gsc/bin/gsc -c RA Decl -r radius -p 2
 */
bool findGSCStar(double RA, double Decl, double epoch, double minDistanceOutput);
