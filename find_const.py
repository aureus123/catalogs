#!/usr/bin/env python3
"""
Script para encontrar la constelación dadas coordenadas ecuatoriales celestes (J2000.0)
Uso: python3 find_const.py RA_h RA_m RA_s Dec_d Dec_m Dec_s
Ejemplo: python3 find_const.py 11 22 05.29 -24 46 39.8
"""

import sys
from astropy.coordinates import SkyCoord
from astropy import units as u


def parse_coordinates(args):
    """
    Parsea los argumentos de línea de comandos para extraer RA y Dec.
    
    Args:
        args: Lista con [RA_h, RA_m, RA_s, Dec_d, Dec_m, Dec_s]
    
    Returns:
        tuple: (ra_string, dec_string) en formato para SkyCoord
    """
    if len(args) != 6:
        print("Error: Se requieren 6 argumentos (RA: h m s, Dec: d m s)", file=sys.stderr)
        print("Uso: python3 find_const.py RA_h RA_m RA_s Dec_d Dec_m Dec_s", file=sys.stderr)
        print("Ejemplo: python3 find_const.py 11 22 05.29 -24 46 39.8", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Construir strings para RA y Dec en formato HMS y DMS
        ra_string = f"{args[0]}h{args[1]}m{args[2]}s"
        dec_string = f"{args[3]}d{args[4]}m{args[5]}s"
        
        return ra_string, dec_string
        
    except ValueError as e:
        print(f"Error: Los argumentos deben ser números válidos: {e}", file=sys.stderr)
        sys.exit(1)


def find_constellation(ra_string, dec_string):
    """
    Encuentra la constelación para las coordenadas dadas.
    
    Args:
        ra_string: Ascensión recta en formato "HhMmSs"
        dec_string: Declinación en formato "DdMmSs"
    
    Returns:
        str: Abreviación de la constelación (3 letras)
    """
    # Crear objeto SkyCoord con las coordenadas en J2000.0 (ICRS)
    coord = SkyCoord(ra_string, dec_string, frame='icrs')
    
    # Obtener la constelación
    constellation = coord.get_constellation(short_name=True)
    
    return constellation


def main():
    # Obtener argumentos de línea de comandos (excluyendo el nombre del script)
    args = sys.argv[1:]
    
    # Parsear coordenadas
    ra_string, dec_string = parse_coordinates(args)
    
    # Encontrar constelación
    constellation = find_constellation(ra_string, dec_string)
    
    # Mostrar resultado
    print(constellation)


if __name__ == "__main__":
    main()
