#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Script to automatically check for duplicated sources in the paKo catalog
of a VLBI session to aid when reducing the data for K/Jy factors.

Written in Python 2.7.9 for compatibility in mrt-lx3.

Set target filename and angular distance for matches in main() method

P. Torne + ChatGPT, IRAM v2025.07.02
'''

import math
import argparse

def sexagesimal_to_degrees_ra(s):
    try:
        parts = s.strip().split(':')
        h = float(parts[0])
        m = float(parts[1])
        s_ = float(parts[2])
        return (h + m/60.0 + s_/3600.0) * 15.0
    except:
        return None

def sexagesimal_to_degrees_dec(s):
    try:
        parts = s.strip().split(':')
        sign = -1 if parts[0].startswith('-') else 1
        d = abs(float(parts[0]))
        m = float(parts[1])
        s_ = float(parts[2])
        return sign * (d + m/60.0 + s_/3600.0)
    except:
        return None

def parse_sources(file_path):
    sources = []

    with open(file_path, 'r') as f:
        for line in f:
            if not line.strip() or line.strip().startswith('!'):
                continue

            parts = line.strip().split()

            # Require at least 5 fields: name, equinox name, equinox value, RA, Dec
            if len(parts) < 5:
                continue

            ra_deg = sexagesimal_to_degrees_ra(parts[3])
            dec_deg = sexagesimal_to_degrees_dec(parts[4])

            if ra_deg is None or dec_deg is None:
                continue

            source = {
                'name': parts[0],
                'ra': ra_deg,
                'dec': dec_deg,
                'line': line.strip()
            }

            sources.append(source)

    return sources

def angular_distance(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    '''
    Computes the angular separation between two points on the celestial sphere.

    Inputs:
        ra1_deg, dec1_deg: Right Ascension and Declination of point 1 (in degrees)
        ra2_deg, dec2_deg: Right Ascension and Declination of point 2 (in degrees)

    Returns:
        Angular distance in degrees.

    Method:
        Uses the Vincenty formula (spherical version), which is:
        
        ΔRA = RA2 - RA1
        x = cos(Dec2) * sin(ΔRA)
        y = cos(Dec1) * sin(Dec2) - sin(Dec1) * cos(Dec2) * cos(ΔRA)
        z = sin(Dec1) * sin(Dec2) + cos(Dec1) * cos(Dec2) * cos(ΔRA)

        angle = atan2(√(x² + y²), z)
        
        This gives the angular distance in radians between the two positions.
        Then it's converted to degrees.

    This method is accurate even for very small separations (< 1 arcsecond).
    '''

    # Convert input coordinates from degrees to radians
    ra1 = math.radians(ra1_deg)
    dec1 = math.radians(dec1_deg)
    ra2 = math.radians(ra2_deg)
    dec2 = math.radians(dec2_deg)

    # Compute the difference in right ascension
    delta_ra = ra2 - ra1

    # Compute the intermediate x, y, and z terms for the Vincenty (spherical) formula
    x = math.cos(dec2) * math.sin(delta_ra)
    y = math.cos(dec1) * math.sin(dec2) - math.sin(dec1) * math.cos(dec2) * math.cos(delta_ra)
    z = math.sin(dec1) * math.sin(dec2) + math.cos(dec1) * math.cos(dec2) * math.cos(delta_ra)

    # atan2 returns the angular separation in radians
    angle_rad = math.atan2(math.sqrt(x**2 + y**2), z)

    # Convert to degrees before returning
    return math.degrees(angle_rad)


def find_duplicates(sources, max_sep_arcsec, ignore_same_name):
    duplicates = []
    max_sep_deg = max_sep_arcsec / 3600.0

    for i in range(len(sources)):
        for j in range(i + 1, len(sources)):
            src1 = sources[i]
            src2 = sources[j]

            if ignore_same_name and src1['name'] == src2['name']:
                continue

            dist_deg = angular_distance(src1['ra'], src1['dec'], src2['ra'], src2['dec'])

            if dist_deg < max_sep_deg:
                duplicates.append((src1, src2, dist_deg))

    return duplicates

def main():
    parser = argparse.ArgumentParser(description="Find duplicate source pairs based on angular distance.")
    parser.add_argument('--file', type=str, required=True, help="Input text file with sources.")
    parser.add_argument('--threshold', type=float, default=3.0, help="Angular distance threshold in arcseconds.")
    parser.add_argument('--ignore_same_name', action='store_true', help="Ignore matches between sources with identical names.")
    args = parser.parse_args()

    sources = parse_sources(args.file)
    duplicates = find_duplicates(sources, args.threshold, args.ignore_same_name)

    if not duplicates:
        print "No duplicates found within %.3f arcsec." % args.threshold
    else:
        print "Duplicate source pairs (closer than %.3f arcsec):\n" % args.threshold
        for src1, src2, dist in duplicates:
            print "Duplicate:\n"
            print "%s" % src1['line']
            print "%s" % src2['line']
            print "Angular separation: %.6f arcsec\n" % (dist * 3600.0)

if __name__ == '__main__':
    main()

