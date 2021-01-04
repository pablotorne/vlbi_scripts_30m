#!/usr/bin/env python
"""
Script to create the Tsys* tables formatted with the standard format
requested by EHT Data Managers in Dec. 2020.

Takes as input:
- a text file containing the calibration information given by
  "read_lastcal_info.py", normally a Field System log
- a text file with the information of the VLBI scans of ONE SCHEDULY,
  typically a .prn file generated from the .vex at the session start.

(This code is based on create_ANTAB.py)

v2: Creates first a file with all the calibration scan metadata and
    the  uses this files to produce the Tsys* table for VLBI scans.
    The strategy in this version is:
    1) If a calibration(s) scan(s) for the science VLBI scan exists
      within 10 minutes of the VLBI scan start time, use the cal. data
      from the CLOSEST cal. scan.
    2) If 1) does not  fullfil,  take the two closest cal. scans on the
       VLBI source before and after the VLBI scan start time and use a
       linear interpolation to derive the values of Tsys*
    3) If 1) or 2) are not possible, use the closest-in-time cal.
       scan on any source and calculate the Tsys* in the elevation
       of the VLBI scan
  
P. Torne, IRAM 03.01.2020
"""

import sys, datetime, os
import numpy as np
import subprocess
import argparse

# Arguments control
parser = argparse.ArgumentParser()

parser.add_argument("-v", "--verbose", action="count", help="outputs detailed execution informaton")
parser.add_argument("calinfo_file", help="input file where to read the calibration information (e.g. a Field System log). This code is written for Pv, and assumes that the text file contains rows with the 'rx:' strings to fetch the calibration data as outputted by the read_lastcal_info.py routine.")
parser.add_argument("prn_file", help="input file where to read the VLBI scans number and timestamps (use a .prn file). The output file name of the Tsys table will be generated from the baseline of this file (which should be the same as the schedule name.)")

args = parser.parse_args()

#sched = args.sched.split(".vex")[0]  # Get the basename of the schedule name
basenm = os.path.basename( os.path.splitext(args.prn_file)[0] )

# Fetch info from input file
try:

    calinfo = subprocess.check_output('grep %s -e "rx:"'%args.calinfo_file, shell=True).splitlines()
    prninfo  = subprocess.check_output('grep -E "no[0-9]{1,}" %s'%args.prn_file, shell=True).splitlines()

except Exception as e:

    print e
    print "\nHalting program."
    sys.exit(1)

# *************************************************
# Output CALIBRATION info in EHT Tsys Table format
# *************************************************

print "\nOpening file %s_PV.CAL.txt to write Tsys Metadata table ..."%basenm
outputfn = open("%s_PV.CAL.txt"%basenm, "w")

# Write the header of the Tsys table file:
outputfn.write("################################################################\n")
outputfn.write("#\n")
outputfn.write("# Tsys/Tsys* table for each CALIBRATION scan from PV: %s\n"%basenm)
outputfn.write("#\n")
outputfn.write("################################################################\n")

all_timestamp = [] # always save last time stamp to avoid duplicates in the output file

print "Reading information from %d calibration scans from %s ..."%(len(calinfo)/4, args.calinfo_file)
for ii in range(0, len(calinfo), 4):  # each cal scan info comes in 4 rows for 4 IFs:
    # Read and format data:
    rx        = [] # Rx Name
    tsys_star = [] # Tsys with airmass correction applied
    tau_z     = [] # Tau at observing frequency at Zenith
    pwv_mm    = [] # Calculation of pwv. For 1.3mm OK, for other bands not so good
    trx       = [] # Trec in K
    timestamp = [] # UTC time
    rxFreq    = [] # Central frequency of the subband
    azimuth   = [] # in deg
    elevation = [] # in deg
    source    = [] 
    gainImage = [] # For 2SB, it is sideband coupling, (Usually for EMIR = -13dB = 5%)
    Tamb      = [] # Physical temperature of the ambient load used in calibration
    Tatm_s    = [] # Physical temperature of the atmopshere in the signal band used in calibration
    DOY       = []

    for jj in range(ii, ii+4): # Loop over four lines starting from ii

        data = calinfo[jj].split(" ; ")

        if args.verbose > 1: print "Reading line: %s"%data

        # Extract DOY and timestamp
        yy = int( data[5].split()[1].split("-")[0] )
        mm = int( data[5].split()[1].split("-")[1] )
        dd = int( data[5].split()[1].split("-")[2] )
 
        dayofyear = datetime.datetime(yy, mm, dd).timetuple().tm_yday
        
        if args.verbose >= 1:
            print "Date from cal data: %s-%s-%s"%(yy, mm, dd)
            print "DOY = datetime.datetime(yy, mm, dd).timetuple().tm_yday = %d"%dayofyear 

        ts = data[5].split("time:")[-1].lstrip()

        DOY.append(dayofyear)
        timestamp.append(ts)

        if jj == ii: # Ony append once to the all_timestamps vector
            all_timestamp.append(ts)

        # Extract Tsys values
        tempsys = float( data[1].split()[1] )
        # Check value
        if tempsys > 9999.9: 
            tsys_star.append(9999.9)
        elif (tempsys > 1 and tempsys <= 9999.9):
            tsys_star.append(tempsys)
        elif tempsys <= 0:
            tsys_star.append("NA")

        # Extract tau
        tau_z.append( float(data[2].split()[1]) )

        try:

            # Extract elevation
            azimuth.append( float(data[7].split()[1]) )

            # Extract elevation
            elevation.append( float(data[8].split()[1]) )

            # Extract the source
            source.append( data[9].split()[1] )

        except Exception as e:
            print "\n* Error: You seem to be passing a log from before Apr2018. Sorry, those are not in a compatible format.\n"
            sys.exit(1)

        # Extract GainImage, Tamb and Tatm_s (and other fields: Tcold, Pcold, Phot, Psky, effForward, effBeam, backend)
        gainImage.append( float(data[10].split()[1]) )
        #Tcold _TBD
        Tamb.append( float(data[12].split()[1]) )
        Tatm_s.append( float(data[13].split()[1]) )
        #Pcold _TBD
        #Phot _TBD
        #Psky _TBD
        #effForward
        #effBeam
        #Backend 

    if args.verbose >= 1:
        #print "DOY = %s"%DOY
        #print "timestamp = %s"%timestamp
        #print "tsys = %s"%tsys
        #print "source = %s"%source
        print "Extracting relevant information for a cal scan on %s"%source[-1]

    # Check:
    if DOY[1:] == DOY[:-1] and timestamp[1:] == timestamp[:-1] and source[1:] == source [:-1]: # these are four lines corresponding to the same scan
        if args.verbose >= 1: print "Writing to calibration info on %s to Tsys* table ..."%source[-1]

        if ii == 0: # First entry, no duplicity check

            # Write one line with the requested metadata:
            # Note: at 64Gbps, we only record LowerInner and UpperInner
            # In EHT/Correlator nomenclature: Tsys* of LI = Tsys* of Band1 *and* Band2 (we do not measure Tsys separately for them)
            #                                 Tsys* of UI = Tsys* of Band3 *and* Band4 (and H=RCP and V=LCP)
            outputfn.write("%s\tCALXML\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n"%(\
                           timestamp[-1], source[-1], azimuth[-1], elevation[-1], \
                           tsys_star[0], tsys_star[2], tsys_star[0], tsys_star[2], \
                           tsys_star[1], tsys_star[3], tsys_star[1], tsys_star[3], \
                           np.mean((tau_z[0], tau_z[2], tau_z[1], tau_z[3])), \
                           np.mean((Tamb[0], Tamb[2], Tamb[1], Tamb[3])), \
                           np.mean((Tatm_s[0], Tatm_s[2], Tatm_s[1], Tatm_s[3])), \
                           ))

        elif ii > 0: # Ignore duplicity check for first entry

            if timestamp[-1] == all_timestamp[-2]:  # Do not write repeated cal info (occurs if no new cal is done when calling getcals in FS)
                print "timestamp = %s | previous = %s"%(timestamp[-1], all_timestamp[-2])
                print "Detected duplicated timestamp. Ignoring line."
                continue

            else:
                outputfn.write("%s\tCALXML\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n"%(\
                           timestamp[-1], source[-1], azimuth[-1], elevation[-1], \
                           tsys_star[0], tsys_star[2], tsys_star[0], tsys_star[2], \
                           tsys_star[1], tsys_star[3], tsys_star[1], tsys_star[3], \
                           np.mean((tau_z[0], tau_z[2], tau_z[1], tau_z[3])), \
                           np.mean((Tamb[0], Tamb[2], Tamb[1], Tamb[3])), \
                           np.mean((Tatm_s[0], Tatm_s[2], Tatm_s[1], Tatm_s[3])), \
                           ))


print "Writing to %s_PV.CAL.txt ..."%basenm

    #if ii >= 4: break

# Write an 'end of file' marker
outputfn.write("/\n")

outputfn.close()

# Now we have a text file containing all the calibration data formatted properly,
# From this and the .prn, we need to create the Tsys* table assigning to each VLBI
# scan the best / corresponding Tsys*, etc.  information.


print "DONE."
