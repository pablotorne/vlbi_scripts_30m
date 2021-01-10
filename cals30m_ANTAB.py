#!/usr/bin/env python
"""
Reads the calibration information from an input file (usually a FS log),
and outputs it in ANTAB format for VLBI correlator.

v0: the cal info must be in the format output by
    read_lastcal_info.py (used in VLBI Field System), i.e.: 
    (now in the version of 2021)

rx: E2HLI ; tsys*: 386.02 ; tau_z: 0.36 ; pwv mm: 6.0 ; trx: 79.0 ; time: 2018-04-20 21:28:29 ; rxFreq: 214.8497 ; azim: 144.3 ; elev: 40.5 ; source: 3C279 ; gainImage: 0.050 ; Tcold: 33.259 ; Tamb: 293.725
 ; Tatm: 258.021 ; Pcold: 149666.71875 ; Phot: 496961.56250 ; Psky: 265534.84375 ; effForward: 0.94 ; effBeam: 0.94 ; backend: NBC
rx: E2HUI ; tsys*: 392.51 ; tau_z: 0.36 ; pwv mm: 5.8 ; trx: 79.9 ; time: 2018-04-20 21:28:29 ; rxFreq: 227.3500 ; azim: 144.3 ; elev: 40.5 ; source: 3C279 ; gainImage: 0.050 ; Tcold: 33.259 ; Tamb: 293.725
 ; Tatm: 258.245 ; Pcold: 147705.71875 ; Phot: 487830.15625 ; Psky: 262296.43750 ; effForward: 0.94 ; effBeam: 0.94 ; backend: NBC
rx: E2VLI ; tsys*: 407.69 ; tau_z: 0.35 ; pwv mm: 5.8 ; trx: 96.2 ; time: 2018-04-20 21:28:29 ; rxFreq: 214.8497 ; azim: 144.3 ; elev: 40.5 ; source: 3C279 ; gainImage: 0.050 ; Tcold: 33.259 ; Tamb: 293.725
 ; Tatm: 258.247 ; Pcold: 198666.14062 ; Phot: 598307.12500 ; Psky: 328366.43750 ; effForward: 0.94 ; effBeam: 0.94 ; backend: NBC
rx: E2VUI ; tsys*: 422.96 ; tau_z: 0.36 ; pwv mm: 5.8 ; trx: 96.7 ; time: 2018-04-20 21:28:29 ; rxFreq: 227.3500 ; azim: 144.3 ; elev: 40.5 ; source: 3C279 ; gainImage: 0.050 ; Tcold: 33.259 ; Tamb: 293.725
 ; Tatm: 258.293 ; Pcold: 196851.14062 ; Phot: 591283.31250 ; Psky: 328991.00000 ; effForward: 0.94 ; effBeam: 0.94 ; backend: NBC

v1: also extracts and writes the weather information
v2: offer output in EHT or GMVA format
v3 (2021.01.10): modified to be compatible with the new output (with more fields) of read_lastcal_info.py after EHT Data Management request in Dec. 2020
                 updates as well the condition to check that the four lines corresponf to the same scan (it was working OK before, but not checking all lines) 

P. Torne, IRAM 23.04.2018, last update 10.01.2021
"""

import sys, datetime, os
import numpy as np
import subprocess
import argparse

# Arguments control
parser = argparse.ArgumentParser()

parser.add_argument("-v", "--verbose", action="count", help="outputs detailed execution informaton")
#parser.add_argument("sched", help="name of the schedule from which calibration is read. Used to name the output file accordingly")
parser.add_argument("calinfo_file", help="input file where to read the calibration and weather information from (e.g. a Field System log). This code is written for Pv, and assumes that the text file contains rows with the 'rx' and 'wx' strings to fetch the calibration and weather data.")
#parser.add_argument("wxinfo_file", help="input file where to read the weather information from")
parser.add_argument("format", help="decides the output format: for GMVA or EHT data", choices=['GMVA', 'EHT'])

args = parser.parse_args()

#sched = args.sched.split(".vex")[0]  # Get the basename of the schedule name
basenm = os.path.basename( os.path.splitext(args.calinfo_file)[0] )

# Functions
def dewtemp(temp, humid):
    '''
    Calculate dewpoint following Thomas 'dewpoint.awk'.
    Calculate saturation vapor pressure(Es) and actual vapor pressure(E) in millibars
    http://www.gorhamschaffler.com/humidity_formulas.htm
    see also https://ag.arizona.edu/azmet/dewpoint.html
    # Formula verified twice. Results are correct (even if numbers look weird).
    '''
    ES = 6.11 * 10.0**(7.5*temp/(237.7+temp))
    E = humid * ES / 100.

    y = np.log(E/6.11) / 7.5 / np.log(10)

    dewtemp = 237.7*y/(1-y)

    # predictable water in mm
    pw = 432.98/(temp+273.16) * np.exp( (1.81 + 17.27 * dewtemp) / (dewtemp + 237.7) )

    return [dewtemp, pw]


# Fetch info from input file
try:

    calinfo = subprocess.check_output('grep %s -e "rx:"'%args.calinfo_file, shell=True).splitlines()
    wxinfo  = subprocess.check_output('grep %s -e "/wx/"'%args.calinfo_file, shell=True).splitlines()
    #wxinfo  = subprocess.check_output('grep %s -e "/wx/"'%args.wxinfo_file, shell=True).splitlines()

except Exception as e:

    print e
    print "\nHalting program."
    sys.exit(1)

# ********************************************
# Output CALIBRATION info in ANTAB format
# ********************************************

print "\nOpening file %s.antab to write ANTAB table ..."%basenm
outputfn = open("%s.antab"%basenm, "w")

if args.format=='EHT': # write header for EHT data (4 IFs)

    outputfn.write("TSYS PV  FT= 1.0  INDEX = 'R1:32', 'L1:32' ,'R1:32', 'L1:32' /\n")
    outputfn.write("!bands             HLI     VLI     HUI     VUI\n")
    outputfn.write("!DOY hh:mm:ss.ss   RCP     LCP     RCP     LCP     !  tau   elv   source\n")

elif args.format=='GMVA': # write header for GMVA data (2 IFs)

    outputfn.write("TSYS PV  FT= 1.0  INDEX = 'R1:8', 'L1:8'' /\n")
    outputfn.write("!bands             HLI     VLI     \n")
    outputfn.write("!DOY hh:mm:ss.ss   RCP     LCP     !  tau   elv   source\n")



all_timestamp = [] # always save last time stamp to avoid duplicates in the output file

print "Reading information from %d calibration scans from %s ..."%(len(calinfo)/4, args.calinfo_file)
for ii in range(0, len(calinfo), 4):  # each cal scan info comes in 4 rows for 4 IFs:
    # Read and format data:
    DOY       = []
    timestamp = []
    tsys      = []
    tau       = []
    elv       = []
    source    = []
    line      = []

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

        ts = data[5].split()[2]

        DOY.append(dayofyear)
        timestamp.append(ts)

        if jj == ii: # Ony append once to the all_timestamps vector
            all_timestamp.append(ts)

        # Extract Tsys values
        tempsys = float( data[1].split()[1] )
        # Check value
        if tempsys > 9999.9: 
            tsys.append(9999.9)
        elif (tempsys > 1 and tempsys <= 9999.9):
            tsys.append(tempsys)
        elif tempsys <= 0:
            tsys.append(np.NaN)

        # Extract tau
        tau.append( float(data[2].split()[1]) )

        try:

            # Extract elevation
            elv.append( float(data[8].split()[1]) ) # it was [7] before 2021 update of read_lastcal_info

            # Extract the source
            source.append(data[9].split()[1]) # it was [8] before 2021 update of read_lastcal_info

        except Exception as e:
            print "\n* Error: You seem to be passing a log from before Apr2018. Sorry, those are not in a compatible format.\n"
            sys.exit(1)
         

    if args.verbose >= 1:
        #print "DOY = %s"%DOY
        #print "timestamp = %s"%timestamp
        #print "tsys = %s"%tsys
        #print "source = %s"%source
        print "Extracting relevant information: DOY, timestamp, Tsys, tau, elev, source. for a cal scan on %s"%source[-1]

    # Check:
    if DOY[0:] == DOY[::-1] and timestamp[0:] == timestamp[::-1] and source[0:] == source [::-1]: # these are four lines corresponding to the same scan
        if args.verbose >= 1: print "Writing to calibration info on %s to ANTAB table ..."%source[-1]

        if ii == 0: # First entry, no duplicity check

            if args.format=="EHT":
                # Write one line with the four Tsys in four columns:
                outputfn.write("%d  %s.00 %7.1f %7.1f %7.1f %7.1f   !  %2.2f  %2.1f  %s\n"%(int(DOY[-1]), timestamp[-1], tsys[0], tsys[2], \
                                                                                tsys[1], tsys[3], np.mean((tau[0], tau[2], tau[1], tau[3])), elv[-1], source[-1]) )
            elif args.format=="GMVA":
                # Write one line with the two Tsys in two columns: (note that the IF order is different from GMVA and EHT! GMVA: E090H, E090V, E230H, E230V)
                outputfn.write("%d  %s.00 %7.1f %7.1f   !  %2.2f  %2.1f  %s\n"%(int(DOY[-1]), timestamp[-1], tsys[0], tsys[1], \
                                                                                     np.mean((tau[0], tau[1])), elv[-1], source[-1]) )

        elif ii > 0: # Ignore duplicity check for first entry

            if timestamp[-1] == all_timestamp[-2]:  # Do not write repeated cal info (occurs if no new cal is done when calling getcals in FS)
                print "timestamp = %s | previous = %s"%(timestamp[-1], all_timestamp[-2])
                print "Detected duplicated timestamp. Ignoring line."
                continue

            else:

                if args.format=="EHT":
                    # Write one line with the four Tsys in four columns:
                    outputfn.write("%d  %s.00 %7.1f %7.1f %7.1f %7.1f   !  %2.2f  %2.1f  %s\n"%(int(DOY[-1]), timestamp[-1], tsys[0], tsys[2], \
                                                                                       tsys[1], tsys[3], np.mean((tau[0], tau[2], tau[1], tau[3])), elv[-1], source[-1]) )
                elif args.format=="GMVA": 
                    # Write one line with the two Tsys in two columns: (note that the IF order is different from GMVA and EHT! GMVA: E090H, E090V, E230H, E230V)
                    outputfn.write("%d  %s.00 %7.1f %7.1f   !  %2.2f  %2.1f  %s\n"%(int(DOY[-1]), timestamp[-1], tsys[0], tsys[1], \
                                                                                          np.mean((tau[0], tau[1])), elv[-1], source[-1]) )
         

print "Writing to %s.antab ..."%basenm

    #if ii >= 4: break

# Write the 'end of file' marker for ANTAB
outputfn.write("/\n")

outputfn.close()

print "DONE."


# ********************************************
# Output WEATHER info in ANTAB format
# ********************************************

print "\nOpening file %s.wx to write ANTAB table ..."%basenm
outputfn = open("%s.wx"%basenm, "w")

# These parameters are written butnot needed for data processing, so set to 0.
wind=0.0
pa=0.0
rain=0.0
gust=0.0

# Header
outputfn.write("! ----- Weather information for Pv -----\n")
outputfn.write("! All values are instantaneous readings at the time indicated.\n")
outputfn.write("!              Temp   Press  DewPt   Wind Spd/Dir    Rain   Gust\n")
outputfn.write("!UT Day-Time    C      mBar    C      m/s    deg      cm     m/s\n")
outputfn.write("WEATHER PV /\n")

for ll in wxinfo:
    # Reformat the date/time:
    doy = int( ll.split()[0].split(".")[1] )
    time = ll.split()[0].split(".")[2]
    temp = float( ll.split()[1].split(",")[0] )
    press = float( ll.split()[2].split(",")[0] )
    humid = float( ll.split(",")[-1] )
    [dewt, pw ] = dewtemp(temp, humid)

    outputfn.write("%d-%s %5.1f   %6.1f %5.1f     %3.1f    %3.1f     %4.2f    %3.1f \n"%(
           doy, time, temp, press, dewt, wind, pa, rain, gust))


# Write the 'end of file' marker for ANTAB
outputfn.write("/\n")

outputfn.close()

print "DONE."

