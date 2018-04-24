#/usr/bin/env python
'''
Reads the calibration information from an input file,
and outputs it in ANTAB format for VLBI correlator.

v0: the cal info must be in the format output by
    read_lastcal_info.py (used in VLBI Field System), i.e.:

CAL info read from: /mrt-lx3/vis/vlbi/vlbireduc/calxmls/iram30m-calibration-NBC-20180421s111.xml
rx: E2HLI ; tsys: 581.70 ; tau: 0.40 ; pwv mm: 6.7 ; trx: 73.5 ; time: 2018-04-21 04:37:24 ; rxFreq: 214.8506 ; elev: 28.3 ; source: Mars
rx: E2HUI ; tsys: 611.80 ; tau: 0.41 ; pwv mm: 6.6 ; trx: 75.0 ; time: 2018-04-21 04:37:24 ; rxFreq: 227.3500 ; elev: 28.3 ; source: Mars
rx: E2VLI ; tsys: 605.85 ; tau: 0.38 ; pwv mm: 6.5 ; trx: 90.6 ; time: 2018-04-21 04:37:24 ; rxFreq: 214.8506 ; elev: 28.3 ; source: Mars
rx: E2VUI ; tsys: 629.44 ; tau: 0.40 ; pwv mm: 6.4 ; trx: 89.8 ; time: 2018-04-21 04:37:24 ; rxFreq: 227.3500 ; elev: 28.3 ; source: Mars

P. Torne, IRAM 23.04.2018
'''

import sys, datetime
import numpy as np
import subprocess
import argparse

# Arguments control
parser = argparse.ArgumentParser()

parser.add_argument("-v", "--verbose", action="store_true", help="outputs detailed execution informaton")
parser.add_argument("sched", help="name of the schedule without extension (e.g. e18c21). Used to name the output file accordingly")
parser.add_argument("calinfo_file", help="input file where to read the calibration information from")
parser.add_argument("wxinfo_file", help="input file where to read the weather information from")

args = parser.parse_args()

if args.verbose == True: verbose = True

outputfn = open("%spv.antab"%args.sched, "w")

# Load input file
try:

    calinfo = subprocess.check_output('grep %s -e "rx:"'%args.calinfo_file, shell=True).splitlines()
    wxinfo  = subprocess.check_output('grep %s -e "/wx/"'%args.wxinfo_file, shell=True).splitlines()

except Exception as e:

    print e
    print "\nHalting program."
    sys.exit(1)

# Output info in ANTAB format
outputfn.write("TSYS PV  FT= 1.0  INDEX = 'R1:32', 'L1:32' ,'R1:32', 'L1:32' /\n")
outputfn.write("!bands             HLI     VLI     HUI     VUI\n")
outputfn.write("!DOY hh:mm:ss.ss   RCP     LCP     RCP     LCP     !  tau   elv   source\n")

for ii in range(0, len(calinfo), 4):  # each cal scan info comes in 4 rows for 4 IFs:
    # Read and format data:
    DOY       = []
    timestamp = []
    tsys      = []
    tau       = []
    elv       = []
    source    = []
    line      = []

    print "ii = %d"%ii

    for jj in range(ii, ii+4): # Loop over four lines starting from ii

        data = calinfo[jj].split(" ; ")

        if verbose: print data

        # Extract DOY and timestamp
        yy = int( data[5].split()[1].split("-")[0] )
        mm = int( data[5].split()[1].split("-")[1] )
        dd = int( data[5].split()[1].split("-")[2] )
 
        dayofyear = datetime.datetime(yy, mm, dd).timetuple().tm_yday
        ts = data[5].split()[2]

        DOY.append(dayofyear)
        timestamp.append(ts)

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

        # Extract elevation
        elv.append( float(data[7].split()[1]) )

        # Extract the source
        source.append(data[8].split()[1])

        if verbose:
            print "ii=%d"%ii
            print "jj=%d"%jj
            #print "DOY Timestamp = %d %s.00 | Tsys = %s\n"%(DOY[-1], timestamp[-1], tsys[-1])
            print "DOY = %s"%DOY
            print "timestamp = %s"%timestamp
            print "tsys = %s"%tsys
            print "source = %s"%source

    # Check:
    if DOY[1:] == DOY[:-1] and timestamp[1:] == timestamp[:-1] and source[1:] == source [:-1]: # these are four lines corresponding to the same scan
        # Write one line with the four Tsys in four columns:
        outputfn.write("%d  %s.00 %7.1f %7.1f %7.1f %7.1f   !  %2.2f  %2.1f   %s\n"%(int(DOY[-1]), timestamp[-1], tsys[0], tsys[1], \
                                                                                   tsys[2], tsys[3], np.mean((tau[0], tau[1], tau[2], tau[3])), elv[-1], source[-1]) )
         

    #if ii >= 4: break

# Write the 'end of file' marker for ANTAB
outputfn.write("/\n")

outputfn.close()
