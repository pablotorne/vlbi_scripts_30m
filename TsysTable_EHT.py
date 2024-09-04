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
v3: - Handles sources observed during VLBI but without any CAL. info
    - Cross-check source names only with first 8 characters to deal with
      sources that appear with a longer name in a CAL scan and a short
      version in the .prn (e.g., j1924-2914 and J1924-29, e21a14 - Apr21)
v4: Need to change some lines to be compatible with the new scan ID format
    introduced in March 2022: DOY-HHMM of start of scan, instad of NoXXXX.
 
P. Torne, IRAM v2: 03.01.2020, v3: 2021.06.05, v4: 2022.07.22
"""

import sys, datetime, os
import numpy as np
import subprocess
import argparse
import matplotlib.pyplot as plt
import re

# Global parameter:
DELTA_TIME = 600 # seconds, or less, to consider a cal. scan info OK to be assigned to a VLBI scan

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
    #prninfo  = subprocess.check_output('grep -E "no[0-9]{1,}" %s'%args.prn_file, shell=True).splitlines()
    # In 2022, they changed the format of the scans ID from no####, to DOY-HHMM of scan start. Need code update:
    prninfo  = subprocess.check_output('grep -E "  :  :  " %s'%args.prn_file, shell=True).splitlines()

except Exception as e:

    print e
    print "\nHalting program."
    sys.exit(1)

# Support functions
def extractTrackInfo(prnfile):
    '''
    Extract Start Date (UT), Experiment name, Station ID
    from PRN file to add to the Tsys* table header
    '''

    print "\nReading information about VLBI scans from %s"%args.prn_file

    # Open and read in all the text file
    f = open(prnfile, 'r')
    expinfo = f.readlines()
    f.close()

    # Search for the first "date =" entry, which gives the Start Date (UT) of the experiment
    for line in expinfo:
        linedata = line.split()
        if args.verbose >= 3:
            print linedata
        if len(linedata) > 0 and linedata[0] == "date" and linedata[1] == "=":
            if args.verbose: print "Found the Start Date from the prn file: %s"%linedata[2]
            startdate = linedata[2]
            break

    # We need two separate loops because the first one must stop after finding the first Date 
    # (more than one Date may be found in one schedule if it goes over 00:00 UTC)
    for line in expinfo:
        linedata = line.split()
        if args.verbose >= 3:
            print linedata
        if len(linedata) > 0 and linedata[0] == "Station:":
            if args.verbose: print "Found the Track name and Station ID: %s, %s"%(linedata[5], linedata[2].lstrip("(").rstrip(")") )
            track = linedata[5]
            stationcode = linedata[2].lstrip("(").rstrip(")")
            break

    return [startdate, track, stationcode]

def find_idx(lst, a):
    '''
    Will return the indices of a list when the values
    match the parameter "a".
    Case insensitive!
    '''
    result = []
    for i, x in enumerate(lst):
        # Truncate 8 chars of the names from the CAL. (Apr2021 - to deal with source with different names in the .prn and the CAL info)
        #print "x.lower() = %s"%x.lower()[:8]
        #print "a.lower() = %s"%a.lower()
        if x.lower()[:8] == a.lower():
            result.append(i)
    return result

def datetime_to_float(d):
    '''
    Convert datetime to float with millisecond precision
    https://stackoverflow.com/questions/35337299/python-datetime-to-float-with-millisecond-precision
    '''
    epoch = datetime.datetime.utcfromtimestamp(0)
    total_seconds =  (d - epoch).total_seconds()
    # total_seconds will be in decimals (millisecond precision)
    return total_seconds


def nearest(items, pivot):
    '''
    Find nearest between two objects of any type supporting comparison.
    Used to find the nearest time of a calibration to a VLBI scan.
    Adapted from
     https://stackoverflow.com/questions/32237862/find-the-closest-date-to-a-given-date/32237949
    -> This function will return the datetime in items which is the closest to the date pivot.
    --> So, order of input parameters is important
    - Modified: also returns delta_t and index of items for nearest element
    '''
    
    delta_dt = []

    for item in items:
        delta_dt.append( abs(item - pivot) )

    delta_t         = min(delta_dt)
    nearest_cal_idx = np.argmin(delta_dt)
    nearest_cal_ts  = items[nearest_cal_idx]

    #return min(items, key=lambda x: abs(x - pivot))

    # TBD: Return the nearest cal. timestamps BEFORE and AFTER
    #      the Pivot. This can be used to make more efficient
    #      the Interpolation later by using only 2 closest values.

    return [nearest_cal_ts, nearest_cal_idx, delta_t]

def checkScanID(line):
    '''
    Check the initial part of the passed line using regular expressions
    to identify which lines correspond to a VLBI scan in the format
    of the .prn from March 2022
    Adapted from https://stackoverflow.com/questions/14966647/check-python-string-format
    '''
    r = re.compile('\d{3}-\d{4}')
    print "line =%s"%str(line)
    # Compare only with the first characters of the line:
    if r.match( str(line) ) is not None:
        print "True"
        return True
    else:
        print "False"
        return False    




# Get the Start Time (UT), Track name and Station Code from the PRN file
startdate, track, stationcode = extractTrackInfo(args.prn_file)

# *************************************************
# Output CALIBRATION info in EHT Tsys Table format
# *************************************************

print "\nOpening file %s_PV.CAL.txt to write Tsys and Cal. Metadata Summary table ..."%basenm
calsummary_table = open("%s_CALs_%s.txt"%(track, stationcode), "w")

# Write the header of the Tsys table file:
calsummary_table.write("################################################################\n")
calsummary_table.write("#\n")
calsummary_table.write("# CALIBRATION DATA SUMMARY from PV for track: %s\n"%track)
calsummary_table.write("#\n")
calsummary_table.write("################################################################\n")
calsummary_table.write("# * NOTE: This is NOT the Tsys table! (that's a different file)\n")
calsummary_table.write("#\n")
calsummary_table.write("## Timestamp(UT)       Scan	Source Pos_Az Pos_El  Tsys_b1r\t_b1l\t_b2r\t_b2l\t_b3r\t_b3l\t_b4r\t_b4l\tTau\tTamb\tTatm\n")
calsummary_table.write("# YYYY-MM-DD HH:MM:SS  (VEX)	       (deg)  (deg)     (K)\t(K)\t(K)\t(K)\t(K)\t(K)\t(K)\t(K)\t(zen)\t(K)\t(K)\n")
calsummary_table.write("################################################################################################################################################################\n")

# Keep all the calibration data in lists to later operate with them, for instance
# to search for nearest calibration and be able to interpolate if necessary
# NOTE: Python lists are ORDERED, so the [0] element of all these lists corresponds to each other
all_timestamps = [] # always save last time stamp to avoid duplicates in the output file
all_sources   = []
all_azimuth   = []  # 
all_elevation = []  # can be used to do a sanity check of the data selected
all_tsys_b1r  = []  # E2HLI (Cal. Line 0 in output of read_lastcal_info.py with usual EMIR setup)
all_tsys_b1l  = []  # E2VLI (Cal. Line 2 in output of read_lastcal_info.py with usual EMIR setup)
all_tsys_b2r  = []  # E2HLI (Cal. Line 0 in output of read_lastcal_info.py with usual EMIR setup)
all_tsys_b2l  = []  # E2VLI (Cal. Line 2 in output of read_lastcal_info.py with usual EMIR setup)
all_tsys_b3r  = []  # E2HUI (Cal. Line 1 in output of read_lastcal_info.py with usual EMIR setup)
all_tsys_b3l  = []  # E2ULI (Cal. Line 3 in output of read_lastcal_info.py with usual EMIR setup)
all_tsys_b4r  = []  # E2HUI (Cal. Line 1 in output of read_lastcal_info.py with usual EMIR setup)
all_tsys_b4l  = []  # E2ULI (Cal. Line 3 in output of read_lastcal_info.py with usual EMIR setup)
all_tauz      = []
all_Tamb      = []
all_Tatm      = []
all_RxFreq    = []

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

        if args.verbose >= 3: print "Reading line: %s"%data

        # Extract DOY and timestamp
        yy = int( data[5].split()[1].split("-")[0] )
        mm = int( data[5].split()[1].split("-")[1] )
        dd = int( data[5].split()[1].split("-")[2] )
 
        dayofyear = datetime.datetime(yy, mm, dd).timetuple().tm_yday
        
        if args.verbose >= 3:
            print "Date from cal data: %s-%s-%s"%(yy, mm, dd)
            print "DOY = datetime.datetime(yy, mm, dd).timetuple().tm_yday = %d"%dayofyear 

        ts = data[5].split("time:")[-1].lstrip()

        DOY.append(dayofyear)
        timestamp.append(ts)

        if jj == ii: # Only append once to the all_timestamps vector
            all_timestamps.append(ts)

        # Extract rxFreq
        rxfrequency = float( data[6].split()[1] )
        rxFreq.append( rxfrequency )



        # Extract Tsys values
        tempsys = float( data[1].split()[1] )
        print "tempsys = %s"%tempsys
        # Check value
        if tempsys > 9999.9: 
            tsys_star.append(9999.9)
        elif (tempsys > 1 and tempsys <= 9999.9):
            tsys_star.append(tempsys)
        elif tempsys <= 0:
            tsys_star.append(-1)
        else:
            tsys_star.append(-1)

        # Extract tau
        tau_z.append( float(data[2].split()[1]) )

        try:

            # Extract elevation
            azimuth.append( float(data[7].split()[1]) )

            if jj == ii: # Only append once per entry to the azimuth vector
                all_azimuth.append(azimuth[-1])

            # Extract elevation
            elevation.append( float(data[8].split()[1]) )

            if jj == ii: # Only append once per entry to the elevation vector
                all_elevation.append(elevation[-1])

            # Extract the source
            source.append( data[9].split()[1] )

            if jj == ii: # Only append once per entry to the sources vector
                all_sources.append(source[-1])

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

    if args.verbose >= 2:
        #print "DOY = %s"%DOY
        #print "timestamp = %s"%timestamp
        #print "tsys = %s"%tsys
        #print "source = %s"%source
        print "Extracting relevant information from a cal scan on %s"%source[-1]

    # Check:
    if DOY[0:] == DOY[::-1] and timestamp[0:] == timestamp[::-1] and source[0:] == source [::-1]: # these are four lines corresponding to the same scan

        # Save all the calibration data into master lists, used later to relate cal.<-> VLBI scans, and interpolate values if needed
        if args.verbose >= 2: print "Saving calibration data on source %s into internal Python lists ..."%source[-1]

        # all_timestamps, all_sources, all_azimuth, all_elevation are appended in the previous loop, 
        # to make sure we only append once per cal. entry

        # Append to the corresponding all_tsys_xxx depending of Rx band:
        all_tsys_b1r.append( tsys_star[0] )
        all_tsys_b1l.append( tsys_star[2] )
        all_tsys_b2r.append( tsys_star[0] )
        all_tsys_b2l.append( tsys_star[2] )

        all_tsys_b3r.append( tsys_star[1] )
        all_tsys_b3l.append( tsys_star[3] )
        all_tsys_b4r.append( tsys_star[1] )
        all_tsys_b4l.append( tsys_star[3] )

        # Append the average Tauz, Tamb and Tatm_s
        all_tauz.append( np.mean(tau_z) )
        all_Tamb.append( np.mean(Tamb)  )
        all_Tatm.append( np.mean(Tatm_s)  )

        # Write all the calibration information formatted as EHT into an output Cal. Summary Table
        if args.verbose >= 2: print "Writing extracted calibration info on %s to CAL. summary table ..."%source[-1]

        if ii == 0: # First entry, no duplicity check

            # Write one line with the requested metadata:
            # Note: at 64Gbps, we only record LowerInner and UpperInner
            # In EHT/Correlator nomenclature: Tsys* of LI = Tsys* of Band1 *and* Band2 (we do not measure Tsys separately for them)
            #                                 Tsys* of UI = Tsys* of Band3 *and* Band4 (and H=RCP and V=LCP)
            calsummary_table.write("%s\tCALXML\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.1f\t%.1f\n"%(\
                           timestamp[-1], source[-1], azimuth[-1], elevation[-1], \
                           tsys_star[0], tsys_star[2], tsys_star[0], tsys_star[2], \
                           tsys_star[1], tsys_star[3], tsys_star[1], tsys_star[3], \
                           np.mean((tau_z[0], tau_z[2], tau_z[1], tau_z[3])), \
                           np.mean((Tamb[0], Tamb[2], Tamb[1], Tamb[3])), \
                           np.mean((Tatm_s[0], Tatm_s[2], Tatm_s[1], Tatm_s[3])), \
                           ))

        elif ii > 0: # Ignore duplicity check for first entry

            if timestamp[-1] == all_timestamps[-2]:  # Do not write repeated cal info (occurs if no new cal is done when calling getcals in FS)
                print "timestamp = %s | previous = %s"%(timestamp[-1], all_timestamps[-2])
                print "Detected duplicated timestamp. Ignoring line."
                continue

            else:
                calsummary_table.write("%s\tCALXML\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.1f\t%.1f\n"%(\
                           timestamp[-1], source[-1], azimuth[-1], elevation[-1], \
                           tsys_star[0], tsys_star[2], tsys_star[0], tsys_star[2], \
                           tsys_star[1], tsys_star[3], tsys_star[1], tsys_star[3], \
                           np.mean((tau_z[0], tau_z[2], tau_z[1], tau_z[3])), \
                           np.mean((Tamb[0], Tamb[2], Tamb[1], Tamb[3])), \
                           np.mean((Tatm_s[0], Tatm_s[2], Tatm_s[1], Tatm_s[3])), \
                           ))


print "Writing to %s_PV.CAL.txt ..."%basenm

    #if ii >= 4: break

calsummary_table.close()
print "Finished writing summary calibration table and saving metadata in memory.\n\n"
# The summary calibration table is an extra, not what they as for.

# We will work now from the PRN file, which gives us the VLBI scan information
# and the Python lists in which we have all the calibration data. From here:
# For each VLBI scan,
# 0) Extract a subset of all the calibration scans only on that source
# 1) Find the nearest calibration on the same source in time, if the cal. scan
#    it taken within delta_t (typically 10 minutes), we use that information
#    to populate the Tsys* table
# 2) If delta_t between the cal. scan and the VLBI scan START TIME is larger,
#    make a linear interpolation on the calibration data and use those to get
#    values at the VLBI start time and populate the Tsys table with them

# Option to improve: ALWAYS interpolate the cal. data and derive values at the
# [TDB]              middle of the VLBI scan time

# REMOVE?
### Get loaded the CAL data properly formatted (note: the file must not contain any header):
#print("Opening and reading in summary and formatted calibration table from %s_PV.CAL.txt\n"%basenm)
#g = open("%s_PV.CAL.txt"%basenm, 'r')
#CALDATAFRMT = g.readlines()
#g.close() 
# Here we should open a new file and write the heade of thw final Tsys table for this track
# BLABLA

# Convert the all_timestamps list of the calibration data time stamps
# to datetime objects to later compare easily with the time stamps of
# the VLBI scans:
all_timestamps_dt = []
if args.verbose > 3: print("Creating datetime objects for cal. time stamps...")
for tstamp in all_timestamps:
    tstamp_dt = datetime.datetime.strptime(tstamp, '%Y-%m-%d %H:%M:%S')
    all_timestamps_dt.append( tstamp_dt )

# Let's loop over ALL the PRN file. For each VLBI scan entry do:
# - Create a proper time stamp with date and time in datetime format
# - Compare with the calibration timestamps and find the CAL. data
#   that corresponds to the VLBI scan
# - Write out the corresponding row in the Tsys table

print("Reading the PRN file and writing out calibration metadata ...\n")

# Open the PRN file and read in
f = open(args.prn_file, 'r')
prntext = f.readlines()
f.close()

# Open the output text file that will be the Tsys metadata table:

tsys_table = open("%s_%s.tsys.txt"%(track, stationcode), 'w')

# Write the header of the Tsys table file:
tsys_table.write("################################################################\n")
tsys_table.write("#\n")
tsys_table.write("# Tsys/Tsys* table for each station and track\n")
tsys_table.write("#\n")
tsys_table.write("################################################################\n")
tsys_table.write("#\n")
tsys_table.write("# Track Start Date (UT): %s\n"%datetime.datetime.strptime(startdate, '%Y%b%d').date())
tsys_table.write("# Experiment name: %s\n"%track)
tsys_table.write("# Station ID: %s\n"%stationcode)
tsys_table.write("# Operators/observer/contact person: Pablo Torne (torne@iram.es)\n")
tsys_table.write("# Tsys* or Tsys (inclusive of opacity or not): Tsys*\n")
#tsys_table.write("# Central observing frequencies in GHz for Tsys measurement (DSB/2SB): %.1f (b1 & b2), %.1f (b3 & b4) (2SB)\n"%(rxFreq[0]+0.25, rxFreq[1]-0.25))
tsys_table.write("# Central observing frequencies in GHz for Tsys measurement (DSB/2SB): %.1f (b1 & b2), %.1f (b3 & b4) (2SB)\n"%(rxFreq[0], rxFreq[1])) #Apr21, remove 0.25 shift, those affect the VLBI recording but should not affect the NBC data
tsys_table.write("# Sideband ratio (if available, for DSB): NA\n")
tsys_table.write("# Sideband coupling coefficient (if available, for 2SB): %.3f\n"%np.mean(gainImage))
tsys_table.write("# Tau observing frequency in GHz: %.1f (average of central observing frequencies)\n"%np.mean(rxFreq))
tsys_table.write("#\n")
tsys_table.write("################################################################\n")
tsys_table.write("# * Timestamp gives the date and time at which Tsys*/Tsys was measured (can be during, before or after each scan)\n")
tsys_table.write("# * Scan column gives the VEX scan number associated with each Tsys*/Tsys measurement\n")
tsys_table.write("# * Tau values should be that measured at zenith\n")
tsys_table.write("# * Empty cells should have an 'NA'\n")
tsys_table.write("# * If the same Tsys/Tsys* values are used for multiple bands, then include repeated values in the appropriate columns\n")
tsys_table.write("# * Tsys*/Tsys values should be provided without sideband correction\n")
tsys_table.write("# * Column Delta_t indicates the absolute time difference between the calibration data and the start time of the VLBI scan.\n")
tsys_table.write("# * If Delta_t > 10 min, we interpolate between the closest calibration values for the source (indicated with Interp.).\n")
tsys_table.write("# * If there is no calibration data for a given scan / source, columns indicate NA and Column Delta_t indicates NOCAL_INFO.\n") 
tsys_table.write("#\n")
tsys_table.write("## Timestamp(UT)       Scan	       Source Pos_Az Pos_El  Tsys_b1r\t_b1l\t_b2r\t_b2l\t_b3r\t_b3l\t_b4r\t_b4l\tTau\tTamb\tTatm\tDelta_t\n")
tsys_table.write("# YYYY-MM-DD HH:MM:SS  (VEX)	              (deg)  (deg)     (K)\t(K)\t(K)\t(K)\t(K)\t(K)\t(K)\t(K)\t(zen)\t(K)\t(K)\t(min)\n")
tsys_table.write("################################################################################################################################################################\n")

# Loop over PRN in search of VLBI scans, work out cal. data, and populate Tsys table:

# Save some numbers for statistics:
direct_assignations = 0
interpolations      = 0
no_cal_data         = 0 

for line in prntext:
    
    if args.verbose > 3: print "line from prntext = %s"%line

    # Detect if the line is a "date =" entry or a VLBI scan
    lineinfo = line.split()
    
    if len(lineinfo) > 0 and lineinfo[0] == "date" and lineinfo[1] == "=": # it is a "date =" entry
        # Save the date in UTC:
        if args.verbose > 2: print("Found a new date line in the PRN file! Saving date to format time stamps.")
        datestr = lineinfo[2] # in str() format

    #elif len(lineinfo) > 0 and lineinfo[0][0:2] == "no" and ( len(lineinfo) == 11 or len(lineinfo) == 12):
    # Previous line needs update to be compatible with the new format of scan ID in the .prn as from March 22 (DOY-HHMM)
    elif len(lineinfo) > 0 and checkScanID(lineinfo[0]) == True and ( len(lineinfo) == 11 or len(lineinfo) == 12): 
        # This is a row with a VLBI scan
        vlbi_scannumber = lineinfo[0]
        source_name     = lineinfo[2]
        az              = float(lineinfo[3])
        el              = float(lineinfo[4])
        if len(lineinfo) == 11: # this is the first VLBI scan of the prn, has one colunm less
            starttime       = lineinfo[6]
        elif len(lineinfo) > 11: # rest of rows
            starttime       = lineinfo[7]
        else:
            print("\n\n\n ERROR when reading VLBI scan info from PRN file. Aborting!")
            sys.exit(1)

        if args.verbose >= 1: print("Found VLBI scan %s in the PRN file (on %s). Matching with CAL. data in memory ..."%(vlbi_scannumber, source_name)) 

        # Convert startime of VLBI scan to a datetime object that
        # we can use to calculate time deltas nicely
        timestamp_frmt = "%s %s"%(datestr, starttime)
        if args.verbose > 2: print "Formatted VLBI scan Date+Time: %s"%timestamp_frmt
        timestamp_vlbi_dt = datetime.datetime.strptime(timestamp_frmt, '%Y%b%d %H:%M:%S')
        if args.verbose > 2: print("Time stamp as datetime object %s"%timestamp_vlbi_dt)
        if args.verbose > 1: print("Obtaining Cal. info Metadata for VLBI Scan %s on %s at %s"%(\
                                 vlbi_scannumber, source_name, timestamp_vlbi_dt))
        if args.verbose > 2: 
            print("Extracted info from the PRN for this VLBI scan:")
            print("%s\t%s\t%s\t%s\t%s:"%(vlbi_scannumber, source_name, az, el, timestamp_vlbi_dt))

        # KEY STEP: Extract the subset of all calibration data on THIS source only:
        source_idx   = find_idx(all_sources, source_name) # extract indices of all_sources for those entries of 'source_name'

        # In April 2021 we found cases in which a source appears with a different name in the .prn and in the CAL data
        # E.g. j1924-2914 in the CALXMLs, and J1924-29 in the VLBI schedules/.prn files
        # We need a clever way to cross-match these cases!
        # Try 1: Truncate first 8 characters of the longest of the names, make all small letters, then compare:
        

        source_ts_dt = [all_timestamps_dt[k] for k in source_idx] # timestamps of cal. scans on this source
        source_calname = [all_sources[k] for k in source_idx]     # the cal. scan source
        source_az       = [all_azimuth[k] for k in source_idx]
        source_el       = [all_elevation[k] for k in source_idx]
        source_tsys_b1r = [float(all_tsys_b1r[k]) for k in source_idx]   # Tsys*_b1r of cal. scans on this source
        source_tsys_b1l = [all_tsys_b1l[k] for k in source_idx]   # Tsys*_b1l of cal. scans on this source
        source_tsys_b2r = [all_tsys_b2r[k] for k in source_idx]   # Tsys*_b2r of cal. scans on this source
        source_tsys_b2l = [all_tsys_b2l[k] for k in source_idx]   # Tsys*_b2l of cal. scans on this source
        source_tsys_b3r = [all_tsys_b3r[k] for k in source_idx]   # Tsys*_b3r of cal. scans on this source
        source_tsys_b3l = [all_tsys_b3l[k] for k in source_idx]   # Tsys*_b3l of cal. scans on this source
        source_tsys_b4r = [all_tsys_b4r[k] for k in source_idx]   # Tsys*_b4r of cal. scans on this source
        source_tsys_b4l = [all_tsys_b4l[k] for k in source_idx]   # Tsys*_b4l of cal. scans on this source
        source_tauz     = [all_tauz[k] for k in source_idx]
        source_Tamb     = [all_Tamb[k] for k in source_idx]
        source_Tatm     = [all_Tatm[k] for k in source_idx]

        if args.verbose >= 3:
            print("Extracting the subset of calibration data in %s for source %s ..."%(args.calinfo_file, source_name))
            print("source_idx = %s"%source_idx)
            print("source_ts_dt = %s"%source_ts_dt)
            print("source_calname = %s"%source_calname)
            print("source_tsys_b1r = %s"%source_tsys_b1r)
            print("source_tsys_b1l = %s"%source_tsys_b1l)
            print("source_tsys_b2r = %s"%source_tsys_b2r)
            print("source_tsys_b2l = %s"%source_tsys_b2l)
            print("source_tsys_b3r = %s"%source_tsys_b3r)
            print("source_tsys_b3l = %s"%source_tsys_b3l)
            print("source_tsys_b4r = %s"%source_tsys_b4r)
            print("source_tsys_b4l = %s"%source_tsys_b4l)
            print("source_tauz = %s"%source_tauz)
            print("source_Tamb = %s"%source_Tamb)
            print("source_Tatm = %s"%source_Tatm)

        # Find nearest calibration time stamp to the VLBI scan time stamp from the subset just extracted,
        # i.e., from the cal. data on this source:

        try:
            nearest_cal_ts_dt, nearest_cal_index, delta_t = nearest(source_ts_dt, timestamp_vlbi_dt)

            print("Nearest-in-time calibration scan found for VLBI scan %s on %s at %s is:"%(vlbi_scannumber, source_name, timestamp_vlbi_dt))
            closest_cal_info = "%s\tCALXML\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.1f\t%.1f\n"%(\
                             source_ts_dt[nearest_cal_index], source_calname[nearest_cal_index],\
                             source_az[nearest_cal_index], source_el[nearest_cal_index],\
                             source_tsys_b1r[nearest_cal_index], source_tsys_b1l[nearest_cal_index],\
                             source_tsys_b2r[nearest_cal_index], source_tsys_b2l[nearest_cal_index],\
                             source_tsys_b3r[nearest_cal_index], source_tsys_b3l[nearest_cal_index],\
                             source_tsys_b4r[nearest_cal_index], source_tsys_b4l[nearest_cal_index],\
                             source_tauz[nearest_cal_index], source_Tamb[nearest_cal_index], source_Tatm[nearest_cal_index])
            print(closest_cal_info)

            # Quick sanity check: double check that cal. scan source == VLBI scan science source:

            # From Apr2021 we check the coincidende of first 9 characters
            if source_calname[nearest_cal_index].lower()[:8] != source_name.lower()[:8]: # case insensitive check
                print("ERROR: Cal. scan source != from VLBI scan source! Halting program.")
                sys.exit(1)
               

            # Strategy 1: Check is delta_t between cal scan and VLBI scan < DELTA_TIME
            if delta_t.seconds <= DELTA_TIME:
                print("This cal. was taken within %.1f minutes of the VLBI scan start time. Accepting cal. metadata -> Writing out to Tsys table..."%(delta_t.seconds/60.))

                tsys_table.write("%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.1f\t%.1f\t%.1f\n"%(\
                             source_ts_dt[nearest_cal_index], vlbi_scannumber, source_calname[nearest_cal_index],\
                             source_az[nearest_cal_index], source_el[nearest_cal_index],\
                             source_tsys_b1r[nearest_cal_index], source_tsys_b1l[nearest_cal_index],\
                             source_tsys_b2r[nearest_cal_index], source_tsys_b2l[nearest_cal_index],\
                             source_tsys_b3r[nearest_cal_index], source_tsys_b3l[nearest_cal_index],\
                             source_tsys_b4r[nearest_cal_index], source_tsys_b4l[nearest_cal_index],\
                             source_tauz[nearest_cal_index], source_Tamb[nearest_cal_index], source_Tatm[nearest_cal_index],\
                             delta_t.seconds/60.)
                            )

                direct_assignations += 1

            # Strategy 2: If the time between the cal. scan and the VLBI start time > DELTA_TIME, interpolate between closest cal. scans:
            else:
                print("Cal. data Delta_t to VLBI scan > %.1f minutes. Interpolating linearly using the closest cal. scans on this source ..."%(delta_t.seconds/60.))

                # To interpolate we need to convert the datetime objects of timestamps to float numbers:
                source_ts_floats = [datetime_to_float(source_ts_dt[k]) for k in range(len(source_ts_dt))]
                vlbiscan_ts_float = datetime_to_float(timestamp_vlbi_dt)
        
                # Create interpolation function:

                # Note: this would be faster if we only interpolate 2-nearest values,...
                #       perhaps nearest could return a SORTED array of cals, ordered by delta_t
                #       but we need to select only the one before and after the VLBI scan.
                #       For the moment interpolate all the array, it is fast.
           
                # Tidy up a bit to allow for the use of a for loop for Tsys interp 
                source_tsys = [source_tsys_b1r, source_tsys_b1l, source_tsys_b2r, source_tsys_b2l, \
                                source_tsys_b3r, source_tsys_b3l, source_tsys_b4r, source_tsys_b4l]            
                tsys_interp = []

                for band in range(8): # We are currently working with 4 bands x 2 polarisations

                    x = np.array(source_ts_floats)
                    y = np.array(source_tsys[band])

                    tsys_interp_value = np.interp(vlbiscan_ts_float, x, y, left=y[0], right=y[-1])

                    # With params left and right we take care of the rare cases in which we have
                    # a VLBI scan that is > DELTA_TIME from the nearest cal. scan and that IS NOT
                    # sorrounded by cal. scans to allow the interpolation with existing data.
                    # This can occur e.g., at the very beginning or end of a track / schedule if
                    # a cal. is not done on the source close enough to the VLBI scan start time.

                    tsys_interp.append( tsys_interp_value )

                # Interpolate azimuth:
                y = np.array(source_az)
                az_interp = np.interp(vlbiscan_ts_float, x, y, left=y[0], right=y[-1])

                # Interpolate elevation:
                y = np.array(source_el)
                el_interp = np.interp(vlbiscan_ts_float, x, y, left=y[0], right=y[-1])

                # Interpolate tauz:
                y = np.array(source_tauz) 
                tauz_interp = np.interp(vlbiscan_ts_float, x, y, left=y[0], right=y[-1])

                # Interpolate Tamb:
                y = np.array(source_Tamb) 
                Tamb_interp = np.interp(vlbiscan_ts_float, x, y, left=y[0], right=y[-1])

                # Interpolate Tatm:
                y = np.array(source_Tatm)
                Tatm_interp = np.interp(vlbiscan_ts_float, x, y, left=y[0], right=y[-1])

                # Write to Tsys table:
                tsys_table.write("%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.1f\t%.1f\t%s\n"%(\
                                 timestamp_vlbi_dt, vlbi_scannumber, source_calname[nearest_cal_index],\
                                 az_interp, el_interp,\
                                 tsys_interp[0], tsys_interp[2], tsys_interp[0], tsys_interp[2],\
                                 tsys_interp[1], tsys_interp[3], tsys_interp[1], tsys_interp[3],\
                                 tauz_interp, Tamb_interp, Tatm_interp,\
                                 "Interp")
                                )

                interpolations += 1

                if args.verbose >= 4:
                    # Make a plot with the interpolated data, and stop the code to avoid creating tens of plots?
                    plt.figure()
                    plt.plot(x, np.array(source_tsys[band]), '.', markersize=15)
                    plt.plot(x, np.array(source_tsys[band]), '-')
                    plt.plot(vlbiscan_ts_float, tsys_interp[band], marker="*", markersize=15)
                    plt.show()
                    #sys.exit(0)



            print("\n")

        except:  # THERE IS NO CALIBRATION DATA FOR THIS SOURCE FOUND. IT CAN HAPPEN! # WRITE "NA" in the columns

            print "\n WARNING: There is NO CALIBRATION DATA found for this source. I cannot produce cal. data for VLBI scan %s\n"%vlbi_scannumber
            print "Adding NA in the result columns for this entry.\n"
            #sys.exit(1)

            tsys_table.write("%s\t%s\t%s\t%.1f\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(\
                             timestamp_vlbi_dt, vlbi_scannumber, source_name,\
                             az, el,\
                             "NA", "NA",\
                             "NA", "NA",\
                             "NA", "NA",\
                             "NA", "NA",\
                             "NA", "NA", "NA",\
                             "NOCAL_INFO")
                            )

            no_cal_data += 1
    
# Close output file          
tsys_table.close()

print "VLBI scans with direct cal. data assigned = %d"%direct_assignations
print "VLBI scans requiring interpolating = %d"%interpolations
print "VLBI scans with NO CAL. DATA = %d"%no_cal_data

print "DONE."
