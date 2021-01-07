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
from scipy.interpolate import interp1d
import subprocess
import argparse

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
    prninfo  = subprocess.check_output('grep -E "no[0-9]{1,}" %s'%args.prn_file, shell=True).splitlines()

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

def singlesource_callist_timestamps_dt(timestamp_list_dt, CALDATAFRMT, source):
    '''
    Loop over our master vector of calibration time stamps and
    save to a temporal vector only those time stamps that
    correspond to calibrations on the source.
    We need to use the CAL. LIST SUMMARY to cross-check sources
    in the calibration scans!
    '''

    singlesource_timestamp_list_dt = []

    for ii in range( len(timestamp_list_dt) ):
        if CALDATAFRMT[ii].split("\t")[2].lower() == source_name.lower():
            if args.verbose > 2: print("Saving a timestamp_dt from a cal. on %s (for VLBI scan on %s)"%(CALDATAFRMT[ii].split("\t")[2].lower(), source_name.lower()))
            singlesource_timestamp_list_dt.append(timestamp_list_dt[ii])

    return singlesource_timestamp_list_dt


def find_idx(lst, a):
    '''
    Will return the indices of a list wher the values
    match the parameter "a".
    Case insensitive!
    '''
    result = []
    for i, x in enumerate(lst):
        if x.lower() == a.lower():
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


# Get the Start Time (UT), Track name and Station Code from the PRN file
startdate, track, stationcode = extractTrackInfo(args.prn_file)

# *************************************************
# Output CALIBRATION info in EHT Tsys Table format
# *************************************************

print "\nOpening file %s_PV.CAL.txt to write Tsys and Cal. Metadata Summary table ..."%basenm
outputfn = open("%s_PV.CAL.txt"%basenm, "w")

# Write the header of the Tsys table file:
#outputfn.write("################################################################\n")
#outputfn.write("#\n")
#outputfn.write("# Tsys/Tsys* table for each CALIBRATION scan from PV: %s\n"%basenm)
#outputfn.write("#\n")
#outputfn.write("################################################################\n")

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

        if jj == ii: # Ony append once to the all_timestamps vector
            all_timestamps.append(ts)

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
        if args.verbose >= 1: print "Saving calibration data on source %s into internal Python lists ..."%source[-1]

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
        if args.verbose >= 1: print "Writing extracted calibration info on %s to CAL. summary table ..."%source[-1]

        if ii == 0: # First entry, no duplicity check

            # Write one line with the requested metadata:
            # Note: at 64Gbps, we only record LowerInner and UpperInner
            # In EHT/Correlator nomenclature: Tsys* of LI = Tsys* of Band1 *and* Band2 (we do not measure Tsys separately for them)
            #                                 Tsys* of UI = Tsys* of Band3 *and* Band4 (and H=RCP and V=LCP)
            outputfn.write("%s\tCALXML\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.3f\t%.1f\t%.1f\n"%(\
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
                outputfn.write("%s\tCALXML\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.3f\t%.1f\t%.1f\n"%(\
                           timestamp[-1], source[-1], azimuth[-1], elevation[-1], \
                           tsys_star[0], tsys_star[2], tsys_star[0], tsys_star[2], \
                           tsys_star[1], tsys_star[3], tsys_star[1], tsys_star[3], \
                           np.mean((tau_z[0], tau_z[2], tau_z[1], tau_z[3])), \
                           np.mean((Tamb[0], Tamb[2], Tamb[1], Tamb[3])), \
                           np.mean((Tatm_s[0], Tatm_s[2], Tatm_s[1], Tatm_s[3])), \
                           ))


print "Writing to %s_PV.CAL.txt ..."%basenm

    #if ii >= 4: break

outputfn.close()
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
print("Opening and reading in summary and formatted calibration table from %s_PV.CAL.txt\n"%basenm)
g = open("%s_PV.CAL.txt"%basenm, 'r')
CALDATAFRMT = g.readlines()
g.close() 
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

print("Reading the PRN file and writing out calibration metadata ...")

# Open the PRN file and read in
f = open(args.prn_file, 'r')
prntext = f.readlines()
f.close()

# Loop and work out
for line in prntext:
    # Detect if the line is a "date =" entry or a VLBI scan
    lineinfo = line.split()
    
    if len(lineinfo) > 0 and lineinfo[0] == "date" and lineinfo[1] == "=": # it is a "date =" entry
        # Save the date in UTC:
        if args.verbose > 2: print("Found a new date line in the PRN file! Saving date to format time stamps.")
        datestr = lineinfo[2] # in str() format

    elif len(lineinfo) > 0 and lineinfo[0][0:2] == "no" and ( len(lineinfo) == 11 or len(lineinfo) == 12):
        # This is a row with a VLBI scan
        if args.verbose > 2: print("Found VLBI scan in the PRN file. Matching with CAL data from %s_PV.CAL.txt"%basenm) 
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

        if args.verbose >= 2:
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
        nearest_cal_ts_dt, nearest_cal_index, delta_t = nearest(source_ts_dt, timestamp_vlbi_dt)

        print("Nearest-in-time calibration scan found for VLBI scan %s on %s at %s is:"%(vlbi_scannumber, source_name, timestamp_vlbi_dt))
        closest_cal_info = "%s\tCALXML\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.3f\t%.1f\t%.1f\n"%(\
                             source_ts_dt[nearest_cal_index], source_calname[nearest_cal_index],\
                             source_az[nearest_cal_index], source_el[nearest_cal_index],\
                             source_tsys_b1r[nearest_cal_index], source_tsys_b1l[nearest_cal_index],\
                             source_tsys_b2r[nearest_cal_index], source_tsys_b2l[nearest_cal_index],\
                             source_tsys_b3r[nearest_cal_index], source_tsys_b3l[nearest_cal_index],\
                             source_tsys_b4r[nearest_cal_index], source_tsys_b4l[nearest_cal_index],\
                             source_tauz[nearest_cal_index], source_Tamb[nearest_cal_index], source_Tatm[nearest_cal_index])
        print(closest_cal_info)

        # Quick sanity check: double check that cal. scan source == VLBI scan science source:

        if source_calname[nearest_cal_index].lower() != source_name.lower(): # case insensitive check
            print("ERROR: Cal. scan source != from VLBI scan source! Halting program.")
            sys.exit(1)

        # Strategy 1: Check is delta_t between cal scan and VLBI scan < DELTA_TIME
        if delta_t.seconds <= DELTA_TIME:
            print("The cal. was taken within %.1f minutes of the VLBI scan start time. Accepting cal. metadata -> Writing out to Tsys table..."%(delta_t.seconds/60.))

        else:
            print("\n\n\n ERROR: Cal. data Delta_t to VLBI scan > 10 minutes. Need interpolation code here :-) .Halting.")
            sys.exit(1)


        # To interpolate we need to convert the dateime objects of timestamps to float numbers:
        source_ts_floats = [datetime_to_float(source_ts_dt[k]) for k in range(len(source_ts_dt))] 
        
        # Create interpolation function:

        # Note: this would be faster if we only interpolate 2-nearest values,...
        #       perhaps nearest could return a SORTED array of cals, ordered by delta_t
        #       but we need to select only the one before and after the VLBI scan.
        #       For the moment interpolate all the array, it is fast.
        x = np.array(source_ts_floats)
        y = np.array(source_tsys_b1r)
        f = interp1d(x, y)



        print("\n")

print "DONE."
