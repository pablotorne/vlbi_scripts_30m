#!/usr/bin/env python
"""
Read the results of the last calibration scan and return
them in the same format as getVlbiCalibration_Xmm.py to
mantain compatibility with possible codes reading from
that output format.

# The information has to be read from XML files
# created after each calibration.
# if .xml contains the string "-calibration-", it's a cal. If not, no.

v0: reads from /ncsServer/mrt/ncs/var/spool/odp/
v1: reads from /ncsServer/mrt/ncs/work/odp/continuum/old/
v2: fetches the xml file from /ncsServer/mrt/ncs/packages/coordinator2009-08-10v1.13/
    and makes a local copy to read from.
v3: added argparser options to have more flexibility. Now this code can be used by
    the VLBI monitor or the VLBI field system by passing the correct -lc parameter.

P. Torne, IRAM, v22.04.2018
"""

import numpy as np
import os, sys
from datetime import datetime, timedelta
import fnmatch
import time
import glob

import argparse

# Arguments control:
parser = argparse.ArgumentParser()

parser.add_argument("-v", "--verbose", action="store_true", help="Outputs detailed execution informaton")
parser.add_argument("-lc", "--localdatapath", help="Copies the xml files to this location and reads the information from there. \nRecommended mode. Make sure you have read/write permission in that folder!\nIf you prefer to not make local copies, do not pass this argument.\n")

args = parser.parse_args()

# Init values
CURRENTDAY = datetime.today().strftime("%Y%m%d") 
verbose = False
makelocalcopy = False
DATAPATH   = "/ncsServer/mrt/ncs/packages/coordinator2009-08-10v1.13/"
LOCALDATAPATH = "/mrt-lx3/vis/vlbi/vlbireduc/CALXMLs/"  # A standard, already created directory.
#LOCALDATAPATH = "/local/users/torne/vlbi_monitor_client/Python/CALXMLs/"
#LOCALDATAPATH = "/mrt-lx3/vis/vlbi/vlbireduc/CALXMLs/"  # version for vlbi Field System

# Format the program variables with the arguments passed:
if args.localdatapath != None: 
    LOCALDATAPATH = args.localdatapath
    makelocalcopy = True

if args.verbose == True: verbose = True

# ---

if verbose: print " "
if verbose: print "Setting the telescope datapath to %s"%DATAPATH
if verbose: print "Setting the localdatapath to %s"%LOCALDATAPATH

# FUNCTIONS
def get_lastmodified_file(glob_pattern, dir):
    '''
    Returns the file with latest modification time from
    a subset of files matching "glob_pattern" inside "dir"
    '''

    try:
        list_of_files = glob.glob(dir+"/"+glob_pattern)
        latest_file = max(list_of_files, key=os.path.getmtime)
    except:
        raise AttributeError('Cannot find latest XML file inside %s'%dir)

    return latest_file


def get_last_modified_folder(dir=DATAPATH):
    '''
    Returns the path of the last-modified subdir inside dir+"/scans",
    that should correspond to the newest scan available.
    '''
    
    # Get the full path of the last-modified subdir
    if verbose: print "Checking DATAPATH = %s ..."%dir
    if verbose: print "Is there a 'scans' subdir?"

    if os.path.isdir(dir+"/scans/"):
        if verbose: print "Yes. Going inside and getting most recent subdir (=scan)."

        scansdir = dir+"/scans/"
        all_subdirs = [scansdir+d for d in os.listdir(scansdir) if os.path.isdir(scansdir+d)]
        latest_subdir = max(all_subdirs, key=os.path.getmtime)

        if verbose: print "Found last modified subdir (=scan %s)! "%latest_subdir

    else:
        raise AttributeError("No 'scans' folder inside %s. Raising exception NotADirectoryError."%dir)

    return latest_subdir


def search_lastCal_XML(scandir):
    '''
    Iterate all scan subdirs backwards starting from scandir (the most recent one)
    until finding the calibration XML of continuum backend
    '''

    # Extract the scan number:    
    last_scan   = int( scandir.split("/")[-1] )
    datapath    = scandir.split("/"+str(last_scan))[-2]
    if verbose: print "last_scan = search_scan = %s"%last_scan


    # Iterate all scan subdirs backwards (latest first)
    # until finding the calibration XML of continuum backend
    search_scan = last_scan
    cal_xml = 'None'

    while cal_xml == 'None':

        if search_scan < 1: # There is no continuum calibration in this DAY, return
            #raise FileNotFoundError( errno.ENOENT, os.strerror(errno.ENOENT), scandir+"/"+search_scan)
            raise AttributeError('No CALs inside %s. FileNotFoundError.'%DATAPATH+"/")
            #sys.exit(1)

        search_folder = datapath+"/"+str(search_scan)
        if verbose: print "Searching continuum cal XML in %s"%search_folder

        try:
            for file in os.listdir(search_folder):
                # Search first NBC calibrations...
                if fnmatch.fnmatch(file, '*-calibration-NBC-*'):
                    if verbose: print "Found an NBC continuum CAL !!!"
                    cal_xml = file
                    backend = 'NBC'
                    break
                # If not NBC cals, search for BBC cals ...
                elif fnmatch.fnmatch(file, '*-calibration-BBC-*'):
                    if verbose: print "Found a BBC continuum CAL !!!"
                    cal_xml = file
                    backend = 'BBC'
                    break
                else:
                    # backend is not identified in filename (can be read from inside xml file if needed)
                    continue
            
            if cal_xml == 'None': # Reduce scan number and search again
                search_scan = search_scan-1

        except Exception as e:
            print "Exception in search_lastCal_XML(%s)."%search_folder
            print e
            # Ignore and go on with previous scan:
            #search_scan = search_scan-1
            sys.exit(1)

    return [datapath+"/"+str(search_scan)+"/"+cal_xml, backend]


# MAIN PROGRAM BODY

if verbose: print "Starting read_lastCal_info.py ..."

# Get the last -calibration-*$DATE*.xml file from 
success_get_cal = False
days_without_cals = 0
DATE = CURRENTDAY

while not success_get_cal:

    if days_without_cals >= 5:
        print "No recent calibration results found for the project. Wait for a new calibration scan to get results."
        #print "5+ days without finding a result of cal in xml. Aborting read_lastCal_info.py ..."
        sys.exit(1)

    try: # Find a cal xml file for $DATE.

        try: # Is there an XML file in coordinator2009/ dir? Yes? -> Get the last one and make a local copy

            if verbose: print "Trying to fetch a cal xml file from %s ..."%DATAPATH
            xml_cal = get_lastmodified_file("*-calibration*%s*.xml"%DATE, DATAPATH)
            if verbose: print "OK. Found %s ..."%xml_cal

            #Yes? -> Get the last one and make a local copy. Check if the file was already copied and if so, skip copy.
            xml_cal_filenm=xml_cal.split("/")[-1]
            if verbose: print "Extracting filenm of xml_cal: %s"%xml_cal_filenm

            if makelocalcopy:

                if not os.path.exists(LOCALDATAPATH+xml_cal_filenm):
                    if verbose: 
                        print "Creating a local copy of the cal xml: cp %s %s"%(xml_cal, LOCALDATAPATH+xml_cal_filenm)
                    status = subprocess.call('cp %s s%'%(xml_cal, LOCALDATAPATH+xml_cal_filenm), shell=True)
                    if status != 0: print "ERROR copying xml file to /local/users/torne/vlbi_monitor_client/Python/CALXMLs/"
                    xml_cal = get_lastmodified_file("*-calibration*%s*.xml"%DATE, LOCALDATAPATH)
                else:
                    if verbose: print "%s already in %s. Skipping copy."%(xml_cal_filenm, LOCALDATAPATH)
                    xml_cal = get_lastmodified_file("*-calibration*%s*.xml"%DATE, LOCALDATAPATH)

            else: # If option makelocalcopy is False, use the file directly from DATAPATH:
                if verbose: print "makelocalcopy = False. Using files in %s"%DATAPATH
                xml_cal = get_lastmodified_file("*-calibration*%s*.xml"%DATE, DATAPATH)     
             
            success_get_cal = True

        except:  # No cal xml file in coordinator2009?

            try: #search for the file in the local folder where we keep the copies (CALXML).

                if makelocalcopy:
   
                    if verbose:
                        print "No XML found in %s"%DATAPATH
                        print "Checking %s"%LOCALDATAPATH

                    xml_cal = get_lastmodified_file("*-calibration*%s*.xml"%DATE, LOCALDATAPATH)

                    if verbose: print "Success. Found %s"%xml_cal
                    
                    success_get_cal = True

                else: # makelocalcopy = False
                    if verbose: print "No calibration results found in %s for %s."%(DATAPATH, DATE)
                    sys.exit(1)

            except:
                if verbose: print "No cal results found for %s. Trying on previous day.\n"%DATE
                sys.exit(1)

    except: # If all the previous try fails, it is because we changed the DATE and there is no CAL this day yet.
            # Reduce the variable DATE by one day and search again for files to fetch the latest cal.
        #raise AttributeError("No -cal*.XML files found on %s"%DATE)
        days_without_cals += 1
        DATE       = (datetime.today()-timedelta(days=days_without_cals)).strftime("%Y%m%d")
        #time.sleep(2)


# 3) Read contents and extract desired info:

if verbose: print "Reading the XML file ..."
f = open(xml_cal)
info = f.readlines()
f.close()

if verbose: print "Extracting calibration results ..."
#if verbose: print "Backend = %s"%backend

nIFs = 8

for ii in range(len(info)):
    if fnmatch.fnmatch(info[ii], '*<RESOURCE name="calibration">\n'):
        if verbose: print "Found the CAL TABLE!"
        start = ii
        # ii now points to the start of the table with the cal results

        # For BBC, the desired data (and line numbers) are:
        # rx name:   ii+30
        # Tsys:      ii+46
        # tau:       ii+50
        # pwv mm:    ii+44
        # trx:       ii+45
        # sky_freq:  ii+32

        # then add 27 lines for IF 2, 27*2 for IF3, and 27*3 for IF4.
        # Be aware that NBC and BBC can be different. NBC 4 channels, BBC 8 channels

        # To-Do: Differenciate BBC and NBC
        #NBC calibration xml files has 340 lines
        #BBC calibration xml files have 540 lines
        # I dont know how many NBC+BBC xml files have.

        # Get RX NAMES           
        rx_arr = []
        for jj in range(30, 30+27*nIFs, 27):

            try:

                rx_name = info[start+jj].split("<TD>")[1].split("</TD>")[-2] 
                rx_arr.append( rx_name )

            except:

                if verbose: print "Exception raised when retrieving rx names"
                rx_name = "Not_tuned"
                rx_arr.append( rx_name )
                #break

        if verbose: print "rx_array = %s"%rx_arr

        # Get Rx SkyFreq.
        rxFreq_arr = []
        for jj in range(32, 32+27*nIFs, 27):

            try:

                rx_skyf = info[start+jj].split("<TD>")[1].split("</TD>")[-2]
                rxFreq_arr.append( rx_skyf )

            except:

               if verbose: print "Exception raised when retrieving rxFreqs"
               rx_skyf= -1
               rxFreq_arr.append( rx_skyf )
               #break

        if verbose: print "rxFreq_array = %s"%rx_arr

    
        # Get Tsys
        tsys_arr = []
        for jj in range(46, 46+27*nIFs, 27):

            try:

                tsys = info[start+jj].split("<TD>")[1].split("</TD>")[-2]
                tsys_arr.append( tsys )

            except:

               if verbose: print "Exception raised when retrieving tsys values"
               tsys = -1
               tsys_arr.append( tsys )
               #break

        if verbose: print "tsys_array = %s"%tsys_arr

        # Get tau     
        tau_arr = []
        for jj in range(50, 50+27*nIFs, 27):

            try:

                tau = info[start+jj].split("<TD>")[1].split("</TD>")[-2]
                tau_arr.append( tau )

            except:

               if verbose: print "Exception raised when retrieving tau values"
               tau = -1
               tau_arr.append( tau )
               #break

        if verbose: print "tau_array = %s"%tau_arr

        # Get pwv mm     
        pwv_arr = []
        for jj in range(44, 44+27*nIFs, 27):

            try:

                pwv = info[start+jj].split("<TD>")[1].split("</TD>")[-2]
                pwv_arr.append( pwv )

            except:

               if verbose: print "Exception raised when retrieving pwv values"
               pwv = -1
               pwv_arr.append( pwv )
               #break

        if verbose: print "pwv_array = %s"%pwv_arr

        # Get Trx
        trx_arr = []
        for jj in range(45, 45+27*nIFs, 27):

            try:

                trx = info[start+jj].split("<TD>")[1].split("</TD>")[-2]
                trx_arr.append( trx )

            except:

               if verbose: print "Exception raised when retrieving trx values"
               trx = -1
               trx_arr.append( trx )
               #break

        if verbose: print "trx_array = %s"%trx_arr


        break # Don't continue searching the XML file if the CALinfo table is found

# Read the time stamp of the data
for ii in range(len(info)):
    if fnmatch.fnmatch(info[ii], '*<PARAM name="timeStamp" value=*'):
        if verbose: print "Found the time stamp!"

        timestamp = info[ii].split('value="')[1].split('" datatype="')[-2]
        # Keep format with getVlbiCals1mm.py:
        timestamp = timestamp.replace("T", " ", 1)[:-4]
        if verbose: print "Timestamp = %s"%timestamp

# Read the elevation of the source:      <PARAM name="elevation" value="
for ii in range(len(info)):
    if fnmatch.fnmatch(info[ii], '*<PARAM name="elevation" value="*'):
        if verbose: print "Found the elevation!"

        elev = info[ii].split('value="')[1].split('" unit="')[-2]
        if verbose: print "elev = %s"%elev

# Read the source name:    <PARAM name="sourceName" value="
for ii in range(len(info)):
    if fnmatch.fnmatch(info[ii], '*<PARAM name="sourceName" value=*'):
        if verbose: print "Found the source!"

        source = info[ii].split('value="')[1].split('" datatype="')[-2]
        if verbose: print "Source = %s"%source

# Output info in same format as getVlbiCalibration1mm.py

# Check that all arrays are same length:
length = len(rx_arr)
if any(len(lst) != length for lst in [tsys_arr, tau_arr, pwv_arr, trx_arr]):
    print "Some data missing. Not broadcasting data."
    sys.exit(1)
else:
    print "CAL info read from: %s"%xml_cal
    for ii in range(len(rx_arr)):
        if rx_arr[ii] != "Not_tuned":
            print "rx: %s ; tsys: %.2f ; tau: %.2f ; pwv mm: %.1f ; trx: %.1f ; time: %s ; rxFreq: %.4f ; elev: %.1f ; source: %s"%(\
             rx_arr[ii], float(tsys_arr[ii]), float(tau_arr[ii]), float(pwv_arr[ii]), float(trx_arr[ii]), timestamp, float(rxFreq_arr[ii]), float(elev), source)
