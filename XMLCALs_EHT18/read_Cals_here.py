#!/usr/bin/env python
"""
Adaptation of ~/src/read_latCal_info.py to digest the xml
cal files of the first 2 days of the EHT 2018 run
for Helge/Thomas.

P. Torne, IRAM, v22.04.2018
"""

import numpy as np
import os, sys
from datetime import datetime, timedelta
import fnmatch
import time
import glob

verbose    = False

if verbose: print " "

def extract_cal_info(xmlfile):
    '''
    Read contents and extract desired info
    from a cal result xml file
    '''

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
        calresults= []
        calresults.append("CAL info read from: %s"%xml_cal)
        for ii in range(len(rx_arr)):
            calresults.append( "rx: %s ; tsys: %.2f ; tau: %.2f ; pwv mm: %.1f ; trx: %.1f ; time: %s ; rxFreq: %.4f ; elev: %.1f ; source: %s"%(\
             rx_arr[ii], float(tsys_arr[ii]), float(tau_arr[ii]), float(pwv_arr[ii]), float(trx_arr[ii]), timestamp, float(rxFreq_arr[ii]), float(elev), source) )

    return calresults


# PROGRAM BODY

if verbose: print "Starting read_Cals_here.py ..."

# Open a file to write results:
out = open("infoCals_EHT2018_firstdays.txt", "w")

# Loop over the xml calibration files, sorted by modification date:
search_dir = "./"
files = filter(os.path.isfile, glob.glob(search_dir + "*-calibration-*.xml"))
files.sort(key=lambda x: os.path.getmtime(x))
#files = files[-1::] # Reverse the list to get the oldest files first

for xml_cal in files:
    calresults = extract_cal_info(xml_cal)
    for line in calresults:
        print line
        out.write(line+"\n")

out.close()
