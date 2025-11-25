#!/usr/bin/env python3
'''
Read in a .vex file, identify the scans for IRAM30m,
and show a table in the terminal with the information
of each scan. The last two columns will show in real
time how much time do we have left until the next
action, and propose a pointing or other commands for
the observer.

P. Torne + ChapGPT, IRAM, v2024.12.30
Bug fixes and improvements in April 2025 (full working version)
NOTE: We use f-strings and thus require python3.6 or later!
(So, code will not work in mrt-lx3) --
Change syntax to "".format() for backwards compatibility!
'''

import argparse
import time
import datetime

# File path (replace with the actual path if running locally)
file_path = './c242a.vex'


def parse_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract scan blocks
    scan_blocks = []
    inside_scan = False
    current_block = []
    contains_station_pv = False  # Tracks if the block contains "station=Pv"

    for line in lines:
        line = line.strip()
        if line.startswith("scan"):
            inside_scan = True
            current_block = [line]
            contains_station_pv = False  # Reset flag for each new block
        elif line.startswith("endscan") and inside_scan:
            current_block.append(line)
            if contains_station_pv:  # Only add the block if "station=Pv" is present
                scan_blocks.append(current_block)
            inside_scan = False
        elif inside_scan:
            if not line.startswith("station="):
                current_block.append(line)
            # Check for "station=Pv" and update the flag
            if line.startswith("station=Pv"):
                contains_station_pv = True
                current_block.append(line)

    return scan_blocks


def extract_scan_data(scan_blocks):
    scans = []
    for block in scan_blocks:
        scan_data = {"scan_number": None, "source_name": None, "start_time": None, "duration": None, "end_time": None}
        for line in block:
            print(f"\n{line}")
            if line.startswith("scan"):
                scan_data["scan_number"] = line.split()[1]
                print(f'Extracting scan_number = {scan_data["scan_number"]}')
            elif "start" in line.lower() and "source" in line.lower():

                scan_data["source_name"] = line.split("=")[-1].strip()
                print(f'Extracting source_name = {scan_data["source_name"]}')

                scan_data["start_time"] = line.split("=")[1].strip().split(";")[0]
                print(f'Extracting start = {scan_data["start_time"]}')

                # Digest the string with the date and time to create a datetime object:
                # Extract components from the string
                year = int(scan_data["start_time"][:4])
                day_of_year = int(scan_data["start_time"][5:8])
                hour = int(scan_data["start_time"][9:11])
                minute = int(scan_data["start_time"][12:14])
                second = int(scan_data["start_time"][15:17])

                # Convert the day of the year to a date
                date = datetime.datetime(year, 1, 1) + datetime.timedelta(days=day_of_year - 1)

                # Create a datetime object combining the date and time
                datetime_obj_start = datetime.datetime.combine(date, datetime.time(hour, minute, second))

                scan_data["start_time"] = datetime_obj_start

                # Display the results
                print(f"Date: {datetime_obj_start.date()}")
                print(f"Time: {datetime_obj_start.time()}")
                print(f"Datetime: {datetime_obj_start}")

            elif "station=pv" in line.lower():
                scan_data["duration"] = float( line.split("sec:")[1].strip() )
                print(line.split("sec:"))
                print(f'Extracting duration = {scan_data["duration"]}')

        # Calculate end time if start time and duration are available
        if scan_data["start_time"] and scan_data["duration"]:
            # Create a timedelta object for the duration
            delta = datetime.timedelta(seconds=scan_data["duration"])
            scan_data["end_time"] = (datetime_obj_start + delta)

        scans.append(scan_data)

    return scans

def action_selector(next_start, current_end, next_end, scan_duration):
    """
    Returns a message depending on how much time is left for the next scan.
    """

    # Use a real-time time_diff to choose action:
    if next_start < datetime.datetime.utcnow():
        # It is a past scan (?). But this is needed to control when to show "PAST SCAN"
        time_diff = next_start - datetime.datetime.utcnow()
    elif next_start >= datetime.datetime.utcnow() and current_end < datetime.datetime.utcnow():
        # We are within a scan gap: update in real-time the action by updating time_diff:
        time_diff = next_start - datetime.datetime.utcnow()
    else:
        # We are in a future scan. Show the action based on the full gap time, until we reach the gap
        time_diff = next_start - current_end

    #print("Time Difference:", time_diff)

    if time_diff == datetime.timedelta(minutes=99):
        return "End of schedule"
    if time_diff >= datetime.timedelta(minutes=20):
        return "Planet"
    elif time_diff < datetime.timedelta(minutes=20) and time_diff >= datetime.timedelta(minutes=9):
        return "p8 || p4 + ff"
    elif time_diff < datetime.timedelta(minutes=9) and time_diff >= datetime.timedelta(minutes=7):
        return "p8 || p4 + ff"
    elif time_diff < datetime.timedelta(minutes=7) and time_diff >= datetime.timedelta(minutes=6):
        return "p8f || p4f + ff"
    elif time_diff < datetime.timedelta(minutes=6) and time_diff >= datetime.timedelta(minutes=4.66):
        return "p8f || p4"
    elif time_diff < datetime.timedelta(minutes=4.66) and time_diff >= datetime.timedelta(minutes=3.5):
        return "p4"
    elif time_diff < datetime.timedelta(minutes=3.5) and time_diff >= datetime.timedelta(minutes=3):
        return "p4f"
    elif time_diff < datetime.timedelta(minutes=3) and time_diff >= datetime.timedelta(minutes=2.6):
        return "p4ff"
    elif time_diff < datetime.timedelta(minutes=2.6) and time_diff >= datetime.timedelta(seconds=40):
        return "cc"
    elif time_diff < datetime.timedelta(seconds=40) and time_diff >= datetime.timedelta(seconds=1):
        return "\033[93mDo nothing (but check correct source)\033[0m"
    elif time_diff >= datetime.timedelta(seconds= -1 * scan_duration):
        #print(next_end)
        #print((datetime.datetime.utcnow() - next_end).total_seconds())
        return f"\033[1;31mWAIT! RECORDING SCAN!\033[0m [{(datetime.datetime.utcnow() - next_end).total_seconds():.0f} sec.]"
    elif time_diff < datetime.timedelta(seconds= -1 * scan_duration):
        return "PAST SCAN"
    else:
        return "???"

def format_time_to_next(delta):
    """Format timedelta to show only minutes and seconds."""
    total_seconds = int(delta.total_seconds())
    minutes, seconds = divmod(total_seconds, 60)
    return f"{minutes:02}:{seconds:02}"


def display_table(scans, vexfile):
    sched = vexfile.split(".vex")[0]
    while True:
        print("\033c", end="")  # Clear terminal
        print(f"\033[1mIRAM 30m VLBI Observer Co-Pilot - Schedule: {sched}\033[0m")
        print("(Dev. by Heino Falcke & Pablo Torne -- Warning: Beta version!)")
        print(f"Current UTC Time: {datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'Scan Number':<15}{'Source':<15}{'Start Time (UTC)':<20}{' End Time':<15}{'Duration':<15}{'Gap[mm:ss]':<15}{'Recommended Action':<15}")
        print("-" * 115)

        for i, scan in enumerate(scans):
            time_to_next = "N/A"
            if i < len(scans) - 1: # If not the last scan, use info of next scan to calculate times
                next_start    = scans[i + 1]["start_time"]
                next_end      = scans[i + 1]["end_time"]
                recording_time = scans[i + 1]["duration"]
            else: # If we are in the last scan, we create fake next_start and next_end values so the code can work well
                #next_start    = scans[i]["start_time"]+datetime.timedelta(seconds=scans[i]["duration"])
                next_start    = scans[i]["end_time"]+datetime.timedelta(seconds=5940) # 5940 sec. are 99min, to mark end of schedule
                next_end      = scans[i]["end_time"]+datetime.timedelta(seconds=6000) # Dummy number, will not be used but parameter need to exist
                recording_time = 5940 # 99 min, placeholder value only

            current_start = scans[i]["start_time"]
            current_end   = scans[i]["end_time"]
            now = datetime.datetime.utcnow()

            # Dynamically show the GAP time, but update it in real time when you reach a given gap:
            if current_end < now:
                time_to_next = next_start - datetime.datetime.utcnow()
            else:
                time_to_next = next_start - current_end

            #print(time_to_next)
            #print(format_time_to_next(time_to_next))

            time_to_next_str = format_time_to_next(time_to_next) if next_start >= datetime.datetime.utcnow() else "Started/Done"
            # Get the date and time in string format for printing:
            start_time = scan['start_time'].strftime("%Y-%m-%d %H:%M:%S")
            end_time   = scan['end_time'].strftime(" %H:%M:%S")
            ##action = action_selector( next_start - datetime.datetime.utcnow(), scan['duration'] )
            #action = action_selector( next_start, current_end, next_start - current_end, scan['duration'] )
            ###action = action_selector( next_start, current_end, next_end, scan['duration'] ) # Bug, message "RECORDING WRONG timing"
            action = action_selector( next_start, current_end, next_end, recording_time )
            # Apply bold formatting if the time to next is greater than 10 minutes
            if time_to_next.total_seconds() >= 600: # Make bold text
                ##time_to_next_str = f"\033[1m{time_to_next_str}\033[0m"
                action           = f"\033[1m{action}\033[0m"

            # Print row; Duration in minutes, rounded to 1 decimal:
            print(f"{scan['scan_number']:<15}{scan['source_name']:<15}{start_time:<20}{end_time:<15}{round(scan['duration']/60.,1):<15}{time_to_next_str:<15}{action:<15}")
            #print(f"{scan['scan_number']:<15}{scan['source_name']:<15}{start_time:<20}{end_time:<15}{round(scan['duration']/60.,1):<15}")
            #print(f"{' '*80}{time_to_next_str:<15}{action:<15}")

        time.sleep(1)



def main():
    # Argument parser
    parser = argparse.ArgumentParser(description="Process .vex and show ASCII table with action suggestions.")
    parser.add_argument("file_path", type=str, help="Path to the .vex file")

    # Parse arguments
    args = parser.parse_args()

    # Use the provided file path
    file_path = args.file_path

    # Parse the file
    scan_blocks = parse_file(file_path)

    # Extract the scan data
    scans = extract_scan_data(scan_blocks)

    # Display the table
    display_table(scans, file_path)

if __name__ == "__main__":
    main()



