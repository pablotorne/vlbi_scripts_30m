# vlbi_scripts_30m
Scripts aiding VLBI observations at the 30-m radio telescope

read_lastcal_info.py: 
- Searches and reads the last results from a calibration scan and prints them formatted properly on the terminal. 
- Used by the VLBI web monitor and the VLBI Field System to monitor the cal. results.
- The results on this script goes into the Field System logs, that will later be used by other scripts to format
  metadata in different formats (e.g., ANTAB, or EHT_metadata_tables).

cals30m_ANTAB.py: 
- Reads files with calibration and weather information form the IRAM 30m (usually a Field System log) and writes the output as ANTAB tables.

TsysTable_EHT.py
- Read files with calibration information from the IRAM 30m (usually a Field System log) and writes the output as EHT Tsys* metadata tables.

check_duplicate_sources.py
- Reads the output of a paKo .sou catalog and detects duplicate sources

---

For master LOCAL code development use mrt-lx3 (vlbi)


