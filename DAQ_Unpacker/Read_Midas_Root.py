import argparse
import os
import struct
import sys

import midas.file_reader
import numpy as np
import ROOT

# Make sure AcquisitionModes is in your PYTHONPATH
from AcquisitionModes import (
    DTQ_COUNT,
    DTQ_SERVICE,
    DTQ_SPECT,
    DTQ_TIMING,
    DTQ_TSPECT,
    Counting,
    Service,
    Spectroscopy,
    Timing,
)

# -------------------------------------------------------------------------------
# CONFIGURATION PATHS
# -------------------------------------------------------------------------------
MIDAS_DATA_DIR = "/home/lorenzo/muEDM_Project/Data/MidasData/"
ROOT_DATA_DIR = "/home/lorenzo/muEDM_Project/Data/RootData/"

# -------------------------------------------------------------------------------
# ARGUMENT PARSING
# -------------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Midas to ROOT Converter (Multi-Board Support)"
)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("run_number", nargs="?", type=int, help="Run number (e.g., 291)")
group.add_argument("-P", "--path", type=str, help="Full path to input file")

args = parser.parse_args()

mfile_string = ""
if args.path:
    mfile_string = args.path
else:
    run_str = f"run{args.run_number:05d}"
    mfile_string = os.path.join(MIDAS_DATA_DIR, f"{run_str}.mid.lz4")

if not os.path.exists(mfile_string):
    print(f"ERROR: Input file does not exist: {mfile_string}")
    sys.exit(1)

print(f"Input file: {mfile_string}")

# -------------------------------------------------------------------------------
# OPEN FILE & PARSE ODB FOR ALL BOARDS
# -------------------------------------------------------------------------------
mfile = midas.file_reader.MidasFile(mfile_string)

# Dictionary to store ODB settings for each board ID
# Format: { 0: {'delay': 10, ...}, 1: {'delay': 20, ...} }
odb_settings_map = {}

print("--- Reading ODB (Configuration) ---")
try:
    odb = mfile.get_bor_odb_dump()
    if odb is not None:
        data = odb.data
        # Navigate to the FersDAQ settings
        # Path: Equipment -> FersDAQ -> Settings -> FersDAQ
        fers_root = (
            data.get("Equipment", {})
            .get("FersDAQ", {})
            .get("Settings", {})
            .get("FersDAQ", {})
        )

        if fers_root:
            # Iterate over all keys (FLab0, FLab1, FLab2, etc.)
            for key, val in fers_root.items():
                if key.startswith("FLab"):
                    try:
                        # Extract ID from string "FLabX"
                        board_id = int(key.replace("FLab", ""))

                        # Store parameters
                        odb_settings_map[board_id] = {
                            "delay": val.get("TimeReferenceDelay_ns", 0),
                            "window": val.get("TimeReferenceWindow_ns", 0),
                            "thr": val.get("TCoarseThr", 0),
                        }
                        print(
                            f"  [ODB] Found config for Board {board_id}: "
                            f"Delay={odb_settings_map[board_id]['delay']}ns, "
                            f"Win={odb_settings_map[board_id]['window']}ns"
                        )
                    except ValueError:
                        pass  # Key looked like FLab but wasn't FLab+Int
        else:
            print("  WARNING: Path 'Equipment/FersDAQ/Settings/FersDAQ' not found.")
    else:
        print("  WARNING: No ODB dump found.")

except Exception as e:
    print(f"  ERROR READING ODB: {e}")

# Reset file pointer
mfile.jump_to_start()

# -------------------------------------------------------------------------------
# PREPARE ROOT OUTPUT
# -------------------------------------------------------------------------------
os.makedirs(ROOT_DATA_DIR, exist_ok=True)
base_name = os.path.basename(mfile_string)
root_name = base_name.replace(".mid.lz4", ".root").replace(".mid", ".root")
if not root_name.endswith(".root"):
    root_name += ".root"

rootfile_string = os.path.join(ROOT_DATA_DIR, root_name)
print(f"Output file: {rootfile_string}")

f = ROOT.TFile(rootfile_string, "RECREATE")

# --- MAIN EVENT TREE ---
tree = ROOT.TTree("raw", "Main event tree")
event_id_val = np.zeros(1, dtype=np.int32)
tree.Branch("EventID", event_id_val, "EventID/I")

# --- CONFIGURATION TREE ---
# This tree will have 1 entry per Board present in the run
config_tree = ROOT.TTree("cfg", "Run Configuration and Board Status")

# Variables for the Config Tree
conf_board_id = np.zeros(1, dtype=np.int32)  # Which board is this row for?
conf_qualifier = np.zeros(1, dtype=np.int32)
conf_tstamp = np.zeros(1, dtype=np.double)  # Service timestamp
conf_time = np.zeros(1, dtype=np.uint64)  # Update time
# ODB variables
conf_delay = np.zeros(1, dtype=np.double)
conf_window = np.zeros(1, dtype=np.double)
conf_thr = np.zeros(1, dtype=np.double)

# Branch definitions
config_tree.Branch("BoardID", conf_board_id, "BoardID/I")
config_tree.Branch("DataQualifier", conf_qualifier, "DataQualifier/I")
config_tree.Branch("serv_tstamp", conf_tstamp, "serv_tstamp/D")
config_tree.Branch("serv_time", conf_time, "serv_time/l")
# ODB Branches
config_tree.Branch("TimeRefDelay", conf_delay, "TimeRefDelay/D")
config_tree.Branch("TimeRefWindow", conf_window, "TimeRefWindow/D")
config_tree.Branch("Threshold", conf_thr, "Threshold/D")

# -------------------------------------------------------------------------------
# PASS 1: SCAN FOR BANKS & FILL CONFIG TREE
# -------------------------------------------------------------------------------
# We need to find which boards exist (FDxx/FSxx) and grab their Service Info
# Then we match them with the ODB settings we loaded earlier.

found_physics = False
found_service = False
boards_found = set()  # Store IDs of boards we've processed

# Temporary storage for Service data before filling tree
# Format: { board_id: ServiceObject }
service_data_map = {}
# Format: { board_id: DataQualifier }
qualifier_map = {}

print("Scanning file for board configuration...")

for event in mfile:
    if event.header.is_midas_internal_event():
        continue

    # Check if we have gathered enough info to stop scanning
    # (Arbitrary logic: stop if we saw both Physics and Service events)
    if found_physics and found_service:
        break

    # -- 1. Analyze Physics Events (ID 1) to define Event Tree branches --
    if event.header.event_id == 1:
        found_physics = True
        for bank_name, bank in event.banks.items():
            if bank_name.startswith("FD") and bank_name[2:].isdigit():
                board_id = int(bank_name[2:])

                # Store Qualifier for Config Tree
                raw = bytes(bank.data)
                dtq = struct.unpack("<I", raw[:4])[0]
                qualifier_map[board_id] = dtq

                # Define Dynamic Branches (Only need to do this once per bank name)
                if bank_name not in [b.GetName() for b in tree.GetListOfBranches()]:
                    print(f"  -> Defining branches for {bank_name} (ID {board_id})")

                    if dtq == DTQ_TIMING:
                        # Create arrays in a dictionary to keep references alive
                        # (Using a global dict 'branches' to store arrays)
                        if "branches" not in locals():
                            branches = {}

                        branches[f"{bank_name}_fine_tstamp"] = np.zeros(
                            1, dtype=np.float64
                        )
                        branches[f"{bank_name}_nhits"] = np.zeros(1, dtype=np.uint32)
                        branches[f"{bank_name}_channel"] = np.zeros(
                            2048, dtype=np.uint8
                        )
                        branches[f"{bank_name}_ToA"] = np.zeros(2048, dtype=np.uint32)
                        branches[f"{bank_name}_ToT"] = np.zeros(2048, dtype=np.uint16)

                        tree.Branch(
                            f"{bank_name}_fine_tstamp",
                            branches[f"{bank_name}_fine_tstamp"],
                            f"{bank_name}_fine_tstamp/D",
                        )
                        tree.Branch(
                            f"{bank_name}_nhits",
                            branches[f"{bank_name}_nhits"],
                            f"{bank_name}_nhits/i",
                        )
                        tree.Branch(
                            f"{bank_name}_channel",
                            branches[f"{bank_name}_channel"],
                            f"{bank_name}_channel[{bank_name}_nhits]/b",
                        )
                        tree.Branch(
                            f"{bank_name}_ToA",
                            branches[f"{bank_name}_ToA"],
                            f"{bank_name}_ToA[{bank_name}_nhits]/i",
                        )
                        tree.Branch(
                            f"{bank_name}_ToT",
                            branches[f"{bank_name}_ToT"],
                            f"{bank_name}_ToT[{bank_name}_nhits]/s",
                        )

                    # Add blocks for SPECT/COUNTING here if needed

    # -- 2. Analyze Service Events (ID 2) to get timestamps --
    if event.header.event_id == 2:
        found_service = True
        for bank_name, bank in event.banks.items():
            if bank_name.startswith("FS") and bank_name[2:].isdigit():
                board_id = int(bank_name[2:])
                # Decode Service
                raw = bytes(bank.data)
                ev = Service(raw)
                service_data_map[board_id] = ev

# --- FILL CONFIG TREE ---
# Now we combine: Service Data + ODB Settings -> RunConfig Tree
print("Filling RunConfig Tree...")

# Get list of all unique Board IDs found in Physics OR Service
all_ids = set(qualifier_map.keys()) | set(service_data_map.keys())

for bid in sorted(all_ids):
    # 1. Board ID
    conf_board_id[0] = bid

    # 2. Data Qualifier (from Physics event)
    conf_qualifier[0] = qualifier_map.get(bid, 0)

    # 3. Service Info (from Service event)
    if bid in service_data_map:
        conf_tstamp[0] = service_data_map[bid].tstamp_us
        conf_time[0] = service_data_map[bid].update_time
    else:
        conf_tstamp[0] = 0
        conf_time[0] = 0

    # 4. ODB Settings (from XML dump)
    if bid in odb_settings_map:
        conf_delay[0] = odb_settings_map[bid]["delay"]
        conf_window[0] = odb_settings_map[bid]["window"]
        conf_thr[0] = odb_settings_map[bid]["thr"]
    else:
        # If ODB doesn't have this board (or ODB failed), set to 0
        conf_delay[0] = 0
        conf_window[0] = 0
        conf_thr[0] = 0

    print(f"  -> Board {bid}: Delay={conf_delay[0]}, Window={conf_window[0]}")
    config_tree.Fill()

# -------------------------------------------------------------------------------
# PASS 2: MAIN EVENT LOOP
# -------------------------------------------------------------------------------
mfile.jump_to_start()
print("Starting main event conversion...")

for j, event in enumerate(mfile):
    if event.header.is_midas_internal_event():
        continue
    if event.header.event_id != 1:
        continue

    # Quick check: does this event have FD banks?
    # (Optimization: avoid processing empty events)
    has_fd = False
    for k in event.banks.keys():
        if k.startswith("FD"):
            has_fd = True
            break
    if not has_fd:
        continue

    if j % 2000 == 0:
        sys.stdout.write(f"Processing event #{j}\r")
        sys.stdout.flush()

    # Clear previous event data
    for val in branches.values():
        if isinstance(val, np.ndarray):
            val[:] = 0

    # Process all banks in this event
    for bank_name, bank in event.banks.items():
        if not (bank_name.startswith("FD") and bank_name[2:].isdigit()):
            continue

        # Determine data type from first 4 bytes
        raw = bytes(bank.data)
        dtq = struct.unpack("<I", raw[:4])[0]

        if dtq == DTQ_TIMING:
            timing_ev = Timing(raw)

            # Fill Header info for this bank
            branches[f"{bank_name}_fine_tstamp"][0] = timing_ev.fine_tstamp
            branches[f"{bank_name}_nhits"][0] = timing_ev.nhits

            # Fill Hits
            limit = min(len(timing_ev.hits), 2048)
            for i in range(limit):
                hit = timing_ev.hits[i]
                branches[f"{bank_name}_channel"][i] = hit["channel"]
                branches[f"{bank_name}_ToA"][i] = hit["tstamp"]
                branches[f"{bank_name}_ToT"][i] = hit["tot"]

        # (Add DTQ_SPECT / DTQ_COUNT logic here if needed)

    tree.Fill()
    event_id_val[0] += 1

print("\nSaving ROOT file...")
config_tree.Write()  # Explicitly write the config tree
tree.Write()  # Write the event tree
f.Close()
print("Conversion Complete.")
