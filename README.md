# Backing_array_analysis_code

📊 Event Building & External Sorting for Timing Data (CSV / DAT → ROOT)

This repository provides high-performance C++ tools to process large-scale timing data (CSV or binary DAT) from CAEN/Janus systems and convert them into structured ROOT files with:

Absolute timestamp reconstruction
External merge sort (memory-efficient for large datasets)
Event building using time-window clustering
Ring mapping (Board_Id, CH_Id → Ring)
ROOT TTree output for physics analysis
🚀 Features
✅ Supports CSV and binary DAT input
✅ Handles very large files via external sort (chunk + k-way merge)

✅ Computes absolute time:

t_abs_us = TStamp_us + 0.0005 × ToA_LSB
✅ Event grouping based on configurable time window
✅ Ring mapping for detector geometry
✅ ROOT output for downstream analysis
⚠️ Known Limitation: DAT Input & Ring Identification
❌ Current Status
Ring identification does NOT work correctly for DAT input
Ring mapping works correctly for CSV input
🔍 Reason

The DAT format currently parsed in this repository does not include (or does not decode) Board_Id.

Current implementation:

r.Board_Id = default_board_id;  // hardcoded (currently = 0)

👉 As a result:

All hits are assigned to a single board
Ring mapping (Board_Id, CH_Id → Ring) becomes incorrect
⚠️ Impact

For DAT-based workflows:

❌ Ring assignments are unreliable
❌ Multi-board detector geometry is lost
❌ Ring-wise histograms and coincidence logic may be incorrect
✅ Recommended Usage
Use Case	Recommended Input
Full analysis (ring-based)	✅ CSV
Timing / debugging	✅ DAT
Event grouping only	✅ DAT
🔧 Workarounds
Prefer CSV input for any ring-based analysis
Use DAT only for:
timing validation
event clustering
debugging raw hits
🚧 Planned Improvements
Extract Board_Id from DAT (if present in header/structure)
Support multi-board reconstruction
Improve DAT parsing performance
📦 Requirements
C++17 compiler
ROOT (tested with ROOT 6.x)

Load ROOT environment:

source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
🔧 Compilation
g++ -O2 -std=c++17 csv_or_dat_to_root_events_external_sort_with_ring.cpp \
    $(root-config --cflags --libs) \
    -o csv_or_dat_to_root_events_external_sort_with_ring
▶️ Usage
./csv_or_dat_to_root_events_external_sort_with_ring \
    <input.(csv|dat)> \
    <output.root> \
    <tol_LSB> \
    [chunk_rows]
Parameters
Argument	Description
input	CSV or DAT file
output.root	Output ROOT file
tol_LSB	Time window in LSB units (1 LSB = 0.5 ns)
chunk_rows	Rows per chunk (default: 300000)
🧪 Example Workflows
🔹 DAT → ROOT (100 ns window)
./csv_or_dat_to_root_events_external_sort_with_ring \
DATA/18_February_2026/Thresh230-5min/Run1_list.dat \
Thresh230-5min_100ns_dat.root \
200 300000

👉 200 LSB = 100 ns

🔹 CSV → ROOT (10 ns window)
../../csv_to_root_events_timestamp_external_sort_with_ring \
../../DATA/4_March_2026/Thresh260/Run1_list.csv \
../../ROOT_files/4_March_2026/Thresh260-5min_10ns_grouping2.root \
20 300000

👉 20 LSB = 10 ns

📊 Histogram Generation: Unique Active Channels per Event (Collision-Free)

After event building, generate histograms:

../../make_hist_df_unique_active_collision_free \
input.root \
output_hist.root \
Events
🎯 Purpose

Analyze detector activity per event:

👉 How many unique channels fired in each event?

🧠 Processing Steps
1. Read Event Data

From TTree Events:

Event_Id
Board_Id
CH_Id
Ring
2. Remove Duplicates (Collision-Free)

Within each event:

(Event_Id, Board_Id, CH_Id, Ring)

Duplicate hits are removed.

👉 Prevents overcounting from repeated hits on same channel.

3. Count Active Channels

For each event:

Count unique active channels
Can be grouped by ring or board
4. Build Histogram
X-axis → Number of active channels per event  
Y-axis → Number of events  
📈 Physical Interpretation
Region	Meaning
Low counts	Noise / single hits
Medium counts	Physical interactions
High counts	Multi-scatter / pile-up
🔬 Why “Collision-Free”?

Without deduplication:

Same channel firing multiple times inflates counts

👉 This step ensures:

True detector occupancy
Correct multiplicity
📁 Output

Produces:

hist_df_*.root

Containing:

Channel multiplicity histograms
Event occupancy distributions
⚠️ Important Note (DAT Input)

Because DAT input lacks correct Board_Id:

❌ Ring-based histograms may be inaccurate
❌ Multi-board distributions unreliable

👉 Histogram is fully reliable only for CSV workflows

🚀 Example Insight
Peak at 2–3 → simple interactions  
Long tail → complex events / pile-up  
🧠 How It Works
Step 1: Chunk Processing
Read input in chunks (chunk_rows)
Compute absolute time
Sort chunks
Write temporary runs
Step 2: External Merge
K-way merge of sorted chunks
Global time ordering
Step 3: Event Building
if (t > event_start + window_us) → new event
Step 4: ROOT Output

TTree: Events

Branch	Description
TStamp_us	Base timestamp
ToA_LSB	Fine timing
t_abs_us	Absolute time
Board_Id	Board index
CH_Id	Channel
Ring	Detector ring
Event_Id	Event grouping
📁 Repository Structure
.
├── csv_or_dat_to_root_events_external_sort_with_ring.cpp
├── csv_to_root_events_timestamp.cpp
├── csv_to_root_events_timestamp_external_sort_with_ring.cpp
├── dat_to_root_events_timestamp.cxx
├── make_hist_df_unique_active_collision_free.cpp
└── README.md
⚡ Performance Notes
Designed for GB-scale datasets
External sort avoids memory overflow
DAT parsing slower due to synchronization scanning
Recommended:
SSD storage
Larger chunk_rows
Avoid network file systems if possible
👤 Author

Pratyush Patel
Northwestern University / Fermilab
SuperCDMS / NEXUS

🚀 Future Improvements
Proper DAT → Board_Id decoding
Faster DAT parsing (optimized scanning)
Parallel chunk sorting
Direct parquet output
Advanced detector topology analysis
