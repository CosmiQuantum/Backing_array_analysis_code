# Backing_array_analysis_code

## Event Building & External Sorting for Timing Data (CSV / DAT to ROOT)

This repository provides high-performance C++ tools to process large-scale timing data from CAEN/Janus systems and convert large CSV or binary DAT files into structured ROOT files for physics analysis.

The main workflow performs:

- Absolute timestamp reconstruction
- External merge sort for large datasets
- Event building using time-window clustering
- Ring mapping using `(Board_Id, CH_Id) -> Ring`
- ROOT TTree output
- Histogram generation for detector activity studies

---

## Features

- Supports CSV and binary DAT input
- Handles very large files using external sort with chunk sorting and k-way merge
- Computes absolute time using in us:

```text
t_abs_us = TStamp_us + 0.0005 * ToA_LSB
```

- Builds events using a configurable time window
- Adds detector ring information for geometry-based analysis
- Writes output to ROOT TTrees
- Supports downstream histogram generation for event multiplicity and occupancy studies

---

## Important Limitation: DAT Input and Ring Identification

### Current Status

Ring identification does not currently work correctly for DAT input.

Ring mapping works correctly for CSV input.

### Reason

The DAT format currently parsed by this repository does not include, or does not yet decode, `Board_Id`.

In the current implementation, the board ID is hardcoded:

```cpp
r.Board_Id = default_board_id;  // hardcoded, currently 0
```

As a result:

- All DAT hits are assigned to a single board
- The mapping `(Board_Id, CH_Id) -> Ring` becomes incorrect
- Multi-board detector geometry is not preserved

### Impact

For DAT-based workflows:

- Ring assignments are unreliable
- Multi-board detector geometry is lost
- Ring-wise histograms may be incorrect
- Coincidence and spatial analysis may be unreliable

### Recommended Usage

| Use Case | Recommended Input |
|---|---|
| Full ring-based physics analysis | CSV |
| Timing validation | DAT |
| Event grouping only | DAT |
| Ring-wise histograms | CSV |

### Workarounds

Use CSV input for any analysis that requires correct ring information.

Use DAT input only for:

- Timing validation
- Event clustering checks
- Debugging raw hit structure
- Non-ring-based studies

### Planned Improvements

- Decode `Board_Id` from DAT if present in the binary structure
- Support multi-board DAT reconstruction
- Improve DAT parsing speed and robustness

---

## Requirements

- C++17 compiler
- ROOT 6.x

---

## Compilation

Compile the main CSV/DAT event builder with:

```bash
g++ -O2 -std=c++17 csv_or_dat_to_root_events_external_sort_with_ring.cpp \
    $(root-config --cflags --libs) \
    -o csv_or_dat_to_root_events_external_sort_with_ring
```

---

## Usage

```bash
./csv_or_dat_to_root_events_external_sort_with_ring \
    <input.csv|input.dat> \
    <output.root> \
    <tol_LSB> \
    [chunk_rows]
```

### Arguments

| Argument | Description |
|---|---|
| `input.csv` or `input.dat` | Input timing data file |
| `output.root` | Output ROOT file |
| `tol_LSB` | Event-building time window in LSB units |
| `chunk_rows` | Number of rows per chunk for external sorting |

One LSB corresponds to 0.5 ns.

Therefore:

```text
window_us = tol_LSB * 0.0005
```

Examples:

| `tol_LSB` | Time Window |
|---|---|
| 20 | 10 ns |
| 200 | 100 ns |

---

## Example Workflows

### DAT to ROOT with 100 ns event window

```bash
./csv_or_dat_to_root_events_external_sort_with_ring \
    DATA/18_February_2026/Thresh230-5min/Run1_list.dat \
    Thresh230-5min_100ns_dat.root \
    200 \
    300000
```

Here:

```text
200 LSB = 100 ns
```

### CSV to ROOT with 10 ns event window

```bash
./csv_to_root_events_timestamp_external_sort_with_ring \
    ../../DATA/4_March_2026/Thresh260/Run1_list.csv \
    ../../ROOT_files/4_March_2026/Thresh260-5min_10ns_grouping2.root \
    20 \
    300000
```

Here:

```text
20 LSB = 10 ns
```

---

## Histogram Generation: Unique Active Channels per Event

After event building, histograms can be generated using:

```bash
./make_hist_df_unique_active_collision_free \
    input.root \
    output_hist.root \
    Events
```

Example:

```bash
./make_hist_df_unique_active_collision_free \
    ../../ROOT_files/4_March_2026/Thresh230-5min_100ns_grouping2.root \
    hist_df_Thresh230-5min_100ns_grouping2.root \
    Events
```

---

## Purpose of the Histogram Step

The histogram step analyzes detector activity per event.

The main question is:

```text
How many unique detector channels fired in each event?
```

This is useful for studying:

- Event multiplicity
- Detector occupancy
- Noise-like events
- Multi-scatter events
- Pile-up behavior
- Ring-wise or board-wise detector activity

---

## Histogram Processing Steps

### 1. Read Event Data

The histogram code reads the `Events` TTree from the ROOT file.

Relevant branches include:

- `Event_Id`
- `Board_Id`
- `CH_Id`
- `Ring`

### 2. Remove Duplicate Hits

Within each event, hits are deduplicated using:

```text
(Event_Id, Board_Id, CH_Id, Ring)
```

This makes the counting collision-free.

A single channel may fire multiple times in one event. Without deduplication, the same physical channel could be counted multiple times.

### 3. Count Unique Active Channels

For each event, the code counts the number of unique active channels.

This can be studied globally or grouped by board or ring.

### 4. Build Histograms

The output histogram represents:

```text
X-axis: Number of active channels per event
Y-axis: Number of events
```

---

## Physical Interpretation of Histograms

| Region | Interpretation |
|---|---|
| Low active-channel counts | Noise-like events or single-site hits |
| Medium active-channel counts | Typical physical interactions |
| High active-channel counts | Multi-scatter events, pile-up, or shower-like activity |

Example interpretation:

```text
Peak at 2-3 channels: mostly simple interactions
Long high-multiplicity tail: complex events or pile-up
```

---

## Why Collision-Free Counting Matters

Without deduplication, a single channel firing multiple times can inflate the channel multiplicity.

Collision-free counting gives a better estimate of:

- True detector occupancy
- Spatial multiplicity
- Event topology
- Ring-wise response

---

## Histogram Output

The histogram step produces a ROOT file such as:

```text
hist_df_*.root
```

The output may contain:

- Channel multiplicity histograms
- Event occupancy distributions
- Ring-wise or board-wise activity histograms, depending on the implementation

---

## Important Note for DAT Input Histograms

Because DAT input currently does not preserve or decode `Board_Id` correctly:

- Ring-based histograms from DAT input may be inaccurate
- Multi-board distributions from DAT input are unreliable
- Total multiplicity may still be useful for timing or debugging studies

For reliable ring-wise histograms, use CSV input.

---

## How the Event Builder Works

### Step 1: Chunk Processing

The input file is read in chunks.

For each row or hit:

1. Read timestamp fields
2. Compute absolute time
3. Assign ring information if possible
4. Store rows in memory until `chunk_rows` is reached

Each chunk is sorted by absolute time and written to a temporary run file.

### Step 2: External Merge

The sorted temporary run files are merged using a k-way merge.

This produces a globally time-ordered stream of hits without loading the full dataset into memory.

### Step 3: Event Building

Events are formed using a fixed time window.

Conceptually:

```text
if current_hit_time > event_start_time + window_us:
    start a new event
else:
    assign hit to current event
```

### Step 4: ROOT Output

The output ROOT file contains a TTree named:

```text
Events
```

Typical branches include:

| Branch | Description |
|---|---|
| `TStamp_us` | Base timestamp |
| `Board_Id` | Board index |
| `Num_Hits` | Number of hits in original block or event group |
| `CH_Id` | Channel ID |
| `DataType` | Data type field |
| `ToA_LSB` | Time-of-arrival in LSB units |
| `ToT_LSB` | Time-over-threshold in LSB units |
| `t_abs_us` | Reconstructed absolute time |
| `Ring` | Detector ring |
| `Event_Id` | Event identifier |

---

## Repository Structure

```text
.
├── csv_or_dat_to_root_events_external_sort_with_ring.cpp
├── csv_to_root_events_timestamp.cpp
├── csv_to_root_events_timestamp_external_sort_with_ring.cpp
├── dat_to_root_events_timestamp.cxx
├── make_hist_df_unique_active_collision_free.cpp
└── README.md
```

---

## Performance Notes

This repository is designed for large datasets where the full input file may not fit comfortably in memory.

External sorting avoids memory overflow by:

1. Sorting smaller chunks
2. Writing temporary sorted run files
3. Merging the sorted runs

### Recommended Practices

- Use SSD storage when possible
- Avoid slow network file systems for temporary files
- Increase `chunk_rows` if enough memory is available
- Use smaller `chunk_rows` if memory usage becomes too high
- Keep temporary files on a fast local disk when possible

### DAT Performance Note

DAT parsing can be slower because the parser may need to synchronize to valid binary blocks.

Future optimized DAT parsing may improve performance significantly.

---

## Future Improvements

Planned or possible improvements include:

- Proper DAT to `Board_Id` decoding
- Faster DAT parsing with optimized scanning
- Parallel chunk sorting
- Direct parquet output
- Automatic ROOT plotting scripts
- Ring-wise occupancy plots
- Detector topology classification
- Batch-job workflow integration for Fermigrid or similar systems

---

## Author

Pratyush Patel, Northwestern University
