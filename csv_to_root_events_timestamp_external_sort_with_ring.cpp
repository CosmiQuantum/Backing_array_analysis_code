// csv_to_root_events_external_merge_with_ring.cpp
//
// External merge sort (chunk sort + k-way merge) for big Janus Timing CSV,
// with absolute time sorting, Event_Id grouping, and Ring mapping.
//
// Input CSV columns:
//   TStamp_us,Board_Id,Num_Hits,CH_Id,DataType,ToA_LSB,ToT_LSB
// plus leading // metadata lines.
//
// Absolute time:
//   t_abs_us = TStamp_us + 0.0005 * ToA_LSB   (0.5 ns/LSB = 0.0005 us)
//
// Output ROOT:
//   TTree "Events": original fields + t_abs_us + Event_Id + Ring
//
// Usage:
//   ./csv_to_root_events_external_merge_with_ring input.csv output.root tol_LSB [chunk_rows]
//
// Example:
//   ./csv_to_root_events_external_merge_with_ring Run1_list.csv Run1_full.root 200 300000

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using ll = long long;

// -----------------------------
// String helpers
// -----------------------------
static inline std::string trim(const std::string &s) {
  size_t b = 0;
  while (b < s.size() && std::isspace((unsigned char)s[b])) b++;
  size_t e = s.size();
  while (e > b && std::isspace((unsigned char)s[e - 1])) e--;
  return s.substr(b, e - b);
}

static inline bool starts_with(const std::string &s, const std::string &p) {
  return s.rfind(p, 0) == 0;
}

static bool split_csv_7(const std::string &line, std::string out[7]) {
  std::stringstream ss(line);
  std::string item;
  int i = 0;
  while (std::getline(ss, item, ',')) {
    if (i >= 7) return false;
    out[i++] = trim(item);
  }
  return (i == 7);
}

// -----------------------------
// Ring lookup (Board_Id, CH_Id) -> Ring
// Translated from your Python mapping.
// -----------------------------
struct PairHash {
  size_t operator()(const std::pair<int,int>& p) const noexcept {
    return ((size_t)p.first << 32) ^ (size_t)p.second;
  }
};

static std::unordered_map<std::pair<int,int>, int, PairHash> build_ring_lookup() {
  std::unordered_map<std::pair<int,int>, int, PairHash> ring_lookup;
  ring_lookup.reserve(4 * 64);

  const int ring_counter_assignment[] = {
    1, 2, 15, 16, 3, 4, 5, 6,
    17, 18, 19, 20, 7, 8, 9, 10,
    21, 22, 23, 24, 11, 12, 13, 14,
    25, 26, 27, 28
  };
  const int N_ASSIGN =
    (int)(sizeof(ring_counter_assignment) / sizeof(ring_counter_assignment[0]));

  int ring_counter = 0;

  auto add_board = [&](int board_id, const std::vector<std::vector<int>>& ring_lists) {
    for (const auto& ring : ring_lists) {
      if (ring_counter >= N_ASSIGN) return;
      const int ring_id = ring_counter_assignment[ring_counter];
      for (int ch : ring) ring_lookup[{board_id, ch}] = ring_id;
      ring_counter++;
    }
  };

  // board 0 (4 rings)
  add_board(0, {
    {47, 45, 43, 41, 46, 44, 42, 40},
    {63, 61, 59, 57, 62, 60, 58, 56},
    {39, 37, 35, 33, 38, 36, 34, 32},
    {55, 53, 51, 49, 54, 52, 50, 48}
  });

  // common 8-ring pattern for boards 1,2,3
  const std::vector<std::vector<int>> common = {
    {15, 14, 13, 12, 11, 10,  9,  8},
    {31, 30, 29, 28, 27, 26, 25, 24},
    {47, 46, 45, 44, 43, 42, 41, 40},
    {63, 62, 61, 60, 59, 58, 57, 56},
    { 7,  6,  5,  4,  3,  2,  1,  0},
    {23, 22, 21, 20, 19, 18, 17, 16},
    {39, 38, 37, 36, 35, 34, 33, 32},
    {55, 54, 53, 52, 51, 50, 49, 48}
  };

  add_board(1, common);
  add_board(2, common);
  add_board(3, common);

  return ring_lookup;
}

static inline int get_ring(
  const std::unordered_map<std::pair<int,int>, int, PairHash>& ring_lookup,
  int board_id,
  int ch_id
) {
  auto it = ring_lookup.find({board_id, ch_id});
  if (it == ring_lookup.end()) return -1;
  return it->second;
}

// -----------------------------
// Row structure
// -----------------------------
struct Row {
  double TStamp_us = 0.0;
  ll Board_Id = 0;
  ll Num_Hits = 0;
  ll CH_Id = 0;
  ll DataType = 0;  // hex like 0x30
  ll ToA_LSB = 0;
  ll ToT_LSB = 0;
  double t_abs_us = 0.0;
  int Ring = -1;
};

// -----------------------------
// CSV parse: compute abs time + ring
// -----------------------------
static bool parse_csv_row(
  const std::string &raw_line,
  Row &r,
  const std::unordered_map<std::pair<int,int>, int, PairHash>& ring_lookup
) {
  std::string line = raw_line;
  if (!line.empty() && line.back() == '\r') line.pop_back(); // Windows CR
  line = trim(line);

  if (line.empty()) return false;
  if (starts_with(line, "//")) return false;
  if (line.find("TStamp_us") != std::string::npos) return false;

  std::string f[7];
  if (!split_csv_7(line, f)) return false;

  try {
    r.TStamp_us = std::stod(f[0]);
    r.Board_Id  = std::stoll(f[1]);
    r.Num_Hits  = std::stoll(f[2]);
    r.CH_Id     = std::stoll(f[3]);
    r.DataType  = std::stoll(f[4], nullptr, 0); // base=0 supports 0x..
    r.ToA_LSB   = std::stoll(f[5]);
    r.ToT_LSB   = std::stoll(f[6]);
  } catch (...) {
    return false;
  }

  r.t_abs_us = r.TStamp_us + 0.0005 * (double)r.ToA_LSB;
  r.Ring = get_ring(ring_lookup, (int)r.Board_Id, (int)r.CH_Id);
  return true;
}

// -----------------------------
// Run file I/O (temporary sorted chunks)
// Store Ring so PASS 2 doesn't need lookup.
// -----------------------------
static bool write_run(const std::string &path, const std::vector<Row> &chunk) {
  std::ofstream out(path);
  if (!out.is_open()) return false;

  // Format (9 cols):
  // t_abs_us,TStamp_us,Board_Id,Num_Hits,CH_Id,DataType,ToA_LSB,ToT_LSB,Ring
  out.setf(std::ios::fixed);
  out.precision(9);

  for (const auto &r : chunk) {
    out << r.t_abs_us << ","
        << r.TStamp_us << ","
        << r.Board_Id  << ","
        << r.Num_Hits  << ","
        << r.CH_Id     << ","
        << r.DataType  << ","
        << r.ToA_LSB   << ","
        << r.ToT_LSB   << ","
        << r.Ring      << "\n";
  }
  return true;
}

static bool read_run_row(std::ifstream &in, Row &r) {
  std::string line;
  if (!std::getline(in, line)) return false;
  if (!line.empty() && line.back() == '\r') line.pop_back();

  std::stringstream ss(line);
  std::string item;
  std::vector<std::string> f;
  while (std::getline(ss, item, ',')) f.push_back(trim(item));
  if (f.size() < 9) return false;

  try {
    r.t_abs_us  = std::stod(f[0]);
    r.TStamp_us = std::stod(f[1]);
    r.Board_Id  = std::stoll(f[2]);
    r.Num_Hits  = std::stoll(f[3]);
    r.CH_Id     = std::stoll(f[4]);
    r.DataType  = std::stoll(f[5], nullptr, 0);
    r.ToA_LSB   = std::stoll(f[6]);
    r.ToT_LSB   = std::stoll(f[7]);
    r.Ring      = std::stoi(f[8]);
  } catch (...) {
    return false;
  }
  return true;
}

// -----------------------------
// K-way merge heap
// -----------------------------
struct HeapItem {
  Row row;
  size_t run_idx = 0;
};

struct HeapCmp {
  bool operator()(const HeapItem &a, const HeapItem &b) const {
    return a.row.t_abs_us > b.row.t_abs_us; // min-heap behavior
  }
};

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr
      << "Usage:\n"
      << "  " << argv[0] << " <input.csv> <output.root> <tol_LSB> [chunk_rows]\n\n"
      << "Example:\n"
      << "  " << argv[0] << " Run1_list.csv out.root 200 500000\n";
    return 1;
  }

  const std::string input_csv   = argv[1];
  const std::string output_root = argv[2];
  const ll tol_LSB = std::stoll(argv[3]);

  // chunk_rows controls memory usage. Start with 200k–1M depending on RAM.
  const size_t chunk_rows = (argc >= 5) ? (size_t)std::stoull(argv[4]) : 300000;

  const double window_us = (double)tol_LSB * 0.0005;

  std::cerr << "Input:  " << input_csv << "\n";
  std::cerr << "Output: " << output_root << "\n";
  std::cerr << "tol_LSB: " << tol_LSB << "  => window_us=" << window_us << " us\n";
  std::cerr << "chunk_rows: " << chunk_rows << "\n";

  const auto ring_lookup = build_ring_lookup();

  // ---------- PASS 1: Create sorted runs ----------
  std::ifstream fin(input_csv);
  if (!fin.is_open()) {
    std::cerr << "ERROR: cannot open input: " << input_csv << "\n";
    return 2;
  }

  std::vector<std::string> run_files;
  run_files.reserve(128);

  std::vector<Row> chunk;
  chunk.reserve(chunk_rows);

  std::string line;
  size_t total_rows = 0;
  size_t run_count = 0;
  size_t bad_lines = 0;

  while (std::getline(fin, line)) {
    Row r;
    if (!parse_csv_row(line, r, ring_lookup)) {
      std::string t = trim(line);
      if (!t.empty() && !starts_with(t, "//") && t.find("TStamp_us") == std::string::npos)
        bad_lines++;
      continue;
    }

    chunk.push_back(r);
    total_rows++;

    if (chunk.size() >= chunk_rows) {
      std::sort(chunk.begin(), chunk.end(),
                [](const Row &a, const Row &b){ return a.t_abs_us < b.t_abs_us; });

      std::string run_path = "run_" + std::to_string(run_count++) + ".csv";
      if (!write_run(run_path, chunk)) {
        std::cerr << "ERROR: failed to write run: " << run_path << "\n";
        return 3;
      }
      run_files.push_back(run_path);
      chunk.clear();

      std::cerr << "Wrote run " << run_files.back()
                << "  (rows so far: " << total_rows << ")\n";
    }
  }
  fin.close();

  // flush final partial chunk
  if (!chunk.empty()) {
    std::sort(chunk.begin(), chunk.end(),
              [](const Row &a, const Row &b){ return a.t_abs_us < b.t_abs_us; });

    std::string run_path = "run_" + std::to_string(run_count++) + ".csv";
    if (!write_run(run_path, chunk)) {
      std::cerr << "ERROR: failed to write run: " << run_path << "\n";
      return 3;
    }
    run_files.push_back(run_path);
    chunk.clear();

    std::cerr << "Wrote run " << run_files.back()
              << "  (rows total: " << total_rows << ")\n";
  }

  if (run_files.empty()) {
    std::cerr << "ERROR: no valid rows parsed. bad_lines=" << bad_lines << "\n";
    return 4;
  }

  std::cerr << "Total parsed rows: " << total_rows << "\n";
  std::cerr << "Runs created: " << run_files.size() << "\n";
  if (bad_lines) std::cerr << "Warning: bad data lines skipped: " << bad_lines << "\n";

  // ---------- PASS 2: K-way merge + event grouping + ROOT output ----------
  std::vector<std::ifstream> run_streams(run_files.size());
  for (size_t i = 0; i < run_files.size(); i++) {
    run_streams[i].open(run_files[i]);
    if (!run_streams[i].is_open()) {
      std::cerr << "ERROR: cannot open run file: " << run_files[i] << "\n";
      return 5;
    }
  }

  std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> pq;

  // prime heap with first row from each run
  for (size_t i = 0; i < run_streams.size(); i++) {
    Row r;
    if (read_run_row(run_streams[i], r)) {
      pq.push(HeapItem{r, i});
    }
  }

  TFile *fout = TFile::Open(output_root.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create ROOT file: " << output_root << "\n";
    return 6;
  }

  TTree *tree = new TTree("Events", "Events grouped by abs time (external merge) + Ring");

  double TStamp_us = 0.0;
  ll Board_Id = 0;
  ll Num_Hits = 0;
  ll CH_Id = 0;
  ll DataType = 0;
  ll ToA_LSB = 0;
  ll ToT_LSB = 0;
  double t_abs_us = 0.0;
  int Ring = -1;
  ll Event_Id = 0;

  tree->Branch("TStamp_us", &TStamp_us);
  tree->Branch("Board_Id", &Board_Id);
  tree->Branch("Num_Hits", &Num_Hits);
  tree->Branch("CH_Id", &CH_Id);
  tree->Branch("DataType", &DataType);
  tree->Branch("ToA_LSB", &ToA_LSB);
  tree->Branch("ToT_LSB", &ToT_LSB);
  tree->Branch("t_abs_us", &t_abs_us);
  tree->Branch("Ring", &Ring);
  tree->Branch("Event_Id", &Event_Id);

  //bool have_last = false;
  //double last_hit_us = 0.0;
  //ll current_event_id = 1;

  //while (!pq.empty()) {
  //  HeapItem top = pq.top();
  //  pq.pop();

  //  const Row &r = top.row;

    // Event grouping on globally sorted stream
   // if (!have_last) {
   //   current_event_id = 1;
   //   last_hit_us = r.t_abs_us;
   //   have_last = true;
   // } else {
   //   if (r.t_abs_us > last_hit_us + window_us) {
   //     current_event_id++;
   //   }
   //   last_hit_us = r.t_abs_us;
   // }
   bool have_event = false;
   double event_start_us = 0.0;
   ll current_event_id = 1;

   while (!pq.empty()) {
     HeapItem top = pq.top();
     pq.pop();

     const Row &r = top.row;

     if (!have_event) {
       current_event_id = 1;
       event_start_us = r.t_abs_us;   // anchor at first hit
       have_event = true;
     } else {
     // start a NEW event only when we cross the window from the FIRST hit
     if (r.t_abs_us > event_start_us + window_us) {
       current_event_id++;
       event_start_us = r.t_abs_us; // re-anchor to this hit
      }
    }

  // assign event id to this row (wherever you store it)
  // out_row.Event_Id = current_event_id;

    TStamp_us = r.TStamp_us;
    Board_Id  = r.Board_Id;
    Num_Hits  = r.Num_Hits;
    CH_Id     = r.CH_Id;
    DataType  = r.DataType;
    ToA_LSB   = r.ToA_LSB;
    ToT_LSB   = r.ToT_LSB;
    t_abs_us  = r.t_abs_us;
    Ring      = r.Ring;
    Event_Id  = current_event_id;

    tree->Fill();

    // pull next row from same run and push
    Row next;
    if (read_run_row(run_streams[top.run_idx], next)) {
      pq.push(HeapItem{next, top.run_idx});
    }
  }

  const auto nEntries = tree->GetEntries();

  std::cerr << "Done. ROOT entries: " << nEntries << "\n";
  std::cerr << "Last Event_Id: " << current_event_id << "\n";

  fout->Write();
  fout->Close();

  // Cleanup temp run files
  for (const auto &p : run_files) std::remove(p.c_str());

  return 0;
}
