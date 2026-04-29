// csv_or_dat_to_root_events_external_merge_with_ring.cpp
//
// External merge sort (chunk sort + k-way merge) for big Janus Timing data,
// with absolute time sorting, Event_Id grouping, and Ring mapping.
//
// Supports input:
//   - CSV/TXT (lines): TStamp_us,Board_Id,Num_Hits,CH_Id,DataType,ToA_LSB,ToT_LSB
//     plus leading // metadata lines.
//
//   - DAT (binary): inferred from your hexdump:
//       [ ... file header ... ]
//       repeated blocks:
//         double   TStamp_us       (little-endian IEEE754)
//         uint16   Num_Hits
//         repeated Num_Hits hits, each 8 bytes:
//           uint8   CH_Id
//           uint8   DataType       (commonly 0x30)
//           uint32  ToA_LSB
//           uint16  ToT_LSB
//
// Absolute time:
//   t_abs_us = TStamp_us + 0.0005 * ToA_LSB   (0.5 ns/LSB = 0.0005 us)
//
// Output ROOT:
//   TTree "Events": original fields + t_abs_us + Event_Id + Ring
//
// Usage:
//   ./csv_or_dat_to_root_events_external_merge_with_ring input.(csv|txt|dat) output.root tol_LSB [chunk_rows]
//
// Example:
//   ./csv_or_dat_to_root_events_external_merge_with_ring Run1_list.csv out.root 200 300000
//   ./csv_or_dat_to_root_events_external_merge_with_ring Run1_list.dat out.root 200 300000

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
#include <cmath>
#include <cstdint>

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
// -----------------------------
struct PairHash {
  size_t operator()(const std::pair<int,int>& p) const noexcept {
    return ((size_t)(uint32_t)p.first << 32) ^ (size_t)(uint32_t)p.second;
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
  ll DataType = 0;
  ll ToA_LSB = 0;
  ll ToT_LSB = 0;
  double t_abs_us = 0.0;
  int Ring = -1;
};

// -----------------------------
// CSV parse
// -----------------------------
static bool parse_csv_row(
  const std::string &raw_line,
  Row &r,
  const std::unordered_map<std::pair<int,int>, int, PairHash>& ring_lookup
) {
  std::string line = raw_line;
  if (!line.empty() && line.back() == '\r') line.pop_back();
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
    r.DataType  = std::stoll(f[4], nullptr, 0);
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
// Temp run file I/O
// -----------------------------
static bool write_run(const std::string &path, const std::vector<Row> &chunk) {
  std::ofstream out(path);
  if (!out.is_open()) return false;

  out.setf(std::ios::fixed);
  out.precision(9);

  // t_abs_us,TStamp_us,Board_Id,Num_Hits,CH_Id,DataType,ToA_LSB,ToT_LSB,Ring
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
// DAT parsing (robust sync + validate block + validate next block)
// -----------------------------
#pragma pack(push, 1)
struct DatHit {
  uint8_t  ch;     // 1
  uint8_t  dtype;  // 1
  uint32_t toa;    // 4
  uint16_t tot;    // 2
};
#pragma pack(pop)

static inline bool plausible_tstamp(double t_us) {
  return std::isfinite(t_us) && t_us >= 0.0 && t_us <= 1.0e15; // very loose
}

static inline bool plausible_numhits(uint16_t n) {
  // Your dumps show ~30. Typical Janus blocks may be <= 4096.
  return n >= 1 && n <= 8192;
}

static inline bool plausible_hit(const DatHit& h) {
  // Strong constraints prevent false sync:
  // - channels are 0..63 (DT5202)
  // - dtype frequently 0x30 (as in your hexdump)
  // - ToT is typically not huge (fits 16-bit but usually < few thousand)
  if (h.ch > 63) return false;
  if (h.dtype != 0x30) return false;
  if (h.tot > 20000) return false; // loose upper guard
  // toa any 32-bit is okay; but guard absurd:
  // if (h.toa > 0xF0000000u) return false;
  return true;
}

// Try to validate a block at an absolute file offset `pos`.
// If valid, returns true and fills (tstamp_us, num_hits).
static bool dat_is_valid_block_at(std::ifstream& fin, std::streampos pos,
                                  double& tstamp_us, uint16_t& num_hits) {
  fin.clear();
  fin.seekg(pos);
  if (!fin.good()) return false;

  double t = 0.0;
  uint16_t nh = 0;

  fin.read(reinterpret_cast<char*>(&t), sizeof(double));
  fin.read(reinterpret_cast<char*>(&nh), sizeof(uint16_t));
  if (!fin.good()) return false;

  if (!plausible_tstamp(t) || !plausible_numhits(nh)) return false;

  // Validate first few hits (up to 4) for strong signature
  const int K = 4;
  int ncheck = std::min<int>(K, (int)nh);
  for (int i = 0; i < ncheck; i++) {
    DatHit h{};
    fin.read(reinterpret_cast<char*>(&h), sizeof(DatHit));
    if (!fin.good()) return false;
    if (!plausible_hit(h)) return false;
  }

  // Also validate that the *next* block header exists at expected end offset
  const std::streampos next_pos = pos + (std::streamoff)(sizeof(double) + sizeof(uint16_t) + (std::streamoff)nh * (std::streamoff)sizeof(DatHit));
  fin.clear();
  fin.seekg(next_pos);
  if (!fin.good()) return false;

  // Peek next header (if we're near EOF this might fail; allow EOF as "valid last block")
  double t2 = 0.0;
  uint16_t nh2 = 0;
  fin.read(reinterpret_cast<char*>(&t2), sizeof(double));
  fin.read(reinterpret_cast<char*>(&nh2), sizeof(uint16_t));
  if (!fin.good()) {
    // Could just be the last block; accept.
    tstamp_us = t;
    num_hits = nh;
    return true;
  }
  if (!plausible_tstamp(t2) || !plausible_numhits(nh2)) return false;

  // Peek 1st hit of next block
  DatHit h2{};
  fin.read(reinterpret_cast<char*>(&h2), sizeof(DatHit));
  if (!fin.good()) return false;
  if (!plausible_hit(h2)) return false;

  tstamp_us = t;
  num_hits = nh;
  return true;
}

// Find next valid block by scanning forward byte-by-byte.
// Leaves stream positioned at start of block if found.
static bool dat_sync_to_next_block(std::ifstream& fin) {
  // We scan up to N bytes each call to avoid infinite loops on corruption.
  const size_t MAX_SCAN = 64 * 1024 * 1024; // 64 MB

  std::streampos start = fin.tellg();
  if (start < 0) return false;

  for (size_t step = 0; step < MAX_SCAN; step++) {
    std::streampos pos = start + (std::streamoff)step;

    double t = 0.0;
    uint16_t nh = 0;
    if (dat_is_valid_block_at(fin, pos, t, nh)) {
      fin.clear();
      fin.seekg(pos);
      return fin.good();
    }
  }
  return false;
}

// Read ONE block at current position (assumes positioned at exact block start).
// Returns false on EOF/bad.
static bool dat_read_block_here(std::ifstream& fin, double& tstamp_us, uint16_t& num_hits, std::vector<DatHit>& hits) {
  hits.clear();
  double t = 0.0;
  uint16_t nh = 0;

  fin.read(reinterpret_cast<char*>(&t), sizeof(double));
  fin.read(reinterpret_cast<char*>(&nh), sizeof(uint16_t));
  if (!fin.good()) return false;

  hits.resize(nh);
  fin.read(reinterpret_cast<char*>(hits.data()), (std::streamsize)nh * (std::streamsize)sizeof(DatHit));
  if (!fin.good()) return false;

  tstamp_us = t;
  num_hits = nh;
  return true;
}

// Reads next synced block and appends hits as Rows
static bool read_dat_block_strict(
  std::ifstream& fin,
  std::vector<Row>& out_rows,
  const std::unordered_map<std::pair<int,int>, int, PairHash>& ring_lookup,
  int default_board_id
) {
  // If we're at EOF, stop
  fin.peek();
  if (!fin.good()) return false;

  // Sync to next block start
  if (!dat_sync_to_next_block(fin)) return false;

  // Now read block at this position
  double tstamp_us = 0.0;
  uint16_t num_hits = 0;
  std::vector<DatHit> hits;
  if (!dat_read_block_here(fin, tstamp_us, num_hits, hits)) return false;

  // Convert hits -> Rows
  for (uint16_t i = 0; i < num_hits; i++) {
    const DatHit& h = hits[i];

    // If the block is valid, these should hold, but keep safe:
    if (!plausible_hit(h)) continue;

    Row r;
    r.TStamp_us = tstamp_us;
    r.Board_Id  = default_board_id;   // If your DAT encodes board id elsewhere, update this
    r.Num_Hits  = (ll)num_hits;
    r.CH_Id     = (ll)h.ch;
    r.DataType  = (ll)h.dtype;
    r.ToA_LSB   = (ll)h.toa;
    r.ToT_LSB   = (ll)h.tot;

    r.t_abs_us  = r.TStamp_us + 0.0005 * (double)r.ToA_LSB;
    r.Ring      = get_ring(ring_lookup, (int)r.Board_Id, (int)r.CH_Id);

    out_rows.push_back(r);
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
    return a.row.t_abs_us > b.row.t_abs_us;
  }
};

static bool has_extension_ci(const std::string& path, const std::string& ext) {
  if (path.size() < ext.size()) return false;
  auto tail = path.substr(path.size() - ext.size());
  auto lower = [](std::string s){
    for (auto &c : s) c = (char)std::tolower((unsigned char)c);
    return s;
  };
  return lower(tail) == lower(ext);
}

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr
      << "Usage:\n"
      << "  " << argv[0] << " <input.(csv|txt|dat)> <output.root> <tol_LSB> [chunk_rows]\n";
    return 1;
  }

  const std::string input_file  = argv[1];
  const std::string output_root = argv[2];
  const ll tol_LSB = std::stoll(argv[3]);
  const size_t chunk_rows = (argc >= 5) ? (size_t)std::stoull(argv[4]) : 300000;

  const double window_us = (double)tol_LSB * 0.0005;

  const bool is_dat = has_extension_ci(input_file, ".dat");
  std::cerr << "Input:  " << input_file << (is_dat ? " (DAT binary)\n" : " (CSV/TXT)\n");
  std::cerr << "Output: " << output_root << "\n";
  std::cerr << "tol_LSB: " << tol_LSB << "  => window_us=" << window_us << " us\n";
  std::cerr << "chunk_rows: " << chunk_rows << "\n";

  const auto ring_lookup = build_ring_lookup();

  // ---------- PASS 1: Create sorted runs ----------
  std::ifstream fin;
  if (is_dat) fin.open(input_file, std::ios::binary);
  else        fin.open(input_file);

  if (!fin.is_open()) {
    std::cerr << "ERROR: cannot open input: " << input_file << "\n";
    return 2;
  }

  std::vector<std::string> run_files;
  run_files.reserve(128);

  std::vector<Row> chunk;
  chunk.reserve(chunk_rows);

  size_t total_rows = 0;
  size_t run_count  = 0;
  size_t bad_lines  = 0;

  if (is_dat) {
    // NOTE: your DAT dump doesn't show explicit board_id in block; using 0.
    // If you want board_id, we can parse it once you show where it is.
    const int default_board_id = 0;

    while (true) {
      const size_t before = chunk.size();
      if (!read_dat_block_strict(fin, chunk, ring_lookup, default_board_id)) break;
      const size_t added = chunk.size() - before;
      total_rows += added;

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
  } else {
    std::string line;
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
  if (!is_dat && bad_lines) std::cerr << "Warning: bad data lines skipped: " << bad_lines << "\n";

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

  bool have_event = false;
  double event_start_us = 0.0;
  ll current_event_id = 1;

  while (!pq.empty()) {
    HeapItem top = pq.top();
    pq.pop();

    const Row &r = top.row;

    if (!have_event) {
      current_event_id = 1;
      event_start_us = r.t_abs_us;
      have_event = true;
    } else {
      if (r.t_abs_us > event_start_us + window_us) {
        current_event_id++;
        event_start_us = r.t_abs_us;
      }
    }

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

    Row next;
    if (read_run_row(run_streams[top.run_idx], next)) {
      pq.push(HeapItem{next, top.run_idx});
    }
  }

  std::cerr << "Done. ROOT entries: " << tree->GetEntries() << "\n";
  std::cerr << "Last Event_Id: " << current_event_id << "\n";

  fout->Write();
  fout->Close();

  for (const auto &p : run_files) std::remove(p.c_str());

  return 0;
}
