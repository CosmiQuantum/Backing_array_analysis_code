#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using ll = long long;

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

struct Row {
  double TStamp_us = 0.0;
  ll Board_Id = 0;
  ll Num_Hits = 0;
  ll CH_Id = 0;
  ll DataType = 0;  // hex like 0x30
  ll ToA_LSB = 0;
  ll ToT_LSB = 0;
  double t_abs_us = 0.0;
};

static bool parse_csv_row(const std::string &raw_line, Row &r) {
  std::string line = raw_line;
  // remove Windows CR if present
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
    r.DataType  = std::stoll(f[4], nullptr, 0); // base=0 supports 0x..
    r.ToA_LSB   = std::stoll(f[5]);
    r.ToT_LSB   = std::stoll(f[6]);
  } catch (...) {
    return false;
  }

  // abs time in microseconds: 0.5 ns/LSB = 0.0005 us/LSB
  r.t_abs_us = r.TStamp_us + 0.0005 * (double)r.ToA_LSB;
  return true;
}

// Write one sorted run to disk (simple text format)
static bool write_run(const std::string &path, const std::vector<Row> &chunk) {
  std::ofstream out(path);
  if (!out.is_open()) return false;

  // Format: t_abs_us,TStamp_us,Board_Id,Num_Hits,CH_Id,DataType,ToA_LSB,ToT_LSB
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
        << r.ToT_LSB   << "\n";
  }
  return true;
}

// Read next Row from a run file
static bool read_run_row(std::ifstream &in, Row &r) {
  std::string line;
  if (!std::getline(in, line)) return false;
  if (!line.empty() && line.back() == '\r') line.pop_back();

  // Expect 8 columns now (t_abs_us first)
  std::stringstream ss(line);
  std::string item;
  std::vector<std::string> f;
  while (std::getline(ss, item, ',')) f.push_back(trim(item));
  if (f.size() < 8) return false;

  try {
    r.t_abs_us  = std::stod(f[0]);
    r.TStamp_us = std::stod(f[1]);
    r.Board_Id  = std::stoll(f[2]);
    r.Num_Hits  = std::stoll(f[3]);
    r.CH_Id     = std::stoll(f[4]);
    r.DataType  = std::stoll(f[5], nullptr, 0);
    r.ToA_LSB   = std::stoll(f[6]);
    r.ToT_LSB   = std::stoll(f[7]);
  } catch (...) {
    return false;
  }
  return true;
}

struct HeapItem {
  Row row;
  size_t run_idx = 0;
};

struct HeapCmp {
  bool operator()(const HeapItem &a, const HeapItem &b) const {
    // min-heap by t_abs_us -> priority_queue is max-heap by default, so invert
    return a.row.t_abs_us > b.row.t_abs_us;
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

  const std::string input_csv  = argv[1];
  const std::string output_root = argv[2];
  const ll tol_LSB = std::stoll(argv[3]);

  // chunk_rows controls memory usage. Start with 200k–1M depending on RAM.
  const size_t chunk_rows = (argc >= 5) ? (size_t)std::stoull(argv[4]) : 300000;

  const double window_us = (double)tol_LSB * 0.0005;

  std::cerr << "Input:  " << input_csv << "\n";
  std::cerr << "Output: " << output_root << "\n";
  std::cerr << "tol_LSB: " << tol_LSB << "  => window_us=" << window_us << " us\n";
  std::cerr << "chunk_rows: " << chunk_rows << "\n";

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
    if (!parse_csv_row(line, r)) {
      // Only count as bad if it's not an ignorable header/comment/blank line
      std::string t = trim(line);
      if (!t.empty() && !starts_with(t, "//") && t.find("TStamp_us")==std::string::npos)
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

      std::cerr << "Wrote run " << run_files.back() << "  (rows so far: " << total_rows << ")\n";
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

    std::cerr << "Wrote run " << run_files.back() << "  (rows total: " << total_rows << ")\n";
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

  TTree *tree = new TTree("Events", "Events grouped by abs time (external merge)");

  double TStamp_us = 0.0;
  ll Board_Id = 0;
  ll Num_Hits = 0;
  ll CH_Id = 0;
  ll DataType = 0;
  ll ToA_LSB = 0;
  ll ToT_LSB = 0;
  double t_abs_us = 0.0;
  ll Event_Id = 0;

  tree->Branch("TStamp_us", &TStamp_us);
  tree->Branch("Board_Id", &Board_Id);
  tree->Branch("Num_Hits", &Num_Hits);
  tree->Branch("CH_Id", &CH_Id);
  tree->Branch("DataType", &DataType);
  tree->Branch("ToA_LSB", &ToA_LSB);
  tree->Branch("ToT_LSB", &ToT_LSB);
  tree->Branch("t_abs_us", &t_abs_us);
  tree->Branch("Event_Id", &Event_Id);

  bool have_last = false;
  double last_hit_us = 0.0;
  ll current_event_id = 1;

  size_t out_rows = 0;

  while (!pq.empty()) {
    HeapItem top = pq.top();
    pq.pop();

    const Row &r = top.row;

    // Event grouping on globally sorted stream
    if (!have_last) {
      current_event_id = 1;
      last_hit_us = r.t_abs_us;
      have_last = true;
    } else {
      if (r.t_abs_us > last_hit_us + window_us) {
        current_event_id++;
      }
      last_hit_us = r.t_abs_us;
    }

    TStamp_us = r.TStamp_us;
    Board_Id  = r.Board_Id;
    Num_Hits  = r.Num_Hits;
    CH_Id     = r.CH_Id;
    DataType  = r.DataType;
    ToA_LSB   = r.ToA_LSB;
    ToT_LSB   = r.ToT_LSB;
    t_abs_us  = r.t_abs_us;
    Event_Id  = current_event_id;

    tree->Fill();
    out_rows++;

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


  // Optional: delete run files to save space
  for (const auto &p : run_files) {
    std::remove(p.c_str());
  }

  return 0;
}
