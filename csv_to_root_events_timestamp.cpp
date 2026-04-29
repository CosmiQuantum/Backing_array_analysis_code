// csv_to_root_events_timestamp.cpp
//
// Reads Janus Timing CSV:
//   // comment lines...
//   TStamp_us,Board_Id,Num_Hits,CH_Id,DataType,ToA_LSB,ToT_LSB
//   0.0190,1,14,54,0x30,4231,87
//   ...
//
// Builds absolute time:
//   t_abs_us = TStamp_us + 0.0005 * ToA_LSB   (since 0.5 ns LSB = 0.0005 us)
//
// Sorts by t_abs_us and assigns Event_Id using a sliding time window.
//
// Usage:
//   ./csv_to_root_events_timestamp input.csv output.root tol_LSB [max_rows]
//
// Example:
//   ./csv_to_root_events_timestamp Run1_list.csv out.root 200

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using ll = long long;

static inline std::string trim(const std::string &s) {
  size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) b++;
  size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) e--;
  return s.substr(b, e - b);
}

static inline bool starts_with(const std::string &s, const std::string &p) {
  return s.rfind(p, 0) == 0;
}

// Split a simple CSV line into exactly 7 fields (no quoted commas expected here).
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
  ll DataType = 0;  // e.g. 0x30
  ll ToA_LSB = 0;
  ll ToT_LSB = 0;

  double t_abs_us = 0.0;
};

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr
      << "Usage:\n"
      << "  " << argv[0] << " <input.csv> <output.root> <tol_LSB> [max_rows]\n\n"
      << "Example:\n"
      << "  " << argv[0] << " Run1_list.csv Run1_list_events.root 200\n";
    return 1;
  }

  const std::string input_path  = argv[1];
  const std::string output_path = argv[2];
  const ll tol_LSB = std::stoll(argv[3]);
  const ll max_rows = (argc >= 5) ? std::stoll(argv[4]) : -1;

  // 1 LSB = 0.5 ns = 0.0005 us
  const double window_us = static_cast<double>(tol_LSB) * 0.0005;

  std::cerr << "Reading: " << input_path << "\n";
  std::cerr << "Writing: " << output_path << "\n";
  std::cerr << "tol_LSB: " << tol_LSB << "  => window_us = " << window_us << " us\n";

  std::ifstream fin(input_path);
  if (!fin.is_open()) {
    std::cerr << "ERROR: cannot open input file: " << input_path << "\n";
    return 2;
  }

  std::vector<Row> rows;
  rows.reserve(1'000'000);

  std::string line;
  ll line_no = 0;

  while (std::getline(fin, line)) {
    line_no++;
    line = trim(line);
    if (line.empty()) continue;

    // Skip Janus metadata lines starting with //
    if (starts_with(line, "//")) continue;

    // Skip ONLY the real CSV header line
    // (Do NOT skip data lines like "...0x30..." which contain letters.)
    if (line.find("TStamp_us") != std::string::npos) continue;

    // Parse CSV with 7 fields
    std::string f[7];
    if (!split_csv_7(line, f)) {
      std::cerr << "WARNING: could not parse line " << line_no << " (skipping): " << line << "\n";
      continue;
    }

    Row r;
    try {
      r.TStamp_us = std::stod(f[0]);
      r.Board_Id  = std::stoll(f[1]);
      r.Num_Hits  = std::stoll(f[2]);
      r.CH_Id     = std::stoll(f[3]);
      // base=0 handles 0x.. hex automatically
      r.DataType  = std::stoll(f[4], nullptr, 0);
      r.ToA_LSB   = std::stoll(f[5]);
      r.ToT_LSB   = std::stoll(f[6]);
    } catch (...) {
      std::cerr << "WARNING: conversion failed on line " << line_no << " (skipping): " << line << "\n";
      continue;
    }

    r.t_abs_us = r.TStamp_us + 0.0005 * static_cast<double>(r.ToA_LSB);
    rows.push_back(r);

    if (max_rows > 0 && static_cast<ll>(rows.size()) >= max_rows) {
      std::cerr << "Reached max_rows=" << max_rows << ", stopping read.\n";
      break;
    }
  }

  fin.close();

  if (rows.empty()) {
    std::cerr << "ERROR: no valid rows read from file.\n";
    return 3;
  }

  std::cerr << "Rows read: " << rows.size() << "\n";

  // Sort by absolute time
  std::sort(rows.begin(), rows.end(),
            [](const Row &a, const Row &b) {
              return a.t_abs_us < b.t_abs_us;
            });

  // ROOT output
  TFile *fout = TFile::Open(output_path.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create ROOT file: " << output_path << "\n";
    return 4;
  }

  TTree *tree = new TTree("Events", "Events grouped by absolute time");

  // Branch variables
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

  // Event grouping (sliding window):
  // New event starts when time gap > window_us from previous hit.
  ll current_event_id = 1;
  double last_hit_us = rows.front().t_abs_us;

  for (size_t i = 0; i < rows.size(); i++) {
    const Row &r = rows[i];

    if (i == 0) {
      current_event_id = 1;
      last_hit_us = r.t_abs_us;
    } else {
      if (r.t_abs_us > last_hit_us + window_us) {
        current_event_id++;
      }
      last_hit_us = r.t_abs_us;
    }

    // Fill branches
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
  }

  std::cerr << "Done. Entries: " << tree->GetEntries() << "\n";
  std::cerr << "Last Event_Id: " << current_event_id << "\n";

  fout->Write();
  fout->Close();

  return 0;
}
