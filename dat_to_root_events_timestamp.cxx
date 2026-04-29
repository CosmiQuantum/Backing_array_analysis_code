#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <cmath>
#include <string>
#include <cstring>

using ll = long long;

// ---------------- Ring lookup (same as before) ----------------
static std::unordered_map<unsigned long long, int> build_ring_lookup() {
  std::unordered_map<unsigned long long, int> ring_lookup;
  ring_lookup.reserve(512);

  auto key = [](ll board, ll ch) -> unsigned long long {
    return (static_cast<unsigned long long>(static_cast<unsigned int>(board)) << 32) |
           (static_cast<unsigned long long>(static_cast<unsigned int>(ch)));
  };

  const int ring_counter_assignment[28] = {
      1, 2, 15, 16, 3, 4, 5, 6,
      17, 18, 19, 20, 7, 8, 9, 10,
      21, 22, 23, 24, 11, 12, 13, 14,
      25, 26, 27, 28};

  std::vector<std::vector<std::vector<int>>> ring_map = {
      {{47,45,43,41,46,44,42,40},{63,61,59,57,62,60,58,56},{39,37,35,33,38,36,34,32},{55,53,51,49,54,52,50,48}},
      {{15,14,13,12,11,10,9,8},{31,30,29,28,27,26,25,24},{47,46,45,44,43,42,41,40},{63,62,61,60,59,58,57,56},
       {7,6,5,4,3,2,1,0},{23,22,21,20,19,18,17,16},{39,38,37,36,35,34,33,32},{55,54,53,52,51,50,49,48}},
      {{15,14,13,12,11,10,9,8},{31,30,29,28,27,26,25,24},{47,46,45,44,43,42,41,40},{63,62,61,60,59,58,57,56},
       {7,6,5,4,3,2,1,0},{23,22,21,20,19,18,17,16},{39,38,37,36,35,34,33,32},{55,54,53,52,51,50,49,48}},
      {{15,14,13,12,11,10,9,8},{31,30,29,28,27,26,25,24},{47,46,45,44,43,42,41,40},{63,62,61,60,59,58,57,56},
       {7,6,5,4,3,2,1,0},{23,22,21,20,19,18,17,16},{39,38,37,36,35,34,33,32},{55,54,53,52,51,50,49,48}},
  };

  int ring_counter = 0;
  for (int board = 0; board < (int)ring_map.size(); board++) {
    for (auto &ring : ring_map[board]) {
      int ring_id = ring_counter_assignment[ring_counter];
      for (int ch : ring) ring_lookup[key(board, ch)] = ring_id;
      ring_counter++;
    }
  }
  return ring_lookup;
}

// ---------------- Little-endian reads ----------------
static bool read_u32(std::ifstream& in, uint32_t& v) { return (bool)in.read((char*)&v, 4); }
static bool read_f64(std::ifstream& in, double& v)   { return (bool)in.read((char*)&v, 8); }

// ---------------- Hit decode ----------------
// packed = (dtype<<24) | (ch<<16) | tot
static inline uint8_t  dtype(uint32_t w) { return uint8_t((w >> 24) & 0xFF); }
static inline uint8_t  ch(uint32_t w)    { return uint8_t((w >> 16) & 0xFF); }
static inline uint16_t tot(uint32_t w)   { return uint16_t(w & 0xFFFF); }

static inline bool hit_word_ok(uint32_t w) {
  if (dtype(w) != 0x30) return false;
  if (ch(w) > 63) return false;
  if (tot(w) > 5000) return false; // loose
  return true;
}

static inline bool ts_ok(double x) {
  // Your CSV has ~67.263; allow broad but reject denorm garbage
  return std::isfinite(x) && x > 1e-9 && x < 1e12;
}

// ============================================================
// Block sync heuristic:
// We scan forward and try to interpret bytes as:
//
//   double   TStamp_us
//   uint32   Board_Id
//   uint32   Num_Hits
//
// Then verify next 2*Num_Hits words are valid hit pairs.
// If verified -> we're aligned.
// ============================================================
static bool try_sync_at(std::ifstream& in, std::streampos pos,
                        double& out_ts, uint32_t& out_board, uint32_t& out_nhits) {
  in.clear();
  in.seekg(pos);

  double ts;
  uint32_t board, nhits;
  if (!read_f64(in, ts)) return false;
  if (!read_u32(in, board)) return false;
  if (!read_u32(in, nhits)) return false;

  if (!ts_ok(ts)) return false;
  if (nhits == 0 || nhits > 500000) return false; // super loose upper bound
  if (board > 1000) return false; // board index should be small

  // Verify following hit pairs
  std::streampos after_hdr = in.tellg();

  for (uint32_t i = 0; i < nhits; i++) {
    uint32_t w1, w2;
    if (!read_u32(in, w1)) return false;
    if (!read_u32(in, w2)) return false; // ToA
    if (!hit_word_ok(w1)) return false;
    // ToA can be any 32-bit positive; no strict check
    (void)w2;
  }

  // success
  out_ts = ts; out_board = board; out_nhits = nhits;

  // restore stream to right after header (so caller can parse)
  in.clear();
  in.seekg(after_hdr);
  return true;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage:\n  " << argv[0] << " input.dat output.root [tol]\n";
    return 1;
  }
  const std::string in_dat = argv[1];
  const std::string out_root = argv[2];
  const ll tol = (argc >= 4) ? std::stoll(argv[3]) : 200;

  std::ifstream in(in_dat, std::ios::binary);
  if (!in) { std::cerr << "ERROR open\n"; return 2; }

  auto ring_lookup = build_ring_lookup();
  auto rkey = [](ll board, ll ch) -> unsigned long long {
    return (static_cast<unsigned long long>(static_cast<unsigned int>(board)) << 32) |
           (static_cast<unsigned long long>(static_cast<unsigned int>(ch)));
  };

  // ROOT output
  TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
  TTree* tree = new TTree("Events", "Events from DAT with Event_Id + Ring");

  double TStamp_us = 0.0;
  ll Board_Id = 0, CH_Id=0, ToA_LSB=0, ToT_LSB=0, Event_Id=0;
  int Ring=-1;

  tree->Branch("TStamp_us",&TStamp_us);
  tree->Branch("Board_Id",&Board_Id);
  tree->Branch("CH_Id",&CH_Id);
  tree->Branch("ToA_LSB",&ToA_LSB);
  tree->Branch("ToT_LSB",&ToT_LSB);
  tree->Branch("Ring",&Ring);
  tree->Branch("Event_Id",&Event_Id);

  // Event assignment (same as your python)
  ll current_event_id = 0;
  std::vector<ll> current_event_toas;
  current_event_toas.reserve(256);

  auto belongs = [&](ll t)->bool{
    for (auto prev: current_event_toas) if (std::llabs(t-prev) <= tol) return true;
    return false;
  };

  // ---------------- SYNC ----------------
  in.seekg(0, std::ios::end);
  std::streampos fsize = in.tellg();
  in.seekg(0);

  // Scan file for a valid block header alignment (step by 4 bytes)
  double ts0=0.0;
  uint32_t b0=0, nh0=0;
  bool synced=false;
  std::streampos sync_pos = 0;

  for (std::streampos p = 0; p + std::streamoff(16) < fsize; p += std::streamoff(4)) {
    if (try_sync_at(in, p, ts0, b0, nh0)) {
      synced = true;
      sync_pos = p;
      break;
    }
  }

  if (!synced) {
    std::cerr << "ERROR: could not sync to DAT block header.\n";
    return 3;
  }

  std::cout << "Synced at byte offset " << (long long)sync_pos
            << "  TStamp_us=" << ts0
            << "  Board=" << b0
            << "  Num_Hits=" << nh0 << "\n";

  // After try_sync_at, stream is positioned right AFTER header of the first block
  // But we still need to remember ts0/b0/nh0 as current block header values.
  TStamp_us = ts0;
  Board_Id = (ll)b0;

  bool first_hit_seen = false;
  long long n_written = 0;

  while (true) {
    // Parse nh0 hit pairs
    for (uint32_t i = 0; i < nh0; i++) {
      uint32_t w1, w2;
      if (!read_u32(in, w1)) goto done;
      if (!read_u32(in, w2)) goto done;

      if (!hit_word_ok(w1)) {
        std::cerr << "ERROR: lost sync inside block (bad hit word)\n";
        goto done;
      }

      CH_Id = (ll)ch(w1);
      ToT_LSB = (ll)tot(w1);
      ToA_LSB = (ll)w2;

      auto it = ring_lookup.find(rkey(Board_Id, CH_Id));
      Ring = (it != ring_lookup.end()) ? it->second : -1;

      if (!first_hit_seen) {
        current_event_id = 1;
        current_event_toas.clear();
        current_event_toas.push_back(ToA_LSB);
        first_hit_seen = true;
      } else if (belongs(ToA_LSB)) {
        current_event_toas.push_back(ToA_LSB);
      } else {
        current_event_id++;
        current_event_toas.clear();
        current_event_toas.push_back(ToA_LSB);
      }
      Event_Id = current_event_id;

      tree->Fill();
      n_written++;
    }

    // Read next block header: (double ts, u32 board, u32 nhits)
    double ts;
    uint32_t board, nh;
    if (!read_f64(in, ts)) break;
    if (!read_u32(in, board)) break;
    if (!read_u32(in, nh)) break;

    if (!ts_ok(ts) || nh == 0 || board > 1000) {
      std::cerr << "ERROR: invalid next block header (lost sync)\n";
      break;
    }

    TStamp_us = ts;
    Board_Id = (ll)board;
    nh0 = nh;
  }

done:
  fout->cd();
  tree->Write();
  fout->Close();

  std::cout << "Wrote " << n_written << " hits to " << out_root << "\n";
  return 0;
}
