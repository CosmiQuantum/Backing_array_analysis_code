#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdint>
#include <string>

using ll = long long;

// Exact channel key within an event-ring
struct ChanKey {
  ll board;
  ll ch;

  bool operator==(const ChanKey& o) const {
    return board == o.board && ch == o.ch;
  }
};

static inline std::size_t hash_combine(std::size_t seed, std::size_t v) {
  seed ^= v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
  return seed;
}

struct ChanKeyHash {
  std::size_t operator()(const ChanKey& k) const {
    std::size_t seed = 0;
    seed = hash_combine(seed, std::hash<ll>{}(k.board));
    seed = hash_combine(seed, std::hash<ll>{}(k.ch));
    return seed;
  }
};

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage:\n"
              << "  " << argv[0] << " input.root output.root [treeName] [tot_enable]\n\n"
              << "Arguments:\n"
              << "  input.root   : input ROOT file\n"
              << "  output.root  : output ROOT file\n"
              << "  treeName     : input TTree name, default = Events\n"
              << "  tot_enable   : 0 = ignore ToT_LSB cut, 1 = require ToT_LSB > 0\n\n"
              << "Examples:\n"
              << "  " << argv[0] << " input.root hist.root Events 0\n"
              << "  " << argv[0] << " input.root hist.root Events 1\n";
    return 1;
  }

  const std::string inpath  = argv[1];
  const std::string outpath = argv[2];
  const std::string tname   = (argc >= 4) ? argv[3] : "Events";
  const bool tot_enable     = (argc >= 5) ? (std::stoi(argv[4]) != 0) : false;

  std::cout << "Input file:  " << inpath << "\n";
  std::cout << "Output file: " << outpath << "\n";
  std::cout << "Tree name:   " << tname << "\n";
  std::cout << "ToT filter:  " << (tot_enable ? "ENABLED, require ToT_LSB > 0" : "DISABLED") << "\n";

  TFile* fin = TFile::Open(inpath.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open input file: " << inpath << "\n";
    return 2;
  }

  TTree* tr = dynamic_cast<TTree*>(fin->Get(tname.c_str()));
  if (!tr) {
    std::cerr << "ERROR: cannot find tree '" << tname << "'\n";
    fin->Close();
    return 3;
  }

  ll Event_Id = 0;
  ll Board_Id = 0;
  ll CH_Id = 0;
  ll ToT_LSB = 0;
  int Ring = -1;

  tr->SetBranchAddress("Event_Id", &Event_Id);
  tr->SetBranchAddress("Board_Id", &Board_Id);
  tr->SetBranchAddress("CH_Id", &CH_Id);
  tr->SetBranchAddress("Ring", &Ring);

  if (tot_enable) {
    if (!tr->GetBranch("ToT_LSB")) {
      std::cerr << "ERROR: tot_enable=1 but branch ToT_LSB does not exist.\n";
      fin->Close();
      return 4;
    }
    tr->SetBranchAddress("ToT_LSB", &ToT_LSB);
  }

  // Output accumulator:
  // per ring counts of events with >= threshold unique channels
  // thresholds: 1 through 8
  std::unordered_map<int, std::array<long long, 8>> hist_by_ring;
  for (int r = 1; r <= 28; r++) {
    hist_by_ring[r] = {0, 0, 0, 0, 0, 0, 0, 0};
  }

  // For current event:
  // ring -> set of unique (Board_Id, CH_Id)
  std::unordered_map<int, std::unordered_set<ChanKey, ChanKeyHash>> uniq_by_ring;

  auto finalize_event = [&]() {
    for (auto& kv : uniq_by_ring) {
      int ring = kv.first;
      int cnt = static_cast<int>(kv.second.size());

      if (hist_by_ring.find(ring) == hist_by_ring.end()) {
        hist_by_ring[ring] = {0, 0, 0, 0, 0, 0, 0, 0};
      }

      for (int thr = 1; thr <= 8; thr++) {
        if (cnt >= thr) {
          hist_by_ring[ring][thr - 1] += 1;
        }
      }
    }

    uniq_by_ring.clear();
  };

  const ll n = tr->GetEntries();

  long long kept_hits = 0;
  long long skipped_bad_ring = 0;
  long long skipped_bad_event = 0;
  long long skipped_tot = 0;

  if (n == 0) {
    std::cerr << "WARNING: input tree is empty.\n";
  } else {
    ll current_event = -1;

    for (ll i = 0; i < n; i++) {
      tr->GetEntry(i);

      if (Event_Id <= 0) {
        skipped_bad_event++;
        continue;
      }

      if (Ring < 0) {
        skipped_bad_ring++;
        continue;
      }

      if (tot_enable && ToT_LSB <= 0) {
        skipped_tot++;
        continue;
      }

      if (current_event < 0) {
        current_event = Event_Id;
      }

      if (Event_Id != current_event) {
        finalize_event();
        current_event = Event_Id;
      }

      auto& s = uniq_by_ring[Ring];
      s.insert(ChanKey{Board_Id, CH_Id});
      kept_hits++;
    }

    finalize_event();
  }

  TFile* fout = TFile::Open(outpath.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: cannot create output file: " << outpath << "\n";
    fin->Close();
    return 5;
  }

  TTree* th = new TTree(
      "hist_df",
      "hist_df: per ring event counts for >= threshold unique channels"
  );

  int out_ring = 0;
  long long ge_1 = 0;
  long long ge_2 = 0;
  long long ge_3 = 0;
  long long ge_4 = 0;
  long long ge_5 = 0;
  long long ge_6 = 0;
  long long ge_7 = 0;
  long long ge_8 = 0;

  th->Branch("Ring", &out_ring);
  th->Branch("ge_1", &ge_1);
  th->Branch("ge_2", &ge_2);
  th->Branch("ge_3", &ge_3);
  th->Branch("ge_4", &ge_4);
  th->Branch("ge_5", &ge_5);
  th->Branch("ge_6", &ge_6);
  th->Branch("ge_7", &ge_7);
  th->Branch("ge_8", &ge_8);

  std::vector<int> rings;
  rings.reserve(hist_by_ring.size());

  for (auto& kv : hist_by_ring) {
    rings.push_back(kv.first);
  }

  std::sort(rings.begin(), rings.end());

  for (int r : rings) {
    out_ring = r;

    const auto& a = hist_by_ring[r];

    ge_1 = a[0];
    ge_2 = a[1];
    ge_3 = a[2];
    ge_4 = a[3];
    ge_5 = a[4];
    ge_6 = a[5];
    ge_7 = a[6];
    ge_8 = a[7];

    th->Fill();
  }

  fout->cd();
  th->Write();
  fout->Close();
  fin->Close();

  std::cout << "\nSummary\n";
  std::cout << "-------\n";
  std::cout << "Input entries:        " << n << "\n";
  std::cout << "Kept hits:            " << kept_hits << "\n";
  std::cout << "Skipped bad Event_Id: " << skipped_bad_event << "\n";
  std::cout << "Skipped bad Ring:     " << skipped_bad_ring << "\n";
  std::cout << "Skipped by ToT cut:   " << skipped_tot << "\n";
  std::cout << "Wrote output:         " << outpath << "\n";

  return 0;
}
