#include "delphys/analysis/interface/TTAllJetsSelector.h"

#include "TString.h"
#include "TSystem.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <iterator>
#include <queue>


TTAllJetsSelector::TTAllJetsSelector(const TString & in_path,
                                     const TString & out_path) :
    BaseAnalyser(in_path, out_path) {
  std::cout << "begin ctor" << std::endl;

  SetBranchAddress();

  std::set<Int_t> channel_codes = {
      0,
      1, 10, 100,
      2, 11, 101, 20, 110, 200
  };

  for (const auto & each : channel_codes) {
    decay_channel_count_[each] = 0;
  }

  std::cout << "end ctor" << std::endl;
}


TTAllJetsSelector::~TTAllJetsSelector() {
  std::cout << "dtor" << std::endl;

  Float_t num_events = static_cast<Float_t>(in_tree_->GetEntries());

  for (const auto & each : decay_channel_count_) {
    switch (each.first) {
      case 0:
        std::cout << "All jets";
        break;
      case 1:
        std::cout << "electron + jets";
        break;
      case 10:
        std::cout << "muon + jets";
        break;
      case 100:
        std::cout << "tau + jets";
        break;
      case 2:
        std::cout << "e + e";
        break;
      case 11:
        std::cout << "e + mu";
        break;
      case 101:
        std::cout << "e + tau";
        break;
      case 20:
        std::cout << "mu + mu";
        break;
      case 110:
        std::cout << "mu + tau";
        break;
      case 200:
        std::cout << "tau + tau";
        break;
      default:
        std::cout << ":(";
    }

    std::cout << ": " << each.second / num_events << std::endl;
  }

}


Bool_t TTAllJetsSelector::IsAllJetsChannel() {

  // Find top quarks
  std::vector<GenParticle*> top_quarks;
  for (Int_t idx = 0; idx <= particles_->GetEntries(); idx++) {
    auto p = dynamic_cast<GenParticle*>(particles_->At(idx));
    if ((std::abs(p->PID) == 6) and (p->Status == 62)) {
      top_quarks.push_back(p);
    }
    // NOTE Assume
    if (top_quarks.size() == 2) break;
  }

  std::queue<Int_t> daughter_indices;
  for (const auto & top : top_quarks) {
    // top daughters
    for (Int_t dau_idx = top->D1; dau_idx <= top->D2; dau_idx++) {
      auto dau = dynamic_cast<GenParticle*>(particles_->At(dau_idx));
      Int_t pid = std::abs(dau->PID);

      if (pid == 24) {
        for (Int_t i = dau->D1; i <= dau->D2; i++) {
          daughter_indices.push(i);
        }
      } else if (pid != 5) {
        std::cout << top->PID << " --> " << dau->PID << std::endl;
      } 
    }
  }


  static const std::set<Int_t> kPassPID = {
      1, 2, 3, 4, 5,
      12, 14, 16
  };


  Int_t num_electrons = 0, num_muons = 0, num_taus = 0;

  while (daughter_indices.size()) {
    Int_t idx = daughter_indices.front();
    daughter_indices.pop();
    auto p = dynamic_cast<GenParticle*>(particles_->At(idx));
    Int_t pid = std::abs(p->PID);

    if (pid == 11) {
      num_electrons++;
    } else if (pid == 13) {
      num_muons++;
    } else if (pid == 15) {
      num_taus++;
    } else if (std::find(kPassPID.begin(), kPassPID.end(), pid) != kPassPID.end()) {
      continue;
    } else {
      for (Int_t dau_idx = p->D1; dau_idx <= p->D2; dau_idx++) {
        daughter_indices.push(dau_idx);
      }
    }
  } // end while loop
 

  Int_t code = num_electrons + 10 * num_muons + 100 * num_taus;
  decay_channel_count_[code]++;

  Bool_t is_all_jets = (code == 0);

  return is_all_jets;
}



Bool_t TTAllJetsSelector::SelectEvent() {
  return IsAllJetsChannel();
}


void TTAllJetsSelector::Analyse() {
  out_tree_->Fill();
}

void TTAllJetsSelector::Loop() {
  const Int_t kNumEntries = in_tree_->GetEntries();
  const Int_t kPrintFreq = kNumEntries / 100;
  const TString kMsgFmt = "[%d/%d (%.2f %)] # of passed events: %d (%.2f %)";

  Int_t num_passed = 0;

  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);

    if (entry % kPrintFreq == 0) {
      Float_t progress = 100 * static_cast<Float_t>(entry + 1) / kNumEntries;
      Float_t eff = 100 * static_cast<Float_t>(num_passed) / (entry + 1);

      auto msg = TString::Format(kMsgFmt, entry + 1, kNumEntries, progress,
                                 num_passed, eff);
      std::cout << msg << std::endl;
    } 

    if (SelectEvent()) {
      Analyse();
      num_passed++;
    }
  }
}



