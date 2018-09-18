#include "delphys/analysis/interface/TTAllJetsFilter.h"

#include "TString.h"
#include "TSystem.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <iterator>
#include <queue>


TTAllJetsFilter::TTAllJetsFilter(const TString & in_path,
                                 const TString & out_path) : BaseAnalyser() {

  std::cout << "start ctor" << std::endl;

  in_file_ = TFile::Open(in_path, "READ");
  in_tree_ = dynamic_cast<TTree*>(in_file_->Get("Delphes"));

  std::cout << "Start to set branch addressess" << std::endl;
  SetBranchAddress();
  std::cout << "Finish off setting branch addressess" << std::endl;

  out_file_ = TFile::Open(out_path, "RECREATE");
  out_tree_ = in_tree_->CloneTree(0);
  out_tree_->SetDirectory(out_file_);

  std::set<Int_t> channel_codes = {
      0,
      1, 10, 100,
      2, 11, 101, 20, 110, 200};

  for (const auto & each : channel_codes) {
    decay_channel_count_[each] = 0;
  }

  std::cout << "end ctor" << std::endl;
}


TTAllJetsFilter::~TTAllJetsFilter() {
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


void TTAllJetsFilter::MakeBranch() {

}


void TTAllJetsFilter::Reset() {

}





Bool_t TTAllJetsFilter::IsAllJetsChannel() {

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
      12, 14, 16};


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
 

  Int_t code = num_electrons + 10*num_muons + 100*num_taus;
  decay_channel_count_[code]++;

  Bool_t is_all_jets = (code == 0);

  return is_all_jets;
}



Bool_t TTAllJetsFilter::SelectEvent() {
  return IsAllJetsChannel();
}


void TTAllJetsFilter::AnalyseEvent() {
  out_tree_->Fill();
}

void TTAllJetsFilter::Loop() {
  Int_t kNumEntries = in_tree_->GetEntries();
  Int_t kPrintFreq = kNumEntries / 100;
  TString kMessageTemplate = "[%d/%d (%.2f %)] # of passed events: %d (%.2f %)\n";

  Int_t num_passed = 0;

  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);

    if (entry % kPrintFreq == 0) {
      TString msg = TString::Format(kMessageTemplate,
                                    entry + 1,
                                    kNumEntries,
                                    100 * static_cast<Float_t>(entry + 1) / kNumEntries,
                                    num_passed,
                                    100 * static_cast<Float_t>(num_passed) / (entry + 1));

      std::cout << msg;
    } 

    if (SelectEvent()) {
      AnalyseEvent();
      num_passed++;
    }
  }
}


int main(int argc, char* argv[]) {
  TString in_path(argv[1]);
  TString out_path(argv[2]);

  TTAllJetsFilter analyser(in_path, out_path);
  analyser.Loop();

  return 0;
}
