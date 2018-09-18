#include "delphys/analysis/interface/ResolvedAnalyser.h"

#include "TString.h"
#include "TSystem.h"

#include <algorithm>
#include <iterator>



ResolvedAnalyser::ResolvedAnalyser(const TString & in_path,
                                   const TString & out_path,
                                   const TString & out_tree_name)
    : BaseAnalyser(in_path, out_path, out_tree_name) {

  SetBranchAddress();
  MakeBranch();

  std::cout << "ctor" << std::endl;
}

ResolvedAnalyser::~ResolvedAnalyser() {
  std::cout << "dtor" << std::endl;
}

void ResolvedAnalyser::MakeBranch() {
  out_tree_->Branch("jet_p4", &jet_p4_);
  out_tree_->Branch("j0_p4", &j0_p4_);
  return ;
}


void ResolvedAnalyser::Reset() {
  jet_p4_.clear();
  j0_p4_.clear();
}





Bool_t ResolvedAnalyser::SelectEvent() {



  // Veto
  // if (VetoElectron()) return kFALSE;
  // if (VetoMuon()) return kFALSE;
  if (jets_->GetEntries() < 6) return kFALSE;

  // NOTE Select objects
  std::vector<Jet> selected_jets = SelectJet();
  if(selected_jets.size() < 6) return kFALSE;

  for (const auto & each : selected_jets) {
    jet_p4_.push_back(each.P4());
  }

  return kTRUE;
}

std::vector<Jet> ResolvedAnalyser::SelectJet() {
  std::vector<Jet> selected_jets;
  for (Int_t i = 0; i < jets_->GetEntries(); i++) {
    auto jet = dynamic_cast<const Jet*>(jets_->At(i));

    if(jet->PT < 45.0f) continue;
    if(std::fabs(jet->Eta) > 2.4f) continue;

    selected_jets.push_back(*jet);
  }
  return selected_jets;
}


void ResolvedAnalyser::AnalyseEvent() {
  out_tree_->Fill();
}

void ResolvedAnalyser::Loop() {
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
                                    100 * static_cast<Float_t>(num_passed) / kNumEntries);

      std::cout << msg;
    }

    // Bool_t is_selected = SelectEvent()
    if (not SelectEvent()) continue;
    num_passed++;


    if (num_passed == 100) break;
    AnalyseEvent();
  }
}

int main(int argc, char* argv[]) {
  TString in_path(argv[1]);
  TString out_path(argv[2]);

  Double_t dr = delphys::ComputeDeltaR(0.4, 0.2, 0.5, 0.2);
  std::cout << "DeltaR: " << dr << std::endl;

  ResolvedAnalyser analyser(in_path, out_path, "test");
  analyser.Loop();


  return 0;
}
