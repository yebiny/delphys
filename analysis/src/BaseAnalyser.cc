#include "delphys/analysis/interface/BaseAnalyser.h"

BaseAnalyser::BaseAnalyser() {
  
}

BaseAnalyser::BaseAnalyser(const TString & in_path,
                           const TString & out_path) {
  in_file_ = TFile::Open(in_path, "READ");
  in_tree_ = dynamic_cast<TTree*>(in_file_->Get("Delphes"));

  out_file_ = TFile::Open(out_path, "RECREATE");
  out_tree_ = in_tree_->CloneTree(0);
  out_tree_->SetDirectory(out_file_);
}



BaseAnalyser::BaseAnalyser(const TString & in_path,
                           const TString & out_path,
                           const TString & out_tree_name) {
  in_file_ = TFile::Open(in_path, "READ");
  in_tree_ = dynamic_cast<TTree*>(in_file_->Get("Delphes"));

  out_file_ = TFile::Open(out_path, "RECREATE");
  out_tree_ = new TTree(out_tree_name, out_tree_name);
  out_tree_->SetDirectory(out_file_);
}


BaseAnalyser::BaseAnalyser(const std::vector<TString> & in_paths,
                           const TString & out_path,
                           const TString & out_tree_name) {
  TChain chain("Delphes");
  for (const auto & each : in_paths) {
    chain.Add(each);
  }

  in_file_ = nullptr;
  in_tree_ = dynamic_cast<TTree*>(&chain);

  out_file_ = TFile::Open(out_path, "RECREATE");
  out_tree_ = new TTree(out_tree_name, out_tree_name);
  out_tree_->SetDirectory(out_file_);
}

BaseAnalyser::~BaseAnalyser() {
  if ((in_file_ != nullptr) and (in_file_->IsOpen())) {
    in_file_->Close();
  }

  if (out_file_->IsOpen()) {
    out_file_->Write();
    out_file_->Close();
  }
}


void BaseAnalyser::InitDelphesBranch() {
  events_ = nullptr;
  particles_ = nullptr;
  gen_jets_ = nullptr;
  gen_mets_ = nullptr;
  tracks_ = nullptr;
  towers_ = nullptr;
  eflow_tracks_ = nullptr;
  eflow_photons_ = nullptr;
  eflow_neutral_hadrons_ = nullptr;
  electrons_ = nullptr;
  muons_ = nullptr;
  photons_ = nullptr;
  jets_ = nullptr;
  fat_jets_ = nullptr;
  mets_ = nullptr;
  scalar_hts_ = nullptr;
  vertices_ = nullptr;
}


void BaseAnalyser::SetBranchAddress(std::set<TString> branches, Bool_t drop) {
  InitDelphesBranch();

  std::set<TString> in_branches;

  if (drop) {
    std::set_difference(kAllDelphesBranches.begin(), kAllDelphesBranches.end(), 
                        branches.begin(), branches.end(),
                        std::inserter(in_branches, in_branches.begin()));
  } else {
    in_branches = std::move(branches);
  }


  for(const auto & each : in_branches) {
    if (each.EqualTo("Event"))
      in_tree_->SetBranchAddress("Event", &events_);
    else if (each.EqualTo("Particle"))
      in_tree_->SetBranchAddress("Particle", &particles_);
    else if (each.EqualTo("GenJet"))
      in_tree_->SetBranchAddress("GenJet", &gen_jets_);
    else if (each.EqualTo("GenMissingET"))
      in_tree_->SetBranchAddress("GenMissingET", &gen_mets_);
    else if (each.EqualTo("Track"))
      in_tree_->SetBranchAddress("Track", &tracks_);
    else if (each.EqualTo("Tower"))
      in_tree_->SetBranchAddress("Tower", &towers_);
    else if (each.EqualTo("EFlowTrack"))
      in_tree_->SetBranchAddress("EFlowTrack", &eflow_tracks_);
    else if (each.EqualTo("EFlowPhoton"))
      in_tree_->SetBranchAddress("EFlowPhoton", &eflow_photons_);
    else if (each.EqualTo("EFlowNeutralHadron"))
      in_tree_->SetBranchAddress("EFlowNeutralHadron", &eflow_neutral_hadrons_);
    else if (each.EqualTo("Electron"))
      in_tree_->SetBranchAddress("Electron", &electrons_);
    else if (each.EqualTo("Photon"))
      in_tree_->SetBranchAddress("Photon", &photons_);
    else if (each.EqualTo("Muon"))
      in_tree_->SetBranchAddress("Muon", &muons_);
    else if (each.EqualTo("Jet"))
      in_tree_->SetBranchAddress("Jet",      &jets_);
    else if (each.EqualTo("FatJet"))
      in_tree_->SetBranchAddress("FatJet", &fat_jets_);
    else if (each.EqualTo("MissingET"))
      in_tree_->SetBranchAddress("MissingET", &mets_);
    else if (each.EqualTo("ScalarHT"))
      in_tree_->SetBranchAddress("ScalarHT", &scalar_hts_);
    else if (each.EqualTo("Vertex"))
      in_tree_->SetBranchAddress("Vertex", &vertices_);
    else
      std::cout << each << std::endl;
  }
}


void BaseAnalyser::SetBranchAddress() {
  SetBranchAddress(kAllDelphesBranches);
}
