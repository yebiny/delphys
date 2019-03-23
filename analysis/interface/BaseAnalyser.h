#ifndef DELPHYS_ANALYSIS_BASEANALYSER_H_
#define DELPHYS_ANALYSIS_BASEANALYSER_H_

#include "delphys/external/interface/MathUtils.h"

#include "classes/DelphesClasses.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TString.h"

#include <set>

// TODO use enum??
// enum class DelphesBranches {
//   kEvent, kParticle, kGenJet, kGenMissingET,
//   kTrack, kTower, kEFlowTrack, kEFlowPhoton, kEFlowNeutralHadron,
//   kElectron, kPhoton, kMuon, kJet, kFatJet, kMissingET, kScalarHT, kVertex,
// };


class BaseAnalyser {
 public:
  BaseAnalyser();

  // constructor to clone Delphes' tree
  BaseAnalyser(const TString & in_path,
               const TString & out_path);

  BaseAnalyser(const TString & in_path,
               const TString & out_path,
               const TString & out_tree_name);

  // constructor for multiple input files
  BaseAnalyser(const std::vector<TString> & in_paths,
               const TString & out_path,
               const TString & out_tree_name);

  ~BaseAnalyser();

  void initDelphesBranch();
  void setBranchAddress(std::set<TString> branches, Bool_t drop=false);
  void setBranchAddress();

  TFile* in_file_,
  TFile* out_file_;
  TTree* in_tree_,
  TFile* out_tree_;

  TClonesArray* events_;
  TClonesArray* particles_; // GenParticle
  TClonesArray* gen_jets_;  // 
  TClonesArray* gen_mets_;
  TClonesArray* tracks_;
  TClonesArray* towers_;
  TClonesArray* eflow_tracks_;
  TClonesArray* eflow_photons_;
  TClonesArray* eflow_neutral_hadrons_;
  TClonesArray* electrons_;
  TClonesArray* muons_;
  TClonesArray* photons_;
  TClonesArray* jets_;
  TClonesArray* fat_jets_;
  TClonesArray* mets_;
  TClonesArray* scalar_hts_;
  TClonesArray* vertices_;

  const std::set<TString> kAllDelphesBranches_ = {
      "Event", "Particle", "GenJet", "GenMissingET",
      "Track", "Tower", "EFlowTrack", "EFlowPhoton", "EFlowNeutralHadron",
      "Electron", "Photon", "Muon",
      "Jet", "FatJet", "MissingET", "ScalarHT", "Vertex"
  };

};


#endif // DELPHYS_ANALYSIS_BASEANALYSER_H_
