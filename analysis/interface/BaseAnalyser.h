#ifndef DELPHYS_ANALYSIS_BASEANALYSER_H_
#define DELPHYS_ANALYSIS_BASEANALYSER_H_

#include "delphys/external/interface/MathUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TString.h"

#include "classes/DelphesClasses.h"

#include <set>



class BaseAnalyser {
 public:
  BaseAnalyser();

  BaseAnalyser(const TString & in_path,
               const TString & out_path);

  BaseAnalyser(const TString & in_path,
               const TString & out_path,
               const TString & out_tree_name);

  BaseAnalyser(const std::vector<TString> & in_paths,
               const TString & out_path,
               const TString & out_tree_name);

  ~BaseAnalyser();

  void InitDelphesBranch();
  void SetBranchAddress(std::set<TString> branches, Bool_t drop=false);
  void SetBranchAddress();

  TFile *in_file_, *out_file_;
  TTree *in_tree_, *out_tree_;

  TClonesArray *events_;
  TClonesArray *particles_; // GenParticle
  TClonesArray *gen_jets_;  // 
  TClonesArray *gen_mets_;
  TClonesArray *tracks_;
  TClonesArray *towers_;
  TClonesArray *eflow_tracks_;
  TClonesArray *eflow_photons_;
  TClonesArray *eflow_neutral_hadrons_;
  TClonesArray *electrons_;
  TClonesArray *muons_;
  TClonesArray *photons_;
  TClonesArray *jets_;
  TClonesArray *fat_jets_;
  TClonesArray *mets_;
  TClonesArray *scalar_hts_;
  TClonesArray *vertices_;

  const std::set<TString> kAllDelphesBranches = {
      "Event", "Particle", "GenJet", "GenMissingET",
      "Track", "Tower", "EFlowTrack", "EFlowPhoton", "EFlowNeutralHadron",
      "Electron", "Photon", "Muon",
      "Jet", "FatJet", "MissingET", "ScalarHT", "Vertex"};

};


#endif // DELPHYS_ANALYSIS_BASEANALYSER_H_
