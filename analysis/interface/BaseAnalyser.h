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
               const TString & out_path,
               const TString & out_tree_name="delphys");

  BaseAnalyser(const std::vector<TString> & in_paths,
               const TString & out_path,
               const TString & out_tree_name="delphys");

  ~BaseAnalyser();

  void InitDelphesBranch();
  void SetBranchAddress(std::set<TString> branches, Bool_t keep=true);
  void SetBranchAddress();

  virtual void MakeBranch() = 0;
  virtual void Reset() = 0;
  virtual Bool_t SelectEvent() = 0;
  virtual void AnalyseEvent() = 0;
  virtual void Loop() = 0;

  TFile *in_file_, *out_file_;
  TTree *in_tree_, *out_tree_;

  TClonesArray *events_, *particles_, *gen_jets_, *gen_mets_;
  TClonesArray *tracks_, *towers_;
  TClonesArray *eflow_tracks_, *eflow_photons_, *eflow_nhads_;
  TClonesArray *electrons_, *muons_, *photons_;
  TClonesArray *jets_, *fat_jets_, *mets_, *scalar_hts_, *vertices_;

  const std::set<TString> kAllDelphesBranches = {
      "Event", "Particle", "GenJet", "GenMissingET",
      "Track", "Tower", "EFlowTrack", "EFlowPhoton", "EFlowNeutralHadron",
      "Electron", "Photon", "Muon",
      "Jet", "FatJet", "MissingET", "ScalarHT", "Vertex"};

};


#endif // DELPHYS_ANALYSIS_BASEANALYSER_H_
