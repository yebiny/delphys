#ifndef DELPHYS_ANALYSIS_TTALLJETSANALYSER_H_
#define DELPHYS_ANALYSIS_TTALLJETSANALYSER_H_

#include "delphys/analysis/interface/BaseAnalyser.h"

#include <numeric>

namespace pdgid {
  static const Int_t kBottom = 5;
  static const Int_t kTop = 6;
  static const Int_t kElectron = 11;
  static const Int_t kMuon = 13;
  static const Int_t kPhoton = 22;
  static const Int_t kWBoson = 24;

  static const Int_t kAntiElectron = -11;
  static const Int_t kAntiMuon = -13;

  // NOTE
  static const Int_t kWrong = std::numeric_limits<Int_t>::max();
  // In TDatabasePDG, 0 means Rootino, which indicates unidentified particle.
  static const Int_t kNeutralHadron = 0;

} // pdgid

// indices for ParticleFlow objects
// consider Make json file and then read from that file
// In deep learning framework like Keras, embedding layaer takes these indices
// as input.
namespace pfindex {
  static const Int_t kElectron = 1;
  static const Int_t kAntiElectron = 2;
  static const Int_t kMuon = 3;
  static const Int_t kAntiMuon = 4;
  static const Int_t kPositivelyChargedHadron = 5;
  static const Int_t kNegativelyChargedHadron = 6;
  static const Int_t kNeutralHadron = 7; 
  static const Int_t kPhoton = 8;
} //pfindex

// FIXME write the comment
namespace p8status {
  static const Int_t kBottom = 23;
} // p8status


class TTAllJetsAnalyser : private BaseAnalyser {
 public:
  TTAllJetsAnalyser(const TString & in_path,
                    const TString & out_path,
                    Int_t label);

  ~TTAllJetsAnalyser();
  void Loop();

 private:
  // inherited
  void makeBranch();
  void resetBranch();
  Bool_t selectEvent();
  void analyse();

  void fillEFlow();
  void fillJetVariables();

  Bool_t trackBottomQuark(const GenParticle* p);
  Float_t getBDaughterRatio(const Jet* jet);

  Int_t getPFIndex(Int_t pid, Int_t charge);

  std::vector<const Jet*> selected_jets_;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Branches
  //////////////////////////////////////////////////////////////////////////////
  // per event
  Int_t label_;

  // unordered set
  Int_t num_eflow_;

  std::vector<Float_t> eflow_pt_;
  std::vector<Float_t> eflow_eta_;
  std::vector<Float_t> eflow_phi_;
  std::vector<Int_t>   eflow_charge_;
  std::vector<Int_t>   eflow_pid_; // PDG id
  std::vector<Int_t>   eflow_type_; //

  Float_t met_;
  Float_t met_eta_;
  Float_t met_phi_;

  // 
  std::vector<Float_t> jet_pt_;
  std::vector<Float_t> jet_eta_;
  std::vector<Float_t> jet_phi_;
  std::vector<Float_t> jet_pid_;

  std::vector<Int_t> jet_num_chad_;
  std::vector<Int_t> jet_num_nhad_;
  std::vector<Int_t> jet_num_electron_;
  std::vector<Int_t> jet_num_muon_;
  std::vector<Int_t> jet_num_photon_;

  std::vector<Float_t> jet_major_axis_;
  std::vector<Float_t> jet_minor_axis_;
  std::vector<Float_t> jet_ptd_;

  std::vector<Bool_t> jet_b_tag_;
  std::vector<Bool_t> jet_b_dr_matching_;
  std::vector<Bool_t> jet_b_tracking_;

  std::vector<std::vector<Float_t> > con_pt_;
  std::vector<std::vector<Float_t> > con_deta_;
  std::vector<std::vector<Float_t> > con_dphi_;
  std::vector<std::vector<Int_t> > con_charge_;
  std::vector<std::vector<Int_t> > con_pid_;
  std::vector<std::vector<Int_t> > con_type_;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE selection cut
  //////////////////////////////////////////////////////////////////////////////
  const Float_t kBMatchingDeltaR_ = 0.3;
  const Float_t kMinJetPT_ = 45.0f;
  const Float_t kMaxJetEta_ = 2.4f;
  const Float_t kMinJetMultiplicity_ = 6;

};



#endif //  DELPHYS_ANALYSIS_TTALLJETSANALYSER_H_
