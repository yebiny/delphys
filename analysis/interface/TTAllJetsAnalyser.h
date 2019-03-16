#ifndef DELPHYS_ANALYSIS_TTALLJETSANALYSER_H_
#define DELPHYS_ANALYSIS_TTALLJETSANALYSER_H_

#include "delphys/analysis/interface/BaseAnalyser.h"

#include <numeric>


class TTAllJetsAnalyser : private BaseAnalyser {
 public:
  TTAllJetsAnalyser(const TString & in_path,
                   const TString & out_path,
                   const TString & out_tree_name);
  ~TTAllJetsAnalyser();
  void Loop();

 private:
  // inherited
  void MakeBranch();
  void ResetBranch();
  Bool_t SelectEvent();
  void Analyse();

  std::vector<const Jet*> SelectJet();

  void FillEFlow();
  void FillJetVariables();

  Bool_t TrackBottomQuark(const GenParticle* p);
  Float_t GetBDaughterRatio(const Jet* jet);



  std::vector<const Jet*> selected_jets_;



  // per event
  Int_t label_;

  // unordered set
  std::vector<Float_t> eflow_pt_;
  std::vector<Float_t> eflow_eta_;
  std::vector<Float_t> eflow_phi_;
  std::vector<Int_t>   eflow_charge_;
  std::vector<Int_t>   eflow_pid_; // PDG id
  std::vector<Int_t>   eflow_type_; //

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


  const Int_t kElectronPID_ = 11;
  const Int_t kMuonPID_ = 13;
  const Int_t kWrongPID_ = std::numeric_limits<Int_t>::max();

  const Int_t kBMatchingDeltaR_ = 0.3;

};



#endif //  DELPHYS_ANALYSIS_TTALLJETSANALYSER_H_
