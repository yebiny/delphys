#ifndef DELPHYS_ANALYSIS_QGJETALYSER_H_
#define DELPHYS_ANALYSIS_QGJETALYSER_H_

#include "delphys/analysis/interface/BaseAnalyser.h"


class QGJetsAnalyser : private BaseAnalyser {
 public:
  QGJetsAnalyser(const TString & in_path,
                 const TString & out_path,
                 const TString & out_tree_name,
                 Bool_t is_dijet,
                 Int_t label);
  ~QGJetsAnalyser();
  void Loop();

 private:
  // inherited
  void makeBranch();
  Bool_t selectEvent();
  void analyse(Int_t entry);
  void resetOnEachJet();
  void resetOnEachEvent();

  Bool_t passZjets(TClonesArray* jets,
                   TClonesArray* muons,
                   TClonesArray* electrons);

  Bool_t isBalanced(TClonesArray* jets);

  void fillDaughters(const Jet* jet);

  Int_t pixelate(Float_t deta, Float_t dphi,
                 Int_t num_deta_bins=33, Int_t num_dphi_bins=33);

  void makeJetImage();

  // Bool_t extractSatellites();

  //////////////////////////////////////////////////////////////////////////////
  // Branches of the output tree
  //////////////////////////////////////////////////////////////////////////////
  Int_t event_;
  Int_t num_jets_;
  Int_t num_good_jets_;
  Int_t num_primary_vertices_;
  Int_t order_;

  Int_t label_;

  Float_t pt_; // jet pt
  Float_t eta_; // jet eta
  Float_t phi_; // jet phi
  Float_t pt_dr_log_;
  Float_t ptd_; // jet energy sharing variable
  Float_t major_axis_;
  Float_t minor_axis_;

  Int_t chad_mult_;
  Int_t nhad_mult_;
  Int_t electron_mult_;
  Int_t muon_mult_;
  Int_t photon_mult_;

  Int_t parton_id_;
  Int_t flavor_id_;
  Int_t flavor_algo_id_;
  Int_t flavor_phys_id_;

  Int_t num_dau_;
  Bool_t matched_;
  // Bool_t balanced;
  // Bool_t pass_zjets;
  Bool_t lepton_overlap_; //

  std::vector<TLorentzVector> dau_p4_;
  std::vector<Float_t>        dau_pt_;
  std::vector<Float_t>        dau_deta_;
  std::vector<Float_t>        dau_dphi_;
  std::vector<Int_t>          dau_charge_;
  std::vector<Int_t>          dau_pid_;
  std::vector<Bool_t>         dau_is_hadronic_;
  std::vector<Bool_t>         dau_is_track_;
  std::vector<Float_t>        dau_eemfrac_;
  std::vector<Float_t>        dau_ehadfrac_;

  // 33*33 = 1089
  Float_t image_chad_pt_33_[1089];
  Float_t image_nhad_pt_33_[1089];
  Float_t image_electron_pt_33_[1089];
  Float_t image_muon_pt_33_[1089];
  Float_t image_photon_pt_33_[1089];
  // multiplicity
  Float_t image_chad_mult_33_[1089];
  Float_t image_nhad_mult_33_[1089];
  Float_t image_electron_mult_33_[1089];
  Float_t image_muon_mult_33_[1089];
  Float_t image_photon_mult_33_[1089];

  //////////////////////////////////////////////////////////////////////////////
  // Internal
  //////////////////////////////////////////////////////////////////////////////
  Bool_t bad_hard_gen_seen_;
  std::map<TString, Int_t> qgjets_stats_;

  //////////////////////////////////////////////////////////////////////////////
  // Constants
  //////////////////////////////////////////////////////////////////////////////
  // Initialized 
  const Bool_t kIsDijet;

  const Int_t kHardProcessPartonStatus = 23;

  const Int_t kElectronPID_ = 11;
  const Int_t kMuonPID_ = 13;
  const Int_t kGluonPID_ = 21;
  const Int_t kPhotonPID_ = 22;

  // Common
  const Float_t kMinJetPT_ = 30.0f;
  const Float_t kMaxJetEta_ = 2.4f;
  const Float_t kDeltaRCut_ = 0.3f;

  // IsBalanced
  const Float_t kDeltaPhiCut_ = 2.5f;

  // PassZJets


  const Float_t kImageDeltaEtaMax_ = 0.4f;
  const Float_t kImageDeltaPhiMax_ = 0.4f;
};


#endif // DELPHYS_ANALYSIS_QGJETALYSER_H_
