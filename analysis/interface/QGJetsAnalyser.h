#ifndef DELPHYS_ANALYSIS_QGJETALYSER_H_
#define DELPHYS_ANALYSIS_QGJETALYSER_H_

#include "delphys/analysis/interface/BaseAnalyser.h"


class QGJetsAnalyser : private BaseAnalyser {
 public:
  QGJetsAnalyser(const TString & in_path,
                 const TString & out_path,
                 const TString & out_tree_name,
                 Bool_t is_dijet);
  ~QGJetsAnalyser();
  void Loop();

 private:
  // inherited
  void MakeBranch();
  Bool_t SelectEvent();
  void Analyse(Int_t entry);

  void ResetOnEachJet();
  void ResetOnEachEvent();

  Bool_t PassZjets(TClonesArray * jets,
                   TClonesArray * muons,
                   TClonesArray * electrons);
  Bool_t IsBalanced(TClonesArray * jets);
  void FillDaughters(const Jet* jet);
  Int_t Pixelate(Float_t deta, Float_t dphi,
                 Float_t deta_max=0.4, Float_t dphi_max=0.4,
                 Int_t num_deta_bins=33, Int_t num_dphi_bins=33);
  void MakeJetImage();
  Bool_t ExtractSatellites();

  // NOTE
  Int_t m_event;
  Int_t m_num_jets;
  Int_t m_num_good_jets;
  Int_t m_num_primary_vertices;
  Int_t m_order;

  Int_t m_label;

  Float_t m_pt; // jet pt
  Float_t m_eta; // jet eta
  Float_t m_phi; // jet phi
  Float_t m_pt_dr_log;
  Float_t m_ptd; // jet energy sharing variable
  Float_t m_major_axis;
  Float_t m_minor_axis;

  Int_t m_chad_mult;
  Int_t m_nhad_mult;
  Int_t m_electron_mult;
  Int_t m_muon_mult;
  Int_t m_photon_mult;

  Int_t m_parton_id;
  Int_t m_flavor_id;
  Int_t m_flavor_algo_id;
  Int_t m_flavor_phys_id;

  Int_t  m_num_dau;
  Bool_t m_matched;
  // Bool_t m_balanced;
  // Bool_t m_pass_zjets;
  Bool_t m_lepton_overlap; //


  std::vector<TLorentzVector> m_dau_p4;
  std::vector<Float_t>        m_dau_pt;
  std::vector<Float_t>        m_dau_deta;
  std::vector<Float_t>        m_dau_dphi;
  std::vector<Int_t>          m_dau_charge;
  std::vector<Int_t>          m_dau_pid;
  std::vector<Bool_t>         m_dau_is_hadronic;
  std::vector<Bool_t>         m_dau_is_track;
  std::vector<Float_t>        m_dau_eemfrac;
  std::vector<Float_t>        m_dau_ehadfrac;


  // 33*33 = 1089
  Float_t m_image_chad_pt_33[1089];
  Float_t m_image_nhad_pt_33[1089];
  Float_t m_image_electron_pt_33[1089];
  Float_t m_image_muon_pt_33[1089];
  Float_t m_image_photon_pt_33[1089];
  // multiplicity
  Float_t m_image_chad_mult_33[1089];
  Float_t m_image_nhad_mult_33[1089];
  Float_t m_image_electron_mult_33[1089];
  Float_t m_image_muon_mult_33[1089];
  Float_t m_image_photon_mult_33[1089];


  // 
  Bool_t bad_hard_gen_seen_;
  UInt_t num_passed_events_;

  // Constants
  const Bool_t kIsDijet;
  const Bool_t kIsQuarkJet;

  const Int_t kHardProcessPartonStatus = 23;

  const Int_t kElectronPID = 11;
  const Int_t kMuonPID = 13;
  const Int_t kGluonPID = 21;
  const Int_t kPhotonPID = 22;

  const Float_t kMinJetPT = 20.0;
  const Float_t kMaxJetEta = 2.4;
  const Float_t kDeltaRCut = 0.3;

  const Float_t kImageDeltaEtaMax = 0.4;
  const Float_t kImageDeltaPhiMax = 0.4;
};


#endif // DELPHYS_ANALYSIS_QGJETALYSER_H_
