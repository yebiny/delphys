#ifndef doubleHiggsAnalyser_H
#define doubleHiggsAnalyser_H
#include "delphys/oxbridgekinetics/src/Mt2/Basic_Mt2_332_Calculator.h"
#include "delphys/oxbridgekinetics/src/Mt2/ChengHanBisect_Mt2_332_Calculator.h"

class doubleHiggsAnalyser {
private :
  TTree *del_tree;

  TFile *out_file;
  TTree *out_tree;
  
  // Output Variables
  //// MT2 variables
  Float_t lester_MT2 = -99;
  Float_t lester_MT2_b = -99;
  Float_t lester_MT2_l = -99;
  Float_t basic_MT2_332 = -99;
  Float_t basic_MT2_332_b = -99;
  Float_t basic_MT2_332_l = -99;
  Float_t ch_bisect_MT2_332 = -99;
  Float_t ch_bisect_MT2_332_b = -99;
  Float_t ch_bisect_MT2_332_l = -99;
  //// lepton kinematic variables
  Float_t lepton_mass = -99;
  Float_t lepton_pt = -99;
  Float_t lepton_px = -99;
  Float_t lepton_py = -99;
  Float_t lepton_deltaR = -99;
  Float_t lepton1_mass = -99;
  Float_t lepton2_mass = -99;
  Float_t lepton1_pt = -99;
  Float_t lepton2_pt = -99;
  //// bottom kinematic variables
  Float_t bottom_mass = -99;
  Float_t bottom_pt = -99;
  Float_t bottom_px = -99;
  Float_t bottom_py = -99;
  Float_t bottom_deltaR = -99;
  Float_t bottom1_mass = -99;
  Float_t bottom2_mass = -99;
  Float_t bottom1_pt = -99;
  Float_t bottom2_pt = -99;
  //// invariant mass of bbll
  Float_t bbll_mass = -99;
  //// missing et
  Float_t missing_et = -99;
  Float_t missing_et_phi = -99;
  //// truth matching variables
  Int_t lep1_mother = 0;
  Int_t lep2_mother = 0;
  Int_t bot1_mother = 0;
  Int_t bot2_mother = 0;
  Int_t lep1_grmother = 0;
  Int_t lep2_grmother = 0;
  Int_t bot1_grmother = 0;
  Int_t bot2_grmother = 0;

  Int_t step = 0;

  // TLorentzVectors
  TLorentzVector lepton1;
  TLorentzVector lepton2;
  TLorentzVector leptonlepton;
  TLorentzVector bottom1;
  TLorentzVector bottom2;
  TLorentzVector bottombottom;
  TLorentzVector bbll;
  TLorentzVector missing;

  // Delphes Variables
  TClonesArray *particles = 0;
  TClonesArray *missings = 0;
  TClonesArray *jets = 0;
  std::map<Float_t, int, std::greater<Float_t>> leptons;
  std::map<Float_t, int, std::greater<Float_t>>::iterator lepton_iter;
  std::map<Float_t, int, std::greater<Float_t>> bottoms;
  std::map<Float_t, int, std::greater<Float_t>>::iterator bottom_iter;

  // MT2 Calculators
  Mt2::Basic_Mt2_332_Calculator basic_mt2_332Calculator;
  Mt2::ChengHanBisect_Mt2_332_Calculator ch_bisect_mt2_332Calculator;

  double invis_mass     =  100; // GeV
  
public :
  static const int Z_PID = 23;
  static const int Higgs_PID = 25;
  static const int Top_PID = 6;
  static const int Bottom_PID = 5;
  static const int Muon_PID = 13;
  static const int Tau_PID = 17;
  static const int Electron_PID = 11;
  // before loop settings
  void MakeOutputBranch(TTree *tree);
  void SetOutput(TString output_file_name);
  void SetBranchAddress();
  void SetNanoBranchAddress();
  void SetDelphesBranchAddress();
  void Initiate(TString output_file_name);
  // during loop
  void ResetVariables();
  bool Analysis();
    // calculating MT2
    double get_MT2();
    double get_Basic_332_MT2();
  // after loop
  void Finalize();
  virtual void Loop();

  doubleHiggsAnalyser(TTree *tree=0, Bool_t isMC = false);
  doubleHiggsAnalyser(TChain *tchain=0, Bool_t isMC = false);
  ~doubleHiggsAnalyser();
};

// Constructors //////
doubleHiggsAnalyser::doubleHiggsAnalyser(TTree *tree, Bool_t isMC) {
  del_tree = tree;
}

doubleHiggsAnalyser::doubleHiggsAnalyser(TChain *tchain, Bool_t isMC) {
  TTree *tree = dynamic_cast<TTree*>(tchain);
  del_tree = tree;
}

// deconstructors
doubleHiggsAnalyser::~doubleHiggsAnalyser() { }

#endif
