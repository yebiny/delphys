#ifndef doubleHiggsAnalyser_H
#define doubleHiggsAnalyser_H
#include "delphys/oxbridgekinetics/src/Mt2/Basic_Mt2_332_Calculator.h"
#include "delphys/oxbridgekinetics/src/Mt2/ChengHanBisect_Mt2_332_Calculator.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCuts.h"
#include "TMinuit.h"
#include "TError.h"

class doubleHiggsAnalyser {
private :
  TTree *del_tree;

  TFile *out_file;
  TTree *out_tree;
  TTree *n_events_tree;
  
  Int_t event_size = 1;
  
  // Output Variables
  //// MT2 variables
  Float_t lester_MT2 = -99;
  Float_t lester_MT2_b = -99;
  Float_t lester_MT2_l = -99;
  Float_t basic_MT2_332_bbll = -99;
  Float_t basic_MT2_332_b = -99;
  Float_t basic_MT2_332_l = -99;
  Float_t basic_MT2_332_blbl = -99;
  Float_t ch_bisect_MT2_332 = -99;
  Float_t ch_bisect_MT2_332_b = -99;
  Float_t ch_bisect_MT2_332_l = -99;
  //// mT = sqrt(2*pt(ll)*MET*[1-cos(phi(ll)-phi(MET))])
  Float_t mT = -99;
  //// lepton kinematic variables
  Float_t ll_deltaR = -99;
  Float_t ll_deltaPhi = -99;
  Float_t ll_Pt = -99;
  Float_t ll_M = -99;
  //// bottom kinematic variables
  Float_t bb_deltaR = -99;
  Float_t bb_deltaPhi = -99;
  Float_t bb_Pt = -99;
  Float_t bb_M = -99;
  //// lepton and bottom kinematic variables
  std::vector<Float_t> bl_deltaR;
  Float_t bl_min_deltaR = -99;
  Float_t bbll_deltaR = -99;
  Float_t bbll_deltaPhi = -99;
  //// truth matching variables
  Int_t lep1_mother = 0;
  Int_t lep2_mother = 0;
  Int_t bot1_mother = 0;
  Int_t bot2_mother = 0;
  Int_t lep1_grmother = 0;
  Int_t lep2_grmother = 0;
  Int_t bot1_grmother = 0;
  Int_t bot2_grmother = 0;
  //// topness and higgsness
  Float_t topness = 0;
  Float_t higgsness = 0;

  Int_t step = 0;

  ////TMVA variables
  //reader
  TMVA::Reader* bdtg_reader;
  Float_t tmva_bdtg_output = 0;

  // TLorentzVectors
  TLorentzVector lepton1;
  TLorentzVector lepton2;
  TLorentzVector leptonlepton;
  TLorentzVector bottom1;
  TLorentzVector bottom2;
  TLorentzVector bottombottom;
  TLorentzVector bottomlepton1;
  TLorentzVector bottomlepton2;
  TLorentzVector bbll;
  TLorentzVector missing;

  // Delphes Variables
  TClonesArray *particles = 0;
  TClonesArray *muons = 0;
  TClonesArray *electrons = 0;
  TClonesArray *missings = 0;
  TClonesArray *jets = 0;
  std::map<Float_t, std::pair<int,int>, std::greater<Float_t>> leptons;
  std::map<Float_t, std::pair<int,int>, std::greater<Float_t>>::iterator lepton_iter;
  std::map<Float_t, int, std::greater<Float_t>> bottoms;
  std::map<Float_t, int, std::greater<Float_t>>::iterator bottom_iter;

  // MT2 Calculators
  Mt2::Basic_Mt2_332_Calculator basic_mt2_332Calculator;
  Mt2::ChengHanBisect_Mt2_332_Calculator ch_bisect_mt2_332Calculator;
  
  TMinuit *ptMinuit;
  TMinuit *ptMinuitT1;
  TMinuit *ptMinuitT2;
  Double_t minuit_arglist[10];
  Int_t ierflag = 0;

  double invis_mass     =  100; // GeV
  
public :
  static const int Z_PID = 23;
  static const int Higgs_PID = 25;
  static const int Top_PID = 6;
  static const int Bottom_PID = 5;
  static const int Electron_PID = 11;
  static const int Muon_PID = 13;
  static const int Tau_PID = 15;
  static const int n_unknown_par = 4;
  static constexpr double minuit_pmin = -5000;
  static constexpr double minuit_pmax = 5000;
  static constexpr float Muon_Mass = 105.6583745; // MeV
  static constexpr float Electron_Mass = 0.5109989461; // MeV

  // before loop settings
  void MakeOutputBranch(TTree *tree);
  void SetOutput(TString output_file_name);
  void SetTMVA(TString weight_file_path);
  void SetBranchAddress();
  void SetNanoBranchAddress();
  void SetDelphesBranchAddress();
  void GetHiggsness();
  void GetTopness();
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
