#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TRefArray.h>
#include <TMatrixDfwd.h>
#include <TVectorD.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include "TMinuit.h"
#include "TError.h"

#include "delphys/analysis/interface/doubleHiggsAnalyser.h"
#include "classes/DelphesClasses.h"
#include "delphys/external/interface/lester_mt2_bisect.h"

#include "delphys/oxbridgekinetics/src/Mt2/Basic_Mt2_332_Calculator.h"
#include "delphys/oxbridgekinetics/src/Mt2/ChengHanBisect_Mt2_332_Calculator.h"

using namespace std;

  const double mtop = 173.;   // <------- top mass
  const double mh = 125;      // <------- Higgs mass
  const double mw = 80.419;   // <------- W mass
  const double mhw = 30;      // <------- Peak of off-shell W* ~ mh - mw

  const double sigmahlep = 2.;  // <------- sigma_{h_l}
  const double sigmaWon = 5.;   // <------- sigma_W
  const double sigmaWoff = 5.;  // <------- sigma_W*
  const double sigmat = 5.;     // <------- sigma_t
  
  TLorentzVector g_missing, g_lepton1, g_lepton2, g_bottom1, g_bottom2;

double  GetMass(double Px, double Py, double Pz, double E)
{
  double f = E*E -Px*Px -Py*Py -Pz*Pz;
  double fout = 0;

  if ( f < 0. ) {
    fout = 0;
  } else {
    fout = sqrt(f);
  }
  return fout;
}

// return mother particle(GenParticle) of particle 'p' among 'particles'.
GenParticle* getMother(TClonesArray* particles, const GenParticle* p){
  if (p->M1==-1) return nullptr;
  auto mom = static_cast<GenParticle *>(particles->At(p->M1));
  while (mom->PID == p->PID){
    mom = static_cast<GenParticle *>(particles->At(mom->M1));
    if (mom->M1==-1) return nullptr;
  }
  return mom;
}

// return the PID of the mother particle of the particle at 'ip' among 'particles'
std::pair<int,int> isFrom(TClonesArray* particles, int ip){
   auto p = static_cast<GenParticle *>(particles->At(ip));
   // check if it's from Higgs
   auto mom = getMother(particles, p); 
   if (mom==nullptr) return std::make_pair(0,0);
   auto grmom = getMother(particles, mom);
   if (grmom==nullptr) return std::make_pair(0,0);
   auto pedigree = std::make_pair(mom->PID, grmom->PID);
   return pedigree;
}

void doubleHiggsAnalyser::MakeOutputBranch(TTree *tree) {
  n_events_tree->Branch("event_size",&event_size,"event_size/I");
  // MT2 variables
  tree->Branch("lester_MT2",&lester_MT2,"lester_MT2/F");
  tree->Branch("basic_MT2_332_bbll",&basic_MT2_332_bbll,"basic_MT2_332_bbll/F");
  tree->Branch("ch_bisect_MT2_332",&ch_bisect_MT2_332,"ch_bisect_MT2_332/F");
  tree->Branch("basic_MT2_332_blbl",&basic_MT2_332_blbl,"basic_MT2_332_blbl/F");
  // MT2(b) variables
  tree->Branch("lester_MT2_b",&lester_MT2_b,"lester_MT2_b/F");
  tree->Branch("basic_MT2_332_b",&basic_MT2_332_b,"basic_MT2_332_b/F");
  tree->Branch("ch_bisect_MT2_332_b",&ch_bisect_MT2_332_b,"ch_bisect_MT2_332_b/F");
  // MT2(l) variables
  tree->Branch("lester_MT2_l",&lester_MT2_l,"lester_MT2_l/F");
  tree->Branch("basic_MT2_332_l",&basic_MT2_332_l,"basic_MT2_332_l/F");
  tree->Branch("ch_bisect_MT2_332_l",&ch_bisect_MT2_332_l,"ch_bisect_MT2_332_l/F");
  // mT = sqrt(2*pt(ll)*MET*[1-cos(phi(ll)-phi(MET))])
  tree->Branch("mT",&mT,"mT/F");
  // lepton kinematic variables
  tree->Branch("lep1","TLorentzVector", &lepton1);
  tree->Branch("lep2","TLorentzVector", &lepton2);
  tree->Branch("ll","TLorentzVector", &leptonlepton);
  tree->Branch("ll_deltaR",&ll_deltaR,"ll_deltaR/F");
  tree->Branch("ll_deltaPhi",&ll_deltaPhi,"ll_deltaPhi/F");
  // bottom kinematic variables
  tree->Branch("bot1","TLorentzVector", &bottom1);
  tree->Branch("bot2","TLorentzVector", &bottom2);
  tree->Branch("bb","TLorentzVector",&bottombottom);
  tree->Branch("bb_deltaR",&bb_deltaR,"bb_deltaR/F");
  tree->Branch("bb_deltaPhi",&bb_deltaPhi,"bb_deltaPhi/F");
  // lepton bottom kinematic variables
  tree->Branch("bl_deltaR","vector<Float_t>",&bl_deltaR);
  tree->Branch("bl_min_deltaR",&bl_min_deltaR,"bl_min_deltaR/F");
  tree->Branch("bbll_deltaR",&bbll_deltaR,"bbll_deltaR/F");
  tree->Branch("bbll_deltaPhi",&bbll_deltaPhi,"bbll_deltaPhi/F");
  // missing et
  tree->Branch("MET","TLorentzVector", &missing);
  // invariant mass of bbll
  tree->Branch("bbll","TLorentzVector",&bbll);
  // topness and higgsness
  tree->Branch("topness",&topness,"topness/F");
  tree->Branch("higgsness",&higgsness,"higgsness/F");

  // truth matching variables
  tree->Branch("lepton1MotherPID",&lep1_mother,"lep1_mother/I");
  tree->Branch("lepton2MotherPID",&lep2_mother,"lep2_mother/I");
  tree->Branch("lepton1GrMotherPID",&lep1_grmother,"lep1_grmother/I");
  tree->Branch("lepton2GrMotherPID",&lep2_grmother,"lep2_grmother/I");
  tree->Branch("bottom1MotherPID",&bot1_mother,"bot1_mother/I");
  tree->Branch("bottom2MotherPID",&bot2_mother,"bot2_mother/I");
  tree->Branch("bottom1GrMotherPID",&bot1_grmother,"bot1_grmother/I");
  tree->Branch("bottom2GrMotherPID",&bot2_grmother,"bot2_grmother/I");

  // cut flow
  tree->Branch("step",&step,"step/I");

  // tmva
  tree->Branch("tmva_bdtg_output",&tmva_bdtg_output,"tmva_bdtg_output/F");
}

void doubleHiggsAnalyser::SetOutput(TString output_file_name) {
  out_file = TFile::Open(output_file_name,"RECREATE");
  out_tree = new TTree("events","events"); 
  n_events_tree = new TTree("nevents","nevents");
}

void doubleHiggsAnalyser::SetTMVA(TString weight_file_path) {
  bdtg_reader = new TMVA::Reader();
  bdtg_reader->AddVariable("ll_deltaR", &ll_deltaR);
  bdtg_reader->AddVariable("ll.Pt()", &ll_Pt);
  bdtg_reader->AddVariable("ll.M()", &ll_M);
  bdtg_reader->AddVariable("bb_deltaR", &bb_deltaR);
  bdtg_reader->AddVariable("bb.Pt()", &bb_Pt);
  bdtg_reader->AddVariable("bb.M()", &bb_M);
  bdtg_reader->AddVariable("bl_min_deltaR", &bl_min_deltaR);
  bdtg_reader->AddVariable("bbll_deltaR", &bbll_deltaR);
  bdtg_reader->AddVariable("bbll_deltaPhi", &bbll_deltaPhi);
  bdtg_reader->AddVariable("mT", &mT);
  bdtg_reader->AddVariable("basic_MT2_332_bbll", &basic_MT2_332_bbll);
  bdtg_reader->BookMVA("BDTG", weight_file_path);
}

void doubleHiggsAnalyser::SetBranchAddress() {
  del_tree->SetBranchAddress("Particle",&particles);
  del_tree->SetBranchAddress("Muon",&muons);
  del_tree->SetBranchAddress("Electron",&electrons);
  del_tree->SetBranchAddress("MissingET",&missings);
  del_tree->SetBranchAddress("Jet",&jets);
}

void fcnT1(  int& npar, double* deriv, double& f, double par[], int flag) {
  double chisq = 0;
  double Nu1Px = par[0];
  double Nu1Py = par[1];
  double Nu1Pz = par[2];
  double Nu1E = sqrt( Nu1Px*Nu1Px + Nu1Py*Nu1Py + Nu1Pz*Nu1Pz );
  
  double Nu2Px = g_missing.Px() - Nu1Px;
  double Nu2Py = g_missing.Py() - Nu1Py;
  double Nu2Pz = par[3];
  double Nu2E = sqrt( Nu2Px*Nu2Px + Nu2Py*Nu2Py + Nu2Pz*Nu2Pz );

  double Lep1E = g_lepton1.E();
  double Lep1Px = g_lepton1.Px();
  double Lep1Py = g_lepton1.Py();
  double Lep1Pz = g_lepton1.Pz();

  double Lep2E = g_lepton2.E();
  double Lep2Px = g_lepton2.Px();
  double Lep2Py = g_lepton2.Py();
  double Lep2Pz = g_lepton2.Pz();

  double Jet1E = g_bottom1.E();
  double Jet1Px = g_bottom1.Px();
  double Jet1Py = g_bottom1.Py();
  double Jet1Pz = g_bottom1.Pz();

  double Jet2E = g_bottom2.E();
  double Jet2Px = g_bottom2.Px();
  double Jet2Py = g_bottom2.Py();
  double Jet2Pz = g_bottom2.Pz();

  double jln111Px = Jet1Px + Lep1Px + Nu1Px;
  double jln111Py = Jet1Py + Lep1Py + Nu1Py;
  double jln111Pz = Jet1Pz + Lep1Pz + Nu1Pz;
  double jln111E = Jet1E + Lep1E + Nu1E;

  double mjln111 = GetMass( jln111Px, jln111Py, jln111Pz, jln111E );

  double deltaJLN111 = (mjln111*mjln111 - mtop*mtop) / (sigmat*sigmat);

  double jln222Px = Jet2Px + Lep2Px + Nu2Px;
  double jln222Py = Jet2Py + Lep2Py + Nu2Py;
  double jln222Pz = Jet2Pz + Lep2Pz + Nu2Pz;
  double jln222E = Jet2E + Lep2E + Nu2E;

  double mjln222 = GetMass( jln222Px, jln222Py, jln222Pz, jln222E );

  double deltaJLN222 = (mjln222*mjln222 - mtop*mtop) / (sigmat*sigmat);

  double ln11Px = Lep1Px + Nu1Px;
  double ln11Py = Lep1Py + Nu1Py;
  double ln11Pz = Lep1Pz + Nu1Pz;
  double ln11E = Lep1E + Nu1E;

  double mln11 = GetMass( ln11Px, ln11Py, ln11Pz, ln11E );
  double deltaLN11on = (mln11*mln11 - mw*mw) / (sigmaWon*sigmaWon);

  double ln22Px = Lep2Px + Nu2Px;
  double ln22Py = Lep2Py + Nu2Py;
  double ln22Pz = Lep2Pz + Nu2Pz;
  double ln22E = Lep2E + Nu2E;

  double mln22 = GetMass( ln22Px, ln22Py, ln22Pz, ln22E );

  double deltaLN22on = (mln22*mln22 - mw*mw) / (sigmaWon*sigmaWon);

  chisq = deltaJLN111*deltaJLN111 + deltaJLN222*deltaJLN222 + deltaLN11on*deltaLN11on + deltaLN22on*deltaLN22on;

    cout << nMuon << endl;
chisq;
}

void fcnT2(  int& npar, double* deriv, double& f, double par[], int flag) {
  double chisq = 0;
  double Nu1Px = par[0];
  double Nu1Py = par[1];
  double Nu1Pz = par[2];
  double Nu1E = sqrt( Nu1Px*Nu1Px + Nu1Py*Nu1Py + Nu1Pz*Nu1Pz );
  
  double Nu2Px = g_missing.Px() - Nu1Px;
  double Nu2Py = g_missing.Py() - Nu1Py;
  double Nu2Pz = par[3];
  double Nu2E = sqrt( Nu2Px*Nu2Px + Nu2Py*Nu2Py + Nu2Pz*Nu2Pz );

  double Lep1E = g_lepton1.E();
  double Lep1Px = g_lepton1.Px();
  double Lep1Py = g_lepton1.Py();
  double Lep1Pz = g_lepton1.Pz();

  double Lep2E = g_lepton2.E();
  double Lep2Px = g_lepton2.Px();
  double Lep2Py = g_lepton2.Py();
  double Lep2Pz = g_lepton2.Pz();

  double Jet1E = g_bottom1.E();
  double Jet1Px = g_bottom1.Px();
  double Jet1Py = g_bottom1.Py();
  double Jet1Pz = g_bottom1.Pz();

  double Jet2E = g_bottom2.E();
  double Jet2Px = g_bottom2.Px();
  double Jet2Py = g_bottom2.Py();
  double Jet2Pz = g_bottom2.Pz();

  double jln121Px = Jet1Px + Lep2Px + Nu1Px;
  double jln121Py = Jet1Py + Lep2Py + Nu1Py;
  double jln121Pz = Jet1Pz + Lep2Pz + Nu1Pz;
  double jln121E = Jet1E + Lep2E + Nu1E;

  double mjln121 = GetMass( jln121Px, jln121Py, jln121Pz, jln121E );

  double deltaJLN121 = (mjln121*mjln121 - mtop*mtop) / (sigmat*sigmat);

  double jln212Px = Jet2Px + Lep1Px + Nu2Px;
  double jln212Py = Jet2Py + Lep1Py + Nu2Py;
  double jln212Pz = Jet2Pz + Lep1Pz + Nu2Pz;
  double jln212E = Jet2E + Lep1E + Nu2E;

  double mjln212 = GetMass( jln212Px, jln212Py, jln212Pz, jln212E );

  double deltaJLN212 = (mjln212*mjln212 - mtop*mtop) / (sigmat*sigmat);
  
  double ln21Px = Lep2Px + Nu1Px;
  double ln21Py = Lep2Py + Nu1Py;
  double ln21Pz = Lep2Pz + Nu1Pz;
  double ln21E = Lep2E + Nu1E;

  double mln21 = GetMass( ln21Px, ln21Py, ln21Pz, ln21E );
  double deltaLN21on = (mln21*mln21 - mw*mw) / (sigmaWon*sigmaWon);

  double ln12Px = Lep1Px + Nu2Px;
  double ln12Py = Lep1Py + Nu2Py;
  double ln12Pz = Lep1Pz + Nu2Pz;
  double ln12E = Lep1E + Nu2E;

  double mln12 = GetMass( ln12Px, ln12Py, ln12Pz, ln12E );

  double deltaLN12on = (mln12*mln12 - mw*mw) / (sigmaWon*sigmaWon);

  chisq = deltaJLN121*deltaJLN121 + deltaJLN212*deltaJLN212 + deltaLN21on*deltaLN21on + deltaLN12on*deltaLN12on;

  f = chisq;
}

void fcnH(  int& npar, double* deriv, double& f, double par[], int flag) {
  double chisq = 0;
  double Nu1Px = par[0];
  double Nu1Py = par[1];
  double Nu1Pz = par[2];
  double Nu1E = sqrt( Nu1Px*Nu1Px + Nu1Py*Nu1Py + Nu1Pz*Nu1Pz );
  
  double Nu2Px = g_missing.Px() - Nu1Px;
  double Nu2Py = g_missing.Py() - Nu1Py;
  double Nu2Pz = par[3];
  double Nu2E = sqrt( Nu2Px*Nu2Px + Nu2Py*Nu2Py + Nu2Pz*Nu2Pz );

  double Lep1E = g_lepton1.E();
  double Lep1Px = g_lepton1.Px();
  double Lep1Py = g_lepton1.Py();
  double Lep1Pz = g_lepton1.Pz();

  double Lep2E = g_lepton2.E();
  double Lep2Px = g_lepton2.Px();
  double Lep2Py = g_lepton2.Py();
  double Lep2Pz = g_lepton2.Pz();

  double llnnPx = Lep1Px + Lep2Px + Nu1Px + Nu2Px;
  double llnnPy = Lep1Py + Lep2Py + Nu1Py + Nu2Py;
  double llnnPz = Lep1Pz + Lep2Pz + Nu1Pz + Nu2Pz;
  double llnnE = Lep1E + Lep2E + Nu1E + Nu2E;

  double mllnn = GetMass( llnnPx, llnnPy, llnnPz, llnnE );
  
  double deltaLLNN = (mllnn*mllnn - mh*mh) / (sigmahlep*sigmahlep);

  double ln11Px = Lep1Px + Nu1Px;
  double ln11Py = Lep1Py + Nu1Py;
  double ln11Pz = Lep1Pz + Nu1Pz;
  double ln11E = Lep1E + Nu1E;
 
  double mln11 = GetMass( ln11Px, ln11Py, ln11Pz, ln11E );

  double deltaLN11on = (mln11*mln11 - mw*mw) / (sigmaWon*sigmaWon);

  double deltaLN11off = (mln11*mln11 - mhw*mhw) / (sigmaWoff*sigmaWoff);

  double ln22Px = Lep2Px + Nu2Px;
  double ln22Py = Lep2Py + Nu2Py;
  double ln22Pz = Lep2Pz + Nu2Pz;
  double ln22E = Lep2E + Nu2E;

  double mln22 = GetMass( ln22Px, ln22Py, ln22Pz, ln22E );

  double deltaLN22on = (mln22*mln22 - mw*mw) / (sigmaWon*sigmaWon);

  double deltaLN22off = (mln22*mln22 - mhw*mhw) / (sigmaWoff*sigmaWoff);

  double ln12Px = Lep1Px + Nu2Px;
  double ln12Py = Lep1Py + Nu2Py;
  double ln12Pz = Lep1Pz + Nu2Pz;
  double ln12E = Lep1E + Nu2E;

  double mln12 = GetMass( ln12Px, ln12Py, ln12Pz, ln12E );

  double deltaLN12on = (mln12*mln12 - mw*mw) / (sigmaWoff*sigmaWoff);
  double deltaLN12off = (mln12*mln12 - mhw*mhw) / (sigmaWoff*sigmaWoff);

  double ln21Px = Lep2Px + Nu1Px;
  double ln21Py = Lep2Py + Nu1Py;
  double ln21Pz = Lep2Pz + Nu2Pz;
  double ln21E = Lep2E + Nu1E;

  double mln21 = GetMass( ln21Px, ln21Py, ln21Pz, ln21E );

  double deltaLN21on = (mln21*mln21 - mw*mw) / (sigmaWon*sigmaWon);
  double deltaLN21off = (mln21*mln21 - mhw*mhw) / (sigmaWoff*sigmaWoff);

 double  deltaCase1 = deltaLN11on*deltaLN11on + deltaLN22off*deltaLN22off;
 double  deltaCase2 = deltaLN12on*deltaLN12on + deltaLN21off*deltaLN21off;
 double  deltaCase3 = deltaLN22on*deltaLN22on + deltaLN11off*deltaLN11off;
 double  deltaCase4 = deltaLN21on*deltaLN21on + deltaLN12off*deltaLN12off;

  double deltaCaseTotal[] = {deltaCase1, deltaCase2, deltaCase3, deltaCase4};
  double deltaCaseMin = *min_element(deltaCaseTotal, deltaCaseTotal+4);

  //double nnPx = Nu1Px + Nu2Px;
  //double nnPy = Nu1Py + Nu2Py;
  //double nnPz = Nu1Pz + Nu2Pz;
  //double nnE = Nu1E + Nu2E;

  //double mnn = GetMass( nnPx, nnPy, nnPz, nnE );

  chisq = deltaLLNN*deltaLLNN + deltaCaseMin;

  f = chisq;
}

void doubleHiggsAnalyser::GetHiggsness() {
  ptMinuit = new TMinuit(n_unknown_par);
  ptMinuit->SetPrintLevel(-1);
  ptMinuit->SetFCN(fcnH);
  minuit_arglist[0] = 1;
  ptMinuit->mnexcm("SET ERR", minuit_arglist, 1, ierflag);
  Double_t vstart[4] = {0.0, 0.0, 0.0, 0.0};
  Double_t step[4] = {0.1, 0.1, 0.1, 0.1};
  TString par_name[4] = {"Nu1Px", "Nu1Py", "Nu1Pz", "Nu2Pz"};
  for (int i=0; i<4; i++) {
    ptMinuit->mnparm(i, par_name[i], vstart[i], step[i], minuit_pmin, minuit_pmax, ierflag);
  }

  minuit_arglist[0] = 500;
  minuit_arglist[1] = 1.;
  
  ptMinuit->mnexcm("MIGRAD", minuit_arglist, 2, ierflag); 

  double outparH[n_unknown_par], errH[n_unknown_par];
  for (int i = 0; i<4; i++) {
    ptMinuit->GetParameter(i, outparH[i], errH[i]);
  }

  Double_t aminH, edmH, errdefH;
  Int_t nvparH, nparxH, icstatH;

  ptMinuit->mnstat(aminH, edmH, errdefH, nvparH, nparxH, icstatH);
  if ( aminH > 0 ) {
    higgsness = log10(aminH);
  }
  delete ptMinuit;
}

void doubleHiggsAnalyser::GetTopness() {
  ptMinuitT1 = new TMinuit(n_unknown_par);
  ptMinuitT1->SetPrintLevel(-1);
  ptMinuitT1->SetFCN(fcnT1);
  minuit_arglist[0] = 1;
  ptMinuitT1->mnexcm("SET ERR", minuit_arglist, 1, ierflag);
  Double_t vstart[4] = {0.0, 0.0, 0.0, 0.0};
  Double_t step[4] = {0.1, 0.1, 0.1, 0.1};
  TString par_name[4] = {"Nu1Px", "Nu1Py", "Nu1Pz", "Nu2Pz"};
  for (int i=0; i<4; i++) {
    ptMinuitT1->mnparm(i, par_name[i], vstart[i], step[i], minuit_pmin, minuit_pmax, ierflag);
  }

  minuit_arglist[0] = 500;
  minuit_arglist[1] = 1.;
  
  ptMinuitT1->mnexcm("MIGRAD", minuit_arglist, 2, ierflag); 

  double outparT1[n_unknown_par], errT1[n_unknown_par];
  for (int i = 0; i<4; i++) {
    ptMinuitT1->GetParameter(i, outparT1[i], errT1[i]);
  }

  Double_t aminT1, edmT1, errdefT1;
  Int_t nvparT1, nparxT1, icstatT1;

  ptMinuitT1->mnstat(aminT1, edmT1, errdefT1, nvparT1, nparxT1, icstatT1);
  
  delete ptMinuitT1;
  
  ptMinuitT2 = new TMinuit(n_unknown_par);
  ptMinuitT2->SetPrintLevel(-1);
  ptMinuitT2->SetFCN(fcnT2);
  minuit_arglist[0] = 1;
  ptMinuitT2->mnexcm("SET ERR", minuit_arglist, 1, ierflag);
  //vstart[4] = {0.0, 0.0, 0.0, 0.0};
  //step[4] = {0.1, 0.1, 0.1, 0.1};
  //par_name[4] = {"Nu1Px", "Nu1Py", "Nu1Pz", "Nu2Pz"};
  for (int i=0; i<4; i++) {
    ptMinuitT2->mnparm(i, par_name[i], vstart[i], step[i], minuit_pmin, minuit_pmax, ierflag);
  }

  minuit_arglist[0] = 500;
  minuit_arglist[1] = 1.;
  
  ptMinuitT2->mnexcm("MIGRAD", minuit_arglist, 2, ierflag); 

  double outparT2[n_unknown_par], errT2[n_unknown_par];
  for (int i = 0; i<4; i++) {
    ptMinuitT2->GetParameter(i, outparT2[i], errT2[i]);
  }

  Double_t aminT2, edmT2, errdefT2;
  Int_t nvparT2, nparxT2, icstatT2;

  ptMinuitT2->mnstat(aminT2, edmT2, errdefT2, nvparT2, nparxT2, icstatT2);
  
  delete ptMinuitT2;

  Double_t aminT;
  aminT = min(aminT1, aminT2);

  if( aminT > 0 )
  {
    topness = log10(aminT);
  }
}

void doubleHiggsAnalyser::Initiate(TString output_file_name) {
  // set output file
  doubleHiggsAnalyser::SetOutput(output_file_name);
  // make output branch
  doubleHiggsAnalyser::MakeOutputBranch(out_tree);
  // set branch address
  doubleHiggsAnalyser::SetBranchAddress();
}

void doubleHiggsAnalyser::ResetVariables() {
  // MT2 variable
  lester_MT2 = -99;
  lester_MT2_b = -99;
  lester_MT2_l = -99;
  basic_MT2_332_bbll = -99;
  basic_MT2_332_b = -99;
  basic_MT2_332_l = -99;
  basic_MT2_332_blbl = -99;
  ch_bisect_MT2_332 = -99;
  ch_bisect_MT2_332_b = -99;
  ch_bisect_MT2_332_l = -99;

  // lepton kinematic variables
  ll_deltaR = -99;
  ll_deltaPhi = -99;

  // bottom kinematic variables
  bb_deltaR = -99;
  bb_deltaPhi = -99;

  // lepton bottom kinematic variables
  bl_deltaR.clear();
  bl_min_deltaR = -99;
  bbll_deltaR = -99;
  bbll_deltaPhi = -99;

  // mT
  mT = -99;

  // truth matching variables 
  lep1_mother = 0;
  lep2_mother = 0;
  bot1_mother = 0;
  bot2_mother = 0;
  lep1_grmother = 0;
  lep2_grmother = 0;
  bot1_grmother = 0;
  bot2_grmother = 0;

  // cut flow
  step = 0;

  topness = std::numeric_limits<double>::quiet_NaN();
  higgsness = std::numeric_limits<double>::quiet_NaN();
  
  // lepton and bottom maps
  leptons.clear();
  bottoms.clear();
}

bool doubleHiggsAnalyser::Analysis() {
    doubleHiggsAnalyser::ResetVariables();
  //map<Float_t, int, greater<Float_t>> leptons : map of <pt,index>:<K,V> of muon sorted by pt.
  //base selections : MissingET > 20, pT(lepton) > 20, deltaR(ll) < 1.0, m(ll) < 65, deltaR(bb) < 1.3, 95 < m(bb) < 140
    
    // Missing ET
    auto m = static_cast<const MissingET *>(missings->At(0)); // There is always one MET object.
    missing.SetPtEtaPhiM(m->MET,0,m->Phi,0);
    // missing_et baseline selection
    if (missing.E()<20) {return false;}
    
    // collet leptons
    // muons
    for (int im = 0; im < muons->GetEntries(); im++){
      auto m = static_cast<const Muon *>(muons->At(im));
      if (fabs(m->Eta) > 2.4 || fabs(m->PT) < 10) continue;
      leptons.insert(make_pair(m->PT,make_pair(doubleHiggsAnalyser::Muon_PID*m->Charge,im)));
    }
    // electrons
    for (int ie = 0; ie < electrons->GetEntries(); ie++){
      auto e = static_cast<const Electron *>(electrons->At(ie));
      if (fabs(e->Eta) > 2.4 || fabs(e->PT) < 10) continue;
      leptons.insert(make_pair(e->PT,make_pair(doubleHiggsAnalyser::Electron_PID*e->Charge,ie)));
    }
    //for (int ip = 0; ip < particles->GetEntries(); ip++){
    //  auto p = static_cast<const GenParticle *>(particles->At(ip));
    //  // if (abs(p->PID)!=doubleHiggsAnalyser::Tau_PID)
    //  if (abs(p->PID)!=doubleHiggsAnalyser::Electron_PID && abs(p->PID)!=doubleHiggsAnalyser::Muon_PID) continue;
    //  if (fabs(p->Eta) > 2.4 || fabs(p->PT) < 20) continue;
    //  leptons.insert(make_pair(p->PT,ip));
    //}

    if (leptons.size()<2) {
      return false;
    }
    
    lepton_iter = leptons.begin();
    auto l1_info = lepton_iter->second;
    int pid1 = l1_info.first;
    int index_l1 = l1_info.second;
    lepton_iter++;
    auto l2_info = lepton_iter->second;
    int pid2 = l2_info.first;
    int index_l2 = l2_info.second;
    lepton_iter++;
    while (pid1*pid2 > 0 && (lepton_iter!=leptons.end())){
      l2_info = lepton_iter->second;
      pid2 = l2_info.first;
      index_l2 = l2_info.second;
      lepton_iter++;
    }
    if (pid1*pid2 > 0) return false;
    if (abs(pid1)==doubleHiggsAnalyser::Electron_PID){
      auto lep1 = static_cast<const Electron *>(electrons->At(index_l1));
      lepton1.SetPtEtaPhiM(lep1->PT,lep1->Eta,lep1->Phi,doubleHiggsAnalyser::Electron_Mass);
    } else {
      auto lep1 = static_cast<const Muon *>(muons->At(index_l1));
      lepton1.SetPtEtaPhiM(lep1->PT,lep1->Eta,lep1->Phi,doubleHiggsAnalyser::Muon_Mass);
    }
    if (abs(pid2)==doubleHiggsAnalyser::Electron_PID){
      auto lep2 = static_cast<const Electron *>(electrons->At(index_l2));
      lepton2.SetPtEtaPhiM(lep2->PT,lep2->Eta,lep2->Phi,doubleHiggsAnalyser::Electron_Mass);
    } else {
      auto lep2 = static_cast<const Muon *>(muons->At(index_l2));
      lepton2.SetPtEtaPhiM(lep2->PT,lep2->Eta,lep2->Phi,doubleHiggsAnalyser::Muon_Mass);
    }
    //auto lep1 = static_cast<const GenParticle *>(particles->At(lepton_iter->second));
    //lepton1.SetPtEtaPhiM(lep1->PT,lep1->Eta,lep1->Phi,lep1->Mass);
    //auto lepton1_pdg = isFrom(particles, lepton_iter->second); // lepton truth matching
    //lep1_mother = abs(lepton1_pdg.first);
    //lep1_grmother = abs(lepton1_pdg.second);
    //++lepton_iter;
    //lepton2.SetPtEtaPhiM(lep2->PT,lep2->Eta,lep2->Phi,lep2->Mass);
    //auto lepton2_pdg = isFrom(particles, lepton_iter->second); // lepton truth matching
    //lep2_mother = abs(lepton2_pdg.first);
    //lep2_grmother = abs(lepton2_pdg.second);
    leptonlepton = lepton1+lepton2;
    ll_M = leptonlepton.M();
    ll_Pt = leptonlepton.Pt();
    ll_deltaR = fabs(lepton1.DeltaR(lepton2));
    ll_deltaPhi = fabs(lepton1.DeltaPhi(lepton2));
     
    // collect b jets
    for (int ij = 0; ij < jets->GetEntries(); ij++){
      auto jet = static_cast<const Jet *>(jets->At(ij));
      if (abs(jet->Flavor) != doubleHiggsAnalyser::Bottom_PID || fabs(jet->PT) < 30 || fabs(jet->Eta) > 2.4) continue;
      bottoms.insert(make_pair(jet->PT,ij));
    }
   
    if (bottoms.size()<2) {
      return false;
    }
    bottom_iter = bottoms.begin();
    auto bot1 = static_cast<const Jet *>(jets->At(bottom_iter->second));
    bottom1.SetPtEtaPhiM(bot1->PT,bot1->Eta,bot1->Phi,bot1->Mass);
    //auto bottom1_pdg = isFrom(particles, bottom_iter->second); // bottom truth matching
    //bot1_mother = abs(bottom1_pdg.first);
    //bot1_grmother = abs(bottom1_pdg.second);
    bottom_iter ++;
    auto bot2 = static_cast<const Jet *>(jets->At(bottom_iter->second));
    bottom2.SetPtEtaPhiM(bot2->PT,bot2->Eta,bot2->Phi,bot2->Mass);
    bottombottom = bottom1+bottom2;
    bb_M = bottombottom.M();
    bb_Pt = bottombottom.Pt();
    //auto bottom2_pdg = isFrom(particles, bottom_iter->second); // bottom truth matching
    //bot2_mother = abs(bottom2_pdg.first);
    //bot2_grmother = abs(bottom2_pdg.second);

    // bottom kinematic variables
    bb_deltaR = fabs(bottom1.DeltaR(bottom2));
    bb_deltaPhi = fabs(bottom1.DeltaPhi(bottom2));
    // bbll
    bbll = bottombottom + leptonlepton;
    // lepton and bottom kinematic variables
    bl_deltaR.push_back(lepton1.DeltaR(bottom1));
    bl_deltaR.push_back(lepton2.DeltaR(bottom1));
    bl_deltaR.push_back(lepton1.DeltaR(bottom2));
    bl_deltaR.push_back(lepton2.DeltaR(bottom2));
    bl_min_deltaR = *min_element(begin(bl_deltaR),end(bl_deltaR));
    bbll_deltaR = leptonlepton.DeltaR(bottombottom);
    bbll_deltaPhi = leptonlepton.DeltaPhi(bottombottom);
    // mT
    mT = sqrt(2*leptonlepton.Pt()*missing.E()*(1-cos(leptonlepton.Phi()-missing.Phi())));

    // lepton baseline selections
    if (ll_deltaR < 0.07 || ll_deltaR > 3.3) {return true;} step++; // pre-MVA selection 1
    if (leptonlepton.M() < 5 || leptonlepton.M() > 100) {return true;} step++; // pre-MVA selection 2
    // bottom baseline selections
    if (bb_deltaR > 5) {return true;} step++; // pre-MVA selection 3
    if (bottombottom.M() < 22) {return true;} step++; // pre-MVA selection 4

    // get TLorentzVector to the global variable for Minuit
    g_missing = missing;
    g_lepton1 = lepton1; g_lepton2 = lepton2;
    g_bottom1 = bottom1; g_bottom2 = bottom2;
    doubleHiggsAnalyser::GetHiggsness();
    doubleHiggsAnalyser::GetTopness();
    
    // MT2
    Mt2::LorentzTransverseVector vis_A(Mt2::TwoVector(leptonlepton.Px(), leptonlepton.Py()), leptonlepton.M());
    Mt2::LorentzTransverseVector vis_B(Mt2::TwoVector(bottombottom.Px(), bottombottom.Py()), bottombottom.M());
    Mt2::TwoVector pT_Miss(missing.Px(), missing.Py());
    //lester_MT2 = asymm_mt2_lester_bisect::get_mT2( // the simpliest way to calculate MT2
    //          leptonlepton.M(),leptonlepton.Px(),leptonlepton.Py(),
    //          bottombottom.M(),bottombottom.Px(),bottombottom.Px(),
    //          missing.Px(),missing.Py(),
    //          missing.M(),missing.M());
    basic_MT2_332_bbll = basic_mt2_332Calculator.mt2_332(vis_A, vis_B, pT_Miss, missing.M());
    ch_bisect_MT2_332 = ch_bisect_mt2_332Calculator.mt2_332(vis_A, vis_B, pT_Miss, missing.M());
    // MT2(bl+bl)
    bottomlepton1 = lepton1 + bottom2;
    bottomlepton2 = lepton2 + bottom1;
    Mt2::LorentzTransverseVector vis_A_blbl(Mt2::TwoVector(bottomlepton1.Px(), bottomlepton1.Py()), bottomlepton1.M());
    Mt2::LorentzTransverseVector vis_B_blbl(Mt2::TwoVector(bottomlepton2.Px(), bottomlepton2.Py()), bottomlepton2.M());
    basic_MT2_332_blbl = basic_mt2_332Calculator.mt2_332(vis_A_blbl, vis_B_blbl, pT_Miss, missing.M());
    // MT2(b)
    Mt2::LorentzTransverseVector vis_A_b(Mt2::TwoVector(bottom1.Px(), bottom1.Py()), bottom1.M());
    Mt2::LorentzTransverseVector vis_B_b(Mt2::TwoVector(bottom2.Px(), bottom2.Py()), bottom2.M());
    Mt2::TwoVector pT_Miss_b(leptonlepton.Px(), leptonlepton.Py());
    
    //lester_MT2_b = asymm_mt2_lester_bisect::get_mT2( // the simpliest way to calculate MT2
    //          bottom1.M(),bottom1.Px(),bottom1.Py(),
    //          bottom2.M(),bottom2.Px(),bottom2.Px(),
    //          leptonlepton.Px(),leptonlepton.Py(),
    //          leptonlepton.M(),leptonlepton.M());
    basic_MT2_332_b = basic_mt2_332Calculator.mt2_332(vis_A_b, vis_B_b, pT_Miss_b, leptonlepton.M());
    ch_bisect_MT2_332_b = ch_bisect_mt2_332Calculator.mt2_332(vis_A_b, vis_B_b, pT_Miss_b, leptonlepton.M());
    // MT2(l)
    Mt2::LorentzTransverseVector vis_A_l(Mt2::TwoVector(lepton1.Px(), lepton1.Py()), lepton1.M());
    Mt2::LorentzTransverseVector vis_B_l(Mt2::TwoVector(lepton2.Px(), lepton2.Py()), lepton2.M());
    Mt2::TwoVector pT_Miss_l(missing.Px(), missing.Py());
    //lester_MT2_l = asymm_mt2_lester_bisect::get_mT2( // the simpliest way to calculate MT2
    //          lepton1.M(),lepton1.Px(),lepton1.Py(),
    //          lepton2.M(),lepton2.Px(),lepton2.Px(),
    //          missing.Px(),missing.Py(),
    //          missing.M(),missing.M());
    basic_MT2_332_l = basic_mt2_332Calculator.mt2_332(vis_A_l, vis_B_l, pT_Miss_l, missing.M());
    ch_bisect_MT2_332_l = ch_bisect_mt2_332Calculator.mt2_332(vis_A_l, vis_B_l, pT_Miss_l, missing.M());
    
    tmva_bdtg_output = bdtg_reader->EvaluateMVA("BDTG");
    return true;
}

void doubleHiggsAnalyser::Loop() {
  TTree *t = 0;
  bool keep = false;
  t = del_tree;
  int nevents = t->GetEntries();
  event_size = 1;
  int proc = 0; int temp = 0;
  int nevent = 0;
  for (int iev = 0; iev < nevents; iev++) {
    t->GetEntry(iev);
    keep = doubleHiggsAnalyser::Analysis();
    if (keep) out_tree->Fill();
    nevent ++; temp = nevent*100/nevents;
    if ( temp != proc ){
        proc ++;
        cout << "###############################" << endl;
        cout << " proceeding : " << proc << " %" << endl;
        cout << "###############################" << endl;
    }
    n_events_tree->Fill();
  }
  cout << "n_events are " << nevents << endl;
}

void doubleHiggsAnalyser::Finalize() {
  out_tree->Write();
  out_file->Close();
}


int main(Int_t argc,Char_t** argv)
{
  //TFile *f = TFile::Open("/cms/scratch/jlee/hh/hh_1.root");
  //TTree *tree;
  //f->GetObject("Delphes", tree);
  TChain *tree = new TChain("Delphes");
  TString output_name;
  TString weight_file_path = "/cms/ldap_home/sunyoung/delPhys/src/delphys/analysis/bin/dataset_HH_SM/weights/DoubleHiggs_BDTG.weights.xml";
  //tree->Add("/xrootd/store/user/jlee/GluGluToHHTo2B2VTo2L2Nu_node_6_13TeV-madgraph-v2/delphes/*.root"); // ui
  //tree->Add("/xrootd/store/user/jlee/GluGluToHHTo2B2VTo2L2Nu_node_6_13TeV-madgraph-v2/delphes/*.root"); // B6 higgs
  //tree->Add("/xrootd/store/user/jlee/GluGluToHHTo2B2VTo2L2Nu_node_11_13TeV-madgraph-v2/delphes/*.root"); // B11 higgs
  //tree->Add("/xrootd/store/user/jlee/GluGluToHHTo2B2VTo2L2Nu_node_sm_13TeV-madgraph-v2/delphes/*.root"); // SM higgs
  for (int i = 1; i<200; i++){
  std::string filename = "/xrootd/store/user/seyang/Data/TopTagging/TT/TT_"+std::to_string(i)+".root";
  tree->Add(filename.c_str());
  }
  //output_name = "HH_B6.root";
  //output_name = "HH_B11.root";
  //output_name = "HH_SM.root";
  output_name = "TT_200.root";
    
  tree->SetBranchStatus("*",true);
    
  doubleHiggsAnalyser ana(tree, true);
  ana.Initiate(output_name); // usage : Initiate("outputfilename.root")
  ana.SetTMVA(weight_file_path);
  ana.Loop(); // Loop through the events and do doubleHIggsAnalyser::Analysis() per event.
  ana.Finalize(); // Write the tree and Close the file.
  return 0;
}
