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

#include "delphys/analysis/interface/doubleHiggsAnalyser.h"
#include "classes/DelphesClasses.h"
#include "delphys/external/interface/lester_mt2_bisect.h"

#include "delphys/oxbridgekinetics/src/Mt2/Basic_Mt2_332_Calculator.h"
#include "delphys/oxbridgekinetics/src/Mt2/ChengHanBisect_Mt2_332_Calculator.h"

using namespace std;

// return mother particle(GenParticle) of particle 'p' among 'particles'.
GenParticle* getMother(TClonesArray* particles, GenParticle* p){
  if (p->M1==-1) return nullptr;
  auto mom = static_cast<GenParticle *>(particles->At(p->M1));
  while (mom->PID == p->PID){
    mom = static_cast<GenParticle *>(particles->At(mom->M1));
    if (mom->M1==-1) return nullptr;
  }
  return mom;
}

// return the PID of the mother particle of the particle at 'ip' among 'particles'
int isFrom(TClonesArray* particles, int ip){
   auto p = static_cast<GenParticle *>(particles->At(ip));
   // check if it's from Higgs
   auto mom = getMother(particles, p); 
   if (mom==nullptr) return 0;
   auto grmom = getMother(particles, mom);
   if (grmom==nullptr) return 0;
   return grmom->PID;
}

void doubleHiggsAnalyser::MakeOutputBranch(TTree *tree) {
  // MT2 variables
  tree->Branch("lester_MT2",&lester_MT2,"lester_MT2/F");
  tree->Branch("basic_MT2_332",&basic_MT2_332,"basic_MT2_332/F");
  tree->Branch("ch_bisect_MT2_332",&ch_bisect_MT2_332,"ch_bisect_MT2_332/F");
  // MT2(b) variables
  tree->Branch("lester_MT2_b",&lester_MT2_b,"lester_MT2_b/F");
  tree->Branch("basic_MT2_332_b",&basic_MT2_332_b,"basic_MT2_332_b/F");
  tree->Branch("ch_bisect_MT2_332_b",&ch_bisect_MT2_332_b,"ch_bisect_MT2_332_b/F");
  // MT2(l) variables
  tree->Branch("lester_MT2_l",&lester_MT2_l,"lester_MT2_l/F");
  tree->Branch("basic_MT2_332_l",&basic_MT2_332_l,"basic_MT2_332_l/F");
  tree->Branch("ch_bisect_MT2_332_l",&ch_bisect_MT2_332_l,"ch_bisect_MT2_332_l/F");

  tree->Branch("lepton1_pt",&lepton1_pt,"lepton1_pt/F");
  tree->Branch("lepton2_pt",&lepton2_pt,"lepton1_pt/F");
  tree->Branch("missing_et",&missing_et,"missing_et/F");
  tree->Branch("mt",&mt,"mt/F");

  // truth matching variables
  tree->Branch("fromHiggs",&fromHiggs,"fromHiggs/I");
  tree->Branch("fromTop",&fromTop,"fromTop/I");
  tree->Branch("fromZ",&fromZ,"fromZ/I");

  tree->Branch("lepton1MotherPID",&from1,"from1/I");
  tree->Branch("lepton2MotherPID",&from2,"from2/I");
}

void doubleHiggsAnalyser::SetOutput(TString output_file_name) {
  out_file = TFile::Open(output_file_name,"RECREATE");
  out_tree = new TTree("events","events"); 
}

void doubleHiggsAnalyser::SetBranchAddress() {
  //del_tree->SetBranchAddress("Particle",&particles);
  del_tree->SetBranchAddress("Particle",&particles);
  del_tree->SetBranchAddress("MissingET",&missings);
  del_tree->SetBranchAddress("Jet",&jets);
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
  lester_MT2 = -99;
  basic_MT2_332 = -99;
  ch_bisect_MT2_332 = -99;
  lepton1_pt = -99;
  lepton2_pt = -99;
  missing_et = -99;
  mt = -99;
  fromHiggs = 0;
  fromTop = 0;
  fromZ = 0;
  from1 = 0;
  from2 = 0;
  
  muons.clear();
  bottoms.clear();
}

bool doubleHiggsAnalyser::Analysis() {
  //map<Float_t, int, greater<Float_t>> muons : map of <pt,index>:<K,V> of muon sorted by pt.
  //base selections : MissingET > 20, pT(lepton) > 20, deltaR(ll) < 1.0, m(ll) < 65, deltaR(bb) < 1.3, 95 < m(bb) < 140
    
    // Missing ET
    auto m = static_cast<const MissingET *>(missings->At(0)); // There is always one MET object.
    if (m->MET<20) return false;
    missing.SetPtEtaPhiM(m->MET,0,m->Phi,0);
    missing_et = m->MET;
 
    // b jets
    for (int ij = 0; ij < jets->GetEntries(); ij++){
      auto jet = static_cast<const Jet *>(jets->At(ij));
      if (abs(jet->Flavor) != doubleHiggsAnalyser::Bottom_PID) continue;
      bottoms.insert(make_pair(jet->PT,ij));
    }
   
    if (bottoms.size()<2) {
      return false;
    }


    bottom_iter = bottoms.begin();
    
    auto bot1 = static_cast<const Jet *>(jets->At(bottom_iter->second));
    bottom1.SetPtEtaPhiM(bot1->PT,bot1->Eta,bot1->Phi,bot1->Mass); 
    bottom_iter ++;
    
    auto bot2 = static_cast<const Jet *>(jets->At(bottom_iter->second));
    bottom2.SetPtEtaPhiM(bot2->PT,bot2->Eta,bot2->Phi,bot2->Mass);
    bottombottom = bottom1+bottom2;
    
    if (bottombottom.M()<95 || bottombottom.M()>140 || bottom1.DeltaR(bottom2)>1.3) {
        return false;
    }
    
    // leptons (muons)
    for (int ip = 0; ip < particles->GetEntries(); ip++){
      auto p = static_cast<const GenParticle *>(particles->At(ip));
      if (abs(p->PID)!=doubleHiggsAnalyser::Muon_PID || p->PT<20 || fabs(p->Eta)>2.5 || p->M1==-1) continue;
      //if (abs(p->PID)!=doubleHiggsAnalyser::Tau_PID || abs(p->PID)!=doubleHiggsAnalyser::Electron_PID || abs(p->PID)!=doubleHiggsAnalyser::Muon_PID || p->PT<20 || fabs(p->Eta)>2.5 || p->M1==-1) continue;
      muons.insert(make_pair(p->PT,ip));
    }
  
    if (muons.size()<2) {
      return false;
    }
    
    muon_iter = muons.begin();
    
    auto mu1 = static_cast<const GenParticle *>(particles->At(muon_iter->second));
    int from1 = isFrom(particles, muon_iter->second);
    from1 = abs(from1);
    if (from1==doubleHiggsAnalyser::Higgs_PID) fromHiggs++;
    if (from1==doubleHiggsAnalyser::Top_PID) fromTop++;
    if (from1==doubleHiggsAnalyser::Z_PID) fromZ++;
    muon1.SetPtEtaPhiM(mu1->PT,mu1->Eta,mu1->Phi,mu1->Mass);
    ++muon_iter;
    auto mu2 = static_cast<const GenParticle *>(particles->At(muon_iter->second));
    int from2 = isFrom(particles, muon_iter->second);
    from2 = abs(from2);
    if (from2==doubleHiggsAnalyser::Higgs_PID) fromHiggs++;
    if (from2==doubleHiggsAnalyser::Top_PID) fromTop++;
    if (from2==doubleHiggsAnalyser::Z_PID) fromZ++;
    muon2.SetPtEtaPhiM(mu2->PT,mu2->Eta,mu2->Phi,mu2->Mass);
    
    lepton1_pt = mu1->PT;
    lepton2_pt = mu2->PT;
    
    auto muonmuon = muon1+muon2;
    if (muonmuon.M() > 65 || muon1.DeltaR(muon2) > 1) {
        return false;
    }
    mt = muonmuon.Mt();

    // MT2
    Mt2::LorentzTransverseVector vis_A(Mt2::TwoVector(muonmuon.Px(), muonmuon.Py()), muonmuon.M());
    Mt2::LorentzTransverseVector vis_B(Mt2::TwoVector(bottombottom.Px(), bottombottom.Py()), bottombottom.M());
    Mt2::TwoVector pT_Miss(missing.Px(), missing.Py());
    
    lester_MT2 = asymm_mt2_lester_bisect::get_mT2( // the simpliest way to calculate MT2
              muonmuon.M(),muonmuon.Px(),muonmuon.Py(),
              bottombottom.M(),bottombottom.Px(),bottombottom.Px(),
              missing.Px(),missing.Py(),
              missing.M(),missing.M());
    basic_MT2_332 = basic_mt2_332Calculator.mt2_332(vis_A, vis_B, pT_Miss, missing.M());
    ch_bisect_MT2_332 = ch_bisect_mt2_332Calculator.mt2_332(vis_A, vis_B, pT_Miss, missing.M());
    // MT2(b)
    Mt2::LorentzTransverseVector vis_A_b(Mt2::TwoVector(bottom1.Px(), bottom1.Py()), bottom1.M());
    Mt2::LorentzTransverseVector vis_B_b(Mt2::TwoVector(bottom2.Px(), bottom2.Py()), bottom2.M());
    Mt2::TwoVector pT_Miss_b(muonmuon.Px(), muonmuon.Py());
    
    lester_MT2_b = asymm_mt2_lester_bisect::get_mT2( // the simpliest way to calculate MT2
              bottom1.M(),bottom1.Px(),bottom1.Py(),
              bottom2.M(),bottom2.Px(),bottom2.Px(),
              muonmuon.Px(),muonmuon.Py(),
              muonmuon.M(),muonmuon.M());
    basic_MT2_332_b = basic_mt2_332Calculator.mt2_332(vis_A_b, vis_B_b, pT_Miss_b, muonmuon.M());
    ch_bisect_MT2_332_b = ch_bisect_mt2_332Calculator.mt2_332(vis_A_b, vis_B_b, pT_Miss_b, muonmuon.M());
    // MT2(l)
    Mt2::LorentzTransverseVector vis_A_l(Mt2::TwoVector(muon1.Px(), muon1.Py()), muon1.M());
    Mt2::LorentzTransverseVector vis_B_l(Mt2::TwoVector(muon2.Px(), muon2.Py()), muon2.M());
    Mt2::TwoVector pT_Miss_l(missing.Px(), missing.Py());
    
    lester_MT2_l = asymm_mt2_lester_bisect::get_mT2( // the simpliest way to calculate MT2
              muon1.M(),muon1.Px(),muon1.Py(),
              muon2.M(),muon2.Px(),muon2.Px(),
              missing.Px(),missing.Py(),
              missing.M(),missing.M());
    basic_MT2_332_l = basic_mt2_332Calculator.mt2_332(vis_A_l, vis_B_l, pT_Miss_l, missing.M());
    ch_bisect_MT2_332_l = ch_bisect_mt2_332Calculator.mt2_332(vis_A_l, vis_B_l, pT_Miss_l, missing.M());
    //if (ch_bisect_MT2_332 < 0.1) return false;
    if (basic_MT2_332==-99) return false;
  
  return true;
}

void doubleHiggsAnalyser::Loop() {
  TTree *t = 0;
  bool keep = false;
  t = del_tree;

  for (int iev = 0; iev < t->GetEntries(); iev++) {
    t->GetEntry(iev);
    keep = doubleHiggsAnalyser::Analysis();
    if (keep) out_tree->Fill();
    doubleHiggsAnalyser::ResetVariables();
  }
}

void doubleHiggsAnalyser::Finalize() {
  out_tree->Write();
  out_file->Close();
}

int main(Int_t argc, Char_t** argv)
{
    //TFile *f = TFile::Open("/cms/scratch/jlee/hh/hh_1.root");
    //TTree *tree;
    //f->GetObject("Delphes", tree);
    TChain *tree = new TChain("Delphes");

    //tree->Add("/cms/scratch/jlee/hh/*.root"); // ui
    tree->Add("/home/scratch/sunyoung/data/nanoAOD/hh/*.root"); // gate

    tree->SetBranchStatus("*",true);
    
    // for delphes format analysis
    doubleHiggsAnalyser ana(tree, true);
    // for nanoAOD format analysis
    // double HiggsAnalyser ana(tree, true);
    ana.Initiate("test.root"); // usage : Initiate("outputfilename.root")
    ana.Loop(); // Loop through the events and do doubleHIggsAnalyser::Analysis() per event.
    ana.Finalize(); // Write the tree and Close the file.

  return 0;
}

