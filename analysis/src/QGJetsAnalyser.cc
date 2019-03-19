#include "delphys/analysis/interface/QGJetsAnalyser.h"

#include "TLorentzVector.h"

#include <algorithm> // std::fill
#include <cmath> // std::abs, std::fabs, std::sqrt
#include <tuple> // tie, ignore
#include <iterator>
#include <numeric> // std::accumulate, std::inner_product


QGJetsAnalyser::QGJetsAnalyser(const TString & in_path,
                               const TString & out_path,
                               const TString & out_tree_name,
                               Bool_t is_dijet,
                               Int_t label) : 
    BaseAnalyser(in_path, out_path, out_tree_name),
    label_(label),
    kIsDijet(is_dijet) {

  std::cout << "ctor begin" << std::endl;

  // Initialize properties
  bad_hard_gen_seen_ = false;
  qgjets_stats_["num_passed_events"] = 0;
  qgjets_stats_["overlap_1st"] = 0;
  qgjets_stats_["overlap_2nd"] = 0;

  //
  setBranchAddress({"Vertex"}, /*drop=*/true);
  MakeBranch();

  //
  std::cout << "Input file: " << in_path << std::endl;
  if (kIsDijet) {
    std::cout << "Dijet --> IsBalanced" << std::endl;
  } else {
    std::cout << "Z+jets --> PassZjets" << std::endl;
  }

  std::cout << "ctor end" << std::endl;
}

QGJetsAnalyser::~QGJetsAnalyser() {
  std::cout << "dtor begin" << std::endl;

  out_tree_->Print();
  out_file_->Write();
  out_file_->Close();

  for (const auto & each : qgjets_stats_) {
    std::cout << ">>> " << each.first << ": " << each.second << std::endl;
  }

  std::cout << "dtor end" << std::endl;
}

void QGJetsAnalyser::MakeBranch() {
  std::cout << "MakeBranch begin" << std::endl;

  ResetOnEachEvent();
  ResetOnEachJet();

  #define BRANCH_(name, suffix) out_tree_->Branch(#name, & name##_, #name "/" #suffix);
  #define BRANCH_I(name) BRANCH_(name, I);
  #define BRANCH_F(name) BRANCH_(name, F);
  #define BRANCH_O(name) BRANCH_(name, O);

  #define BRANCH_A_(name, size, suffix) out_tree_->Branch(#name, & name##_, #name"["#size"]/"#suffix);
  #define BRANCH_AF(name, size)  BRANCH_A_(name, size, F);

  #define BRANCH_V_(name, type) out_tree_->Branch(#name, "vector<"#type">", & name##_ );
  #define BRANCH_VF(name) BRANCH_V_(name, Float_t);
  #define BRANCH_VI(name) BRANCH_V_(name, Int_t);
  #define BRANCH_VO(name) BRANCH_V_(name, Bool_t);

  BRANCH_I(event)
  BRANCH_I(num_jets)
  BRANCH_I(num_good_jets)
  BRANCH_I(num_primary_vertices)
  BRANCH_I(order)

  BRANCH_I(label);

  BRANCH_F(pt)
  BRANCH_F(eta)
  BRANCH_F(phi)
  BRANCH_F(pt_dr_log)
  BRANCH_F(ptd)
  BRANCH_F(major_axis)
  BRANCH_F(minor_axis)

  BRANCH_I(chad_mult)
  BRANCH_I(nhad_mult)
  BRANCH_I(electron_mult)
  BRANCH_I(muon_mult)
  BRANCH_I(photon_mult)

  BRANCH_I(parton_id)
  BRANCH_I(flavor_id)
  BRANCH_I(flavor_algo_id)
  BRANCH_I(flavor_phys_id)


  out_tree_->Branch("dau_p4", &dau_p4_);
  BRANCH_VF(dau_pt)
  BRANCH_VF(dau_deta)
  BRANCH_VF(dau_dphi)

  BRANCH_VI(dau_charge)
  BRANCH_VI(dau_pid)

  BRANCH_VO(dau_is_hadronic)
  BRANCH_VO(dau_is_track)

  BRANCH_VF(dau_eemfrac)
  BRANCH_VF(dau_ehadfrac)

  BRANCH_I(num_dau)
  BRANCH_O(matched)
  // BRANCH_O(balanced)
  // BRANCH_O(pass_zjets)
  BRANCH_O(lepton_overlap)

  // transverse momentum
  BRANCH_AF(image_chad_pt_33,     1089)
  BRANCH_AF(image_nhad_pt_33,     1089)
  BRANCH_AF(image_electron_pt_33, 1089)
  BRANCH_AF(image_muon_pt_33,     1089)
  BRANCH_AF(image_photon_pt_33,   1089)
  // multiplicity
  BRANCH_AF(image_chad_mult_33,     1089)
  BRANCH_AF(image_nhad_mult_33,     1089)
  BRANCH_AF(image_electron_mult_33, 1089)
  BRANCH_AF(image_muon_mult_33,     1089)
  BRANCH_AF(image_photon_mult_33,   1089)

  std::cout << "MakeBranch end" << std::endl;
}

void QGJetsAnalyser::ResetOnEachEvent() {
  event_ = 0; // nEvent
  num_jets_ = 0; // nJets
  num_good_jets_ = 0; // nGoodJets
  num_primary_vertices_ = 0; // nPriVtxs
  order_ = 0;
}


void QGJetsAnalyser::ResetOnEachJet() {
  // NOTE label do not need to be reset.

  pt_ = 0.0f;
  eta_ = 0.0f;
  phi_ = 0.0f;
  pt_dr_log_ = 0.0f;
  ptd_ = 0.0f; // jet energy sharing variable
  major_axis_ = 0.0f;
  minor_axis_ = 0.0f;

  // multiplicity
  chad_mult_ = 0; // charged hadrons
  nhad_mult_ = 0; // neutral hadrons
  electron_mult_ = 0;
  muon_mult_ = 0;
  photon_mult_ = 0;

  parton_id_ = 0;
  flavor_id_ = 0;
  flavor_algo_id_ = 0;
  flavor_phys_id_ = 0;

  num_dau_ = 0;
  matched_ = false;
  lepton_overlap_ = false;

  dau_p4_.clear();
  dau_pt_.clear();
  dau_deta_.clear();
  dau_dphi_.clear();
  dau_charge_.clear();
  dau_pid_.clear();
  dau_is_hadronic_.clear();
  dau_is_track_.clear();
  dau_eemfrac_.clear();
  dau_ehadfrac_.clear();

  #define FILL_ZERO(array, size) std::fill(array, array + size, 0.0);
  FILL_ZERO(image_electron_pt_33_,   1089)
  FILL_ZERO(image_muon_pt_33_,       1089)
  FILL_ZERO(image_photon_pt_33_,     1089)
  FILL_ZERO(image_chad_pt_33_,       1089)
  FILL_ZERO(image_nhad_pt_33_,       1089)
  FILL_ZERO(image_electron_mult_33_, 1089)
  FILL_ZERO(image_muon_mult_33_,     1089)
  FILL_ZERO(image_photon_mult_33_,   1089)
  FILL_ZERO(image_chad_mult_33_,     1089)
  FILL_ZERO(image_nhad_mult_33_,     1089)
}


Bool_t QGJetsAnalyser::SelectEvent() {
  if (kIsDijet) return IsBalanced(jets_);
  else          return PassZjets(jets_, muons_, electrons_);
}


void QGJetsAnalyser::Analyse(Int_t entry) {
  event_ = entry; 
  num_jets_ = jets_->GetEntries();
  // TODO m_good_jets
  num_primary_vertices_ = vertices_ ? vertices_->GetEntries() : 1;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE [GenParticle] Find hard process partons
  //////////////////////////////////////////////////////////////////////////////
  std::vector<const GenParticle*> hard_partons;

  for (Int_t idx_gen = 0; idx_gen < particles_->GetEntries(); idx_gen++) {
    auto p = dynamic_cast<const GenParticle*>(particles_->At(idx_gen));
    if (p->Status != kHardProcessPartonStatus) continue;
    // Consider light quarks and gluons only
    if (std::abs(p->PID) > 5 and p->PID != kGluonPID_) continue;
    hard_partons.push_back(p);
    if (hard_partons.size() == 2) break;
  }

  if (not bad_hard_gen_seen_ and (hard_partons.size() != 2)) {
    std::cout << "hard partons: " << hard_partons.size() << std::endl;
    bad_hard_gen_seen_ = true;
  }

  //////////////////////////////////////////////////////////////////////////////
  // NOTE [Jet]
  //////////////////////////////////////////////////////////////////////////////
  Int_t max_jet_idx = kIsDijet ? 2 : 1;
  for (Int_t idx_jet = 0; idx_jet < max_jet_idx; idx_jet++) {
    ResetOnEachJet();

    // if (kIsDijet and idx_jet >= 2) continue;
    // if ((not kIsDijet) and idx_jet >= 1) continue;

    auto jet = dynamic_cast<Jet*>(jets_->At(idx_jet));
    // IsBalanced, PassZjets
    // if (jet->PT < kMinJetPT) continue;


    ////////////////////////////////////////////////////////////////////////////
    // NOTE Match to hard process in Pythia8
    ////////////////////////////////////////////////////////////////////////////
    const GenParticle* matched_parton = nullptr;
    for (const auto & parton : hard_partons) {
      Float_t delta_r = jet->P4().DeltaR(parton->P4());

      if (delta_r < kDeltaRCut_) {
        matched_ = true;
        matched_parton = parton;
        break;
      }
    }

    ///////////////////////////////////////////////
    // NOTE check overlapping jets
    ///////////////////////////////////////////////////////
    Bool_t overlap = false;
    for (Int_t idx_other = 0; idx_other < jets_->GetEntries(); idx_other++) {
      if (idx_jet == idx_other) continue;
      auto other_jet = dynamic_cast<Jet*>(jets_->At(idx_other));
      Float_t delta_r = jet->P4().DeltaR(other_jet->P4()); 
      if (delta_r < kDeltaRCut_) {
        overlap = true;
        break;
      }
    }

    if (overlap) {
      if (idx_jet == 0)      qgjets_stats_["overlap_1st"]++;
      else if (idx_jet == 1) qgjets_stats_["overlap_2nd"]++;
      continue;
    }


    //////////////////////////////////////////////////
    // NOTE check overlap with lepton
    ///////////////////////////////////////////////////
    for (Int_t idx_elec = 0; idx_elec < electrons_->GetEntries(); ++idx_elec) {
      auto electron = dynamic_cast<Electron*>(electrons_->At(idx_elec));
      Float_t delta_r = jet->P4().DeltaR(electron->P4());
      if (delta_r < kDeltaRCut_) {
        lepton_overlap_ = true;
        break;
      }
    }

    if (not lepton_overlap_) {
      for (Int_t idx_mu = 0; idx_mu < muons_->GetEntries(); ++idx_mu) {
        auto muon = dynamic_cast<Muon*>(muons_->At(idx_mu));
        Float_t delta_r = jet->P4().DeltaR(muon->P4());
        if (delta_r < kDeltaRCut_) {
          lepton_overlap_ = true;
          break;
        }      
      }
    }

    /////////////////////////////////////////////
    // NOTE Fill tree
    ///////////////////////////////////////////////
    pt_ = jet->PT;
    eta_ = jet->Eta;
    phi_ = jet->Phi;

    parton_id_ = matched_parton ? matched_parton->PID : 0;

    flavor_id_ = jet->Flavor;
    flavor_algo_id_ = jet->FlavorAlgo;
    flavor_phys_id_ = jet->FlavorPhys;

    FillDaughters(jet);

    ////////////////////////////////////////
    // 
    //////////////////////////////////////////////
    std::tie(major_axis_, minor_axis_, std::ignore) = delphys::ComputeAxes(
        dau_deta_, dau_dphi_, dau_pt_);

    // jet energy sharing variable
    Float_t sum_pt = std::accumulate(dau_pt_.begin(), dau_pt_.end(), 0.0f);
    Float_t sum_pt_squared = std::inner_product(
        dau_pt_.begin(), dau_pt_.end(), dau_pt_.begin(), 0.0f);

    ptd_ = std::sqrt(sum_pt_squared) / sum_pt;

    MakeJetImage();

    // ExtractSatellites();

    order_++;

    out_tree_->Fill();
  } // end loop over jets

}

/* Is the event balanced according to the criteria of pg 13 of
   http://cds.cern.ch/record/2256875/files/JME-16-003-pas.pdf */
Bool_t QGJetsAnalyser::IsBalanced(TClonesArray* jets) {
  if (jets->GetEntries() < 2) return false;

  auto jet1 = dynamic_cast<Jet*>(jets->At(0));
  auto jet2 = dynamic_cast<Jet*>(jets->At(1));

  // 2 jets of 30 GeV
  if (jet1->PT < kMinJetPT_) return false;
  if (jet2->PT < kMinJetPT_) return false;

  /***************************************************************************** 
   * The eta cut is there to match the tracker range (or should be around that),
   * if we discard one of the jets because of it, then, yes, the event should be
   * discarded.
   ****************************************************************************/
  if (std::fabs(jet1->Eta) > kMaxJetEta_) return false;
  if (std::fabs(jet2->Eta) > kMaxJetEta_) return false;

  // that are back-to-back
  if (jet1->P4().DeltaPhi(jet2->P4()) < kDeltaPhiCut_) return false;

  // and any 3rd jet requires pt < 30% of the avg. of the 2 leading jets
  if (jets->GetEntries() > 2) {
    auto jet3 = dynamic_cast<Jet*>(jets->At(2));
    return (jet3->PT < 0.3 * (0.5 * (jet1->PT  + jet2->PT)));
  } else {
    return true;
  }
}

/* Does the event pass the Zjets criteria according to the criteria of pg 11-12
   of http://cds.cern.ch/record/2256875/files/JME-16-003-pas.pdf */
Bool_t QGJetsAnalyser::PassZjets(TClonesArray* jets,
                                 TClonesArray* muons,
                                 TClonesArray* electrons) {

  Bool_t pass_zjets = false;

  if (jets->GetEntries() < 1) return false;

  // FIXME m_nGoodJets = 0;

  Float_t max_pt = 0.0f;
  Int_t idx_max_pt = -1;
  for (Int_t idx_jet = 0; idx_jet < jets_->GetEntries(); ++idx_jet) {
    auto jet = dynamic_cast<const Jet*>(jets->At(idx_jet));
    if (jet->PT < kMinJetPT_) continue;
    if (std::fabs(jet->Eta) > kMaxJetEta_) continue;

    // FIXME m_num_good_jets++;
    if (max_pt < jet->PT) {
      max_pt = jet->PT;
      idx_max_pt = idx_jet;
    }
  }

  if (idx_max_pt < 0) return false;

  // check for Z event
  TLorentzVector dimuon;
  for (Int_t idx_mu = 0; idx_mu < muons->GetEntries(); idx_mu++) {

    auto mu = dynamic_cast<Muon*>(muons->At(idx_mu));

    if (mu->PT < 20.0f) continue;

    for (Int_t idx_other = idx_mu; idx_other < muons->GetEntries(); idx_other++) {
      auto mu2 = dynamic_cast<Muon*>(muons->At(idx_other));
      if (mu2->PT < 20.0f) continue;
      if (mu->Charge * mu2->Charge > 0) continue;

      // dimuon candidate
      TLorentzVector candidate = (mu->P4() + mu2->P4());
      if (candidate.M() < 70.0f or candidate.M() > 110.0f) continue;

      pass_zjets = true;
      dimuon = candidate;
      break;
    }
  }


  if (pass_zjets) {
    auto hardest_jet = dynamic_cast<const Jet *>(jets->At(idx_max_pt));

    // require them to be back to back
    if (dimuon.DeltaPhi(hardest_jet->P4()) < 2.1) pass_zjets = false;

    for (Int_t idx_jet = 0; idx_jet < jets->GetEntries(); ++idx_jet) {
      if (idx_jet == idx_max_pt) continue;
      auto jet = dynamic_cast<const Jet *>(jets->At(idx_jet));
      if (jet->PT > 0.3 * dimuon.Pt())
        pass_zjets = false;
    }
  }

  return pass_zjets;
}

void QGJetsAnalyser::FillDaughters(const Jet* jet) {
  TLorentzVector jet_p4 = jet->P4();

  Int_t num_daughters = jet->Constituents.GetEntries();
  Float_t deta, dphi;
  for (Int_t idx_dau = 0; idx_dau < num_daughters; idx_dau++) {
    TObject* daughter = jet->Constituents.At(idx_dau);

    if (auto tower = dynamic_cast<Tower*>(daughter)) {
      num_dau_++;
      deta = tower->Eta - jet->Eta;
      dphi = tower->P4().DeltaPhi(jet_p4);

      dau_p4_.push_back(tower->P4());
      dau_pt_.push_back(tower->ET);
      dau_deta_.push_back(deta);
      dau_dphi_.push_back(dphi);
      dau_charge_.push_back(0);

      if (tower->Eem == 0.0) {
        nhad_mult_++;
        dau_is_hadronic_.push_back(true);
        dau_pid_.push_back(0); // Neutral hadron
      } else if (tower->Ehad == 0.0) {
        photon_mult_++;
        dau_is_hadronic_.push_back(false); // leptonic
        dau_pid_.push_back(kPhotonPID_);
      } else {
          std::cout << "ERROR: Tower with Had " << tower->Ehad
                    << " and EM " << tower->Eem << " energy" << std::endl;
      }    

      dau_is_track_.push_back(false);
      dau_eemfrac_.push_back(tower->Eem / tower->E);
      dau_ehadfrac_.push_back(tower->Ehad / tower->E);

    } else if (auto track = dynamic_cast<Track*>(daughter)) {
      num_dau_++;
      deta = track->Eta - jet->Eta;
      dphi = track->P4().DeltaPhi(jet_p4);

      dau_p4_.push_back(track->P4());
      dau_pt_.push_back(track->PT);
      dau_deta_.push_back(deta);
      dau_dphi_.push_back(dphi);
      dau_charge_.push_back(track->Charge);
      dau_pid_.push_back(track->PID);
      dau_is_hadronic_.push_back(false);
      dau_is_track_.push_back(true);
      dau_eemfrac_.push_back(0.0);
      dau_ehadfrac_.push_back(0.0);

      Int_t abs_pid = std::abs(track->PID);
      if (abs_pid == kElectronPID_)  electron_mult_++;
      else if (abs_pid == kMuonPID_) muon_mult_++;
      else                           chad_mult_++;

    } else {
      std::cout << "BAD DAUGHTER! " << daughter << std::endl;
    }

  }

}

Int_t QGJetsAnalyser::Pixelate(Float_t deta,
                               Float_t dphi,
                               Int_t num_deta_bins,
                               Int_t num_dphi_bins) {
  Int_t eta_idx = Int_t(num_deta_bins * (deta + kImageDeltaEtaMax_) / (2 * kImageDeltaEtaMax_));
  Int_t phi_idx = Int_t(num_dphi_bins * (dphi + kImageDeltaPhiMax_) / (2 * kImageDeltaPhiMax_));
  Int_t idx = num_deta_bins * eta_idx + phi_idx;
  return idx;
}

void QGJetsAnalyser::MakeJetImage() {
  Float_t deta, dphi, pt;
  Int_t pixel;
  for (Int_t idx_dau = 0; idx_dau < num_dau_; idx_dau++) {
    deta = dau_deta_.at(idx_dau);
    dphi = dau_dphi_.at(idx_dau);

    if (std::fabs(deta) >= kImageDeltaEtaMax_) continue;
    if (std::fabs(dphi) >= kImageDeltaPhiMax_) continue;

    pixel = Pixelate(deta, dphi);

    pt = dau_pt_.at(idx_dau);
    Int_t abs_pid = std::abs(dau_pid_.at(idx_dau));

    if (abs_pid == kElectronPID_) {
      image_electron_pt_33_[pixel] += pt;
      image_electron_mult_33_[pixel] += 1.0f;
    } else if (abs_pid == kMuonPID_) {
      image_muon_pt_33_[pixel] += pt;
      image_muon_mult_33_[pixel] += 1.0f;
    } else if (abs_pid == kPhotonPID_) {
      image_photon_pt_33_[pixel] += pt;
      image_photon_mult_33_[pixel] += 1.0f;
    } else if (abs_pid == 0) {
      image_nhad_pt_33_[pixel] += pt;
      image_nhad_mult_33_[pixel] += 1.0f;
    } else {
      image_chad_pt_33_[pixel] += pt;
      image_chad_mult_33_[pixel] += 1.0f;
    }
  } // end loop over daughters
}


void QGJetsAnalyser::Loop() {
  const Int_t kNumEntries = in_tree_->GetEntries();
  const Int_t kPrintFreq = kNumEntries / 20;
  TString kMsg = "[%d/%d (%.2f %)]";

  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);

    if (entry % kPrintFreq == 0) {

      std::cout << TString::Format(kMsg, entry + 1, kNumEntries,
                                   100 * Float_t(entry + 1) / kNumEntries)
                << std::endl;
        
    }


    if (not SelectEvent()) continue;

    qgjets_stats_["num_passed_events"]++;
    ResetOnEachEvent();
    Analyse(entry);
  }
}
