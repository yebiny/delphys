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
    m_label(label),
    kIsDijet(is_dijet) {

  std::cout << "ctor begin" << std::endl;

  // Initialize properties
  bad_hard_gen_seen_ = false;
  num_passed_events_ = 0;

  //
  SetBranchAddress({"Vertex"}, /*drop=*/true);
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

  std::cout << "Passed / Total = "
            << num_passed_events_ << " / " << in_tree_->GetEntries()
            << std::endl;

  std::cout << "dtor end" << std::endl;
}

void QGJetsAnalyser::MakeBranch() {
  std::cout << "MakeBranch begin" << std::endl;

  ResetOnEachEvent();
  ResetOnEachJet();

  #define BRANCH_(name, suffix) out_tree_->Branch(#name, &m_##name, #name "/" #suffix);
  #define BRANCH_I(name) BRANCH_(name, I);
  #define BRANCH_F(name) BRANCH_(name, F);
  #define BRANCH_O(name) BRANCH_(name, O);

  #define BRANCH_A_(name, size, suffix) out_tree_->Branch(#name, &m_##name, #name"["#size"]/"#suffix);
  #define BRANCH_AF(name, size)  BRANCH_A_(name, size, F);

  #define BRANCH_VF(name) out_tree_->Branch(#name, "vector<Float_t>", &m_##name);
  #define BRANCH_VI(name) out_tree_->Branch(#name, "vector<Int_t>",   &m_##name);  
  #define BRANCH_VO(name) out_tree_->Branch(#name, "vector<Bool_t>",  &m_##name);  

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

  out_tree_->Branch("dau_p4", &m_dau_p4);
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
  m_event = 0; // nEvent
  m_num_jets = 0; // nJets
  m_num_good_jets = 0; // nGoodJets
  m_num_primary_vertices = 0; // nPriVtxs
  m_order = 0;
}


void QGJetsAnalyser::ResetOnEachJet() {
  // NOTE m_label do not need to be reset.

  m_pt = 0.0f;
  m_eta = 0.0f;
  m_phi = 0.0f;
  m_pt_dr_log = 0.0f;
  m_ptd = 0.0f; // jet energy sharing variable
  m_major_axis = 0.0f;
  m_minor_axis = 0.0f;

  // multiplicity
  m_chad_mult = 0; // charged hadrons
  m_nhad_mult = 0; // neutral hadrons
  m_electron_mult = 0;
  m_muon_mult = 0;
  m_photon_mult = 0;

  m_parton_id = 0;
  m_flavor_id = 0;
  m_flavor_algo_id = 0;
  m_flavor_phys_id = 0;

  m_num_dau = 0;
  m_matched = false;
  m_lepton_overlap = false;

  m_dau_p4.clear();
  m_dau_pt.clear();
  m_dau_deta.clear();
  m_dau_dphi.clear();
  m_dau_charge.clear();
  m_dau_pid.clear();
  m_dau_is_hadronic.clear();
  m_dau_is_track.clear();
  m_dau_eemfrac.clear();
  m_dau_ehadfrac.clear();

  #define FILL_ZERO(array, size) std::fill(array, array + size, 0.0);
  FILL_ZERO(m_image_electron_pt_33,   1089)
  FILL_ZERO(m_image_muon_pt_33,       1089)
  FILL_ZERO(m_image_photon_pt_33,     1089)
  FILL_ZERO(m_image_chad_pt_33,       1089)
  FILL_ZERO(m_image_nhad_pt_33,       1089)
  FILL_ZERO(m_image_electron_mult_33, 1089)
  FILL_ZERO(m_image_muon_mult_33,     1089)
  FILL_ZERO(m_image_photon_mult_33,   1089)
  FILL_ZERO(m_image_chad_mult_33,     1089)
  FILL_ZERO(m_image_nhad_mult_33,     1089)
}


Bool_t QGJetsAnalyser::SelectEvent() {
  if (kIsDijet) {
    return IsBalanced(jets_);
  } else {
    return PassZjets(jets_, muons_, electrons_);
  }
}



void QGJetsAnalyser::Analyse(Int_t entry) {
  m_event = entry; 
  m_num_jets = jets_->GetEntries();
  // TODO m_good_jets
  m_num_primary_vertices = vertices_ ? vertices_->GetEntries() : 1;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE [GenParticle] Find hard process partons
  //////////////////////////////////////////////////////////////////////////////
  std::vector<const GenParticle*> hard_partons;

  for (Int_t idx_gen = 0; idx_gen < particles_->GetEntries(); idx_gen++) {
    auto p = dynamic_cast<const GenParticle*>(particles_->At(idx_gen));

    if (p->Status != kHardProcessPartonStatus) continue;
    // Consider light quarks and gluons only
    if (std::abs(p->PID) > 5 and p->PID != kGluonPID) continue;

    hard_partons.push_back(p);

    if (hard_partons.size() == 2) break;
  }

  if (not bad_hard_gen_seen_ and (hard_partons.size() != 2)) {
    std::cout << "hard partons: " << hard_partons.size() << std::endl;
    bad_hard_gen_seen_ = true;
  }

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Jet
  //////////////////////////////////////////////////////////////////////////////

  for (Int_t idx_jet = 0; idx_jet < jets_->GetEntries(); idx_jet++) {
    ResetOnEachJet();

    if (kIsDijet and idx_jet >= 2) {
      // m_balanced = false;
      continue;
    }

    if (not kIsDijet and idx_jet >= 1) {
      // m_pass_zjets = false;
      continue;
    }


    auto jet = dynamic_cast<Jet*>(jets_->At(idx_jet));
    if (jet->PT < kMinJetPT) continue;
    if (std::fabs(jet->Eta) > kMaxJetEta) continue;

    ////////////////////////////////////////////////////////////////////////////
    // NOTE Match to hard process in Pythia8
    ////////////////////////////////////////////////////////////////////////////
    const GenParticle* matched_parton = nullptr;

    for (const auto & parton : hard_partons) {
      Float_t delta_r = jet->P4().DeltaR(parton->P4());

      if (delta_r < kDeltaRCut) {
        m_matched = true;
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
      if (delta_r < kDeltaRCut) {
        overlap = true;
        break;
      }
    }

    if (overlap) continue;


    //////////////////////////////////////////////////
    // NOTE check overlap with lepton
    ///////////////////////////////////////////////////
    for (Int_t idx_elec = 0; idx_elec < electrons_->GetEntries(); ++idx_elec) {
      auto electron = dynamic_cast<Electron*>(electrons_->At(idx_elec));
      Float_t delta_r = jet->P4().DeltaR(electron->P4());
      if (delta_r < kDeltaRCut) {
        m_lepton_overlap = true;
        break;
      }
    }

    if (not m_lepton_overlap) {
      for (Int_t idx_mu = 0; idx_mu < muons_->GetEntries(); ++idx_mu) {
        auto muon = dynamic_cast<Muon*>(muons_->At(idx_mu));
        Float_t delta_r = jet->P4().DeltaR(muon->P4());
        if (delta_r < kDeltaRCut) {
          m_lepton_overlap = true;
          break;
        }      
      }
    }

    /////////////////////////////////////////////
    // NOTE
    ///////////////////////////////////////////////
    m_pt = jet->PT;
    m_eta = jet->Eta;
    m_phi = jet->Phi;
    m_parton_id = matched_parton ? matched_parton->PID : 0;

    FillDaughters(jet);

    ////////////////////////////////////////
    // 
    //////////////////////////////////////////////
    std::tie(m_major_axis, m_minor_axis, std::ignore) = delphys::ComputeAxes(
        m_dau_deta, m_dau_dphi, m_dau_pt);

    // jet energy sharing variable
    Float_t sum_pt = std::accumulate(m_dau_pt.begin(), m_dau_pt.end(), 0.0f);
    Float_t sum_pt_squared = std::inner_product(
        m_dau_pt.begin(), m_dau_pt.end(),
        m_dau_pt.begin(), 0.0f);

    m_ptd = std::sqrt(sum_pt_squared) / sum_pt;

    MakeJetImage();

    // ExtractSatellites();

    out_tree_->Fill();

    m_order++;
  } // end loop over jets

}

/* Is the event balanced according to the criteria of pg 13 of
   http://cds.cern.ch/record/2256875/files/JME-16-003-pas.pdf */
Bool_t QGJetsAnalyser::IsBalanced(TClonesArray * jets) {
  if (jets->GetEntries() < 2) return false;

  auto jet1 = dynamic_cast<Jet*>(jets->At(0));
  auto jet2 = dynamic_cast<Jet*>(jets->At(1));

  // 2 jets of 30 GeV
  if (jet1->PT < 30.0) return false;
  if (jet2->PT < 30.0) return false;

  // that are back-to-back
  if (jet1->P4().DeltaPhi(jet2->P4()) < 2.5) return false;

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
Bool_t QGJetsAnalyser::PassZjets(TClonesArray * jets,
                                 TClonesArray * muons,
                                 TClonesArray * electrons) {

  Bool_t pass_zjets = false;

  if (jets->GetEntries() < 1) return false;

  // FIXME m_nGoodJets = 0;

  Float_t max_pt = 0.0f;
  Int_t idx_max_pt = -1;
  for (Int_t idx_jet = 0; idx_jet < jets_->GetEntries(); ++idx_jet) {
    auto jet = dynamic_cast<const Jet*>(jets->At(idx_jet));
    if (jet->PT < kMinJetPT) continue;
    if (std::fabs(jet->Eta) > kMaxJetEta) continue;

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
      m_num_dau++;
      deta = tower->Eta - jet->Eta;
      dphi = tower->P4().DeltaPhi(jet_p4);

      m_dau_p4.push_back(tower->P4());
      m_dau_pt.push_back(tower->ET);
      m_dau_deta.push_back(deta);
      m_dau_dphi.push_back(dphi);
      m_dau_charge.push_back(0);

      if (tower->Eem == 0.0) {
        m_nhad_mult++;
        m_dau_is_hadronic.push_back(true);
        m_dau_pid.push_back(0); // Neutral hadron
      } else if (tower->Ehad == 0.0) {
        m_photon_mult++;
        m_dau_is_hadronic.push_back(false); // leptonic
        m_dau_pid.push_back(kPhotonPID);
      } else {
          std::cout << "ERROR: Tower with Had " << tower->Ehad
                    << " and EM " << tower->Eem << " energy" << std::endl;
      }    

      m_dau_is_track.push_back(false);
      m_dau_eemfrac.push_back(tower->Eem / tower->E);
      m_dau_ehadfrac.push_back(tower->Ehad / tower->E);

    } else if (auto track = dynamic_cast<Track*>(daughter)) {
      m_num_dau++;
      deta = track->Eta - jet->Eta;
      dphi = track->P4().DeltaPhi(jet_p4);

      m_dau_p4.push_back(track->P4());
      m_dau_pt.push_back(track->PT);
      m_dau_deta.push_back(deta);
      m_dau_dphi.push_back(dphi);
      m_dau_charge.push_back(track->Charge);
      m_dau_pid.push_back(track->PID);
      m_dau_is_hadronic.push_back(false);
      m_dau_is_track.push_back(true);
      m_dau_eemfrac.push_back(0.0);
      m_dau_ehadfrac.push_back(0.0);


      Int_t abs_pid = std::abs(track->PID);
      if (abs_pid == kElectronPID) {
        m_electron_mult++;
      } else if (abs_pid == kMuonPID) {
        m_muon_mult++;
      } else {
        m_chad_mult++;
      }

    } else {
      std::cout << "BAD DAUGHTER! " << daughter << std::endl;
    }

  }

}

Int_t QGJetsAnalyser::Pixelate(Float_t deta,
                               Float_t dphi,
                               Float_t deta_max,
                               Float_t dphi_max, 
                               Int_t num_deta_bins,
                               Int_t num_dphi_bins) {
  Int_t eta_idx = static_cast<Int_t>(num_deta_bins * (deta + deta_max) / (2 * deta_max));
  Int_t phi_idx = static_cast<Int_t>(num_dphi_bins * (dphi + dphi_max) / (2 * dphi_max));
  Int_t idx = num_deta_bins * eta_idx + phi_idx;
  return idx;
}

void QGJetsAnalyser::MakeJetImage() {
  Float_t deta, dphi, pt;
  Int_t pixel;
  for (Int_t idx_dau = 0; idx_dau < m_num_dau; idx_dau++) {
    deta = m_dau_deta.at(idx_dau);
    dphi = m_dau_dphi.at(idx_dau);

    if (std::fabs(deta) >= kImageDeltaEtaMax) continue;
    if (std::fabs(dphi) >= kImageDeltaPhiMax) continue;

    pixel = Pixelate(deta, dphi);

    pt = m_dau_pt.at(idx_dau);
    Int_t abs_pid = std::abs(m_dau_pid.at(idx_dau));

    if (abs_pid == kElectronPID) {
      m_image_electron_pt_33[pixel] += pt;
      m_image_electron_mult_33[pixel] += 1.0;
    } else if (abs_pid == kMuonPID) {
      m_image_muon_pt_33[pixel] += pt;
      m_image_muon_mult_33[pixel] += 1.0;
    } else if (abs_pid == kPhotonPID) {
      m_image_photon_pt_33[pixel] += pt;
      m_image_photon_mult_33[pixel] += 1.0;
    } else if (abs_pid == 0) {
      m_image_nhad_pt_33[pixel] += pt;
      m_image_nhad_mult_33[pixel] += 1.0;
    } else {
      m_image_chad_pt_33[pixel] += pt;
      m_image_chad_mult_33[pixel] += 1.0;
    }
  } // end loop over daughters
}


void QGJetsAnalyser::Loop() {
  Int_t kNumEntries = in_tree_->GetEntries();
  Int_t kPrintFreq = kNumEntries / 20;
  TString kMsg = "[%d/%d (%.2f %)]";

  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);

    if (entry % kPrintFreq == 0) {

      std::cout << TString::Format(kMsg, entry + 1, kNumEntries,
                                   100 * Float_t(entry + 1) / kNumEntries)
                << std::endl;
        
    }


    if (SelectEvent()) {
      num_passed_events_++;
    } else {
      continue;
    }

    ResetOnEachEvent();
    Analyse(entry);
  }
}


int main(int argc, char* argv[]) {
  TString in_path(argv[1]);
  TString out_path(argv[2]);

  Bool_t is_dijet = in_path.Contains("qq") or in_path.Contains("gg");
  Int_t label = (in_path.Contains("qq") or in_path.Contains("zq")) ? 1 : 0;

  QGJetsAnalyser analyser(in_path, out_path, "jetAnalyser", is_dijet, label);
  analyser.Loop();

  return 0;
}
