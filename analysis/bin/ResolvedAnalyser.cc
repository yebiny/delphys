#include "delphys/analysis/interface/ResolvedAnalyser.h"

#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"

#include <algorithm>
#include <numeric>
#include <iterator>
#include <memory>



ResolvedAnalyser::ResolvedAnalyser(const TString & in_path,
                                   const TString & out_path,
                                   const TString & out_tree_name)
    : BaseAnalyser(in_path, out_path, out_tree_name) {

  SetBranchAddress({"Vertex"}, /*drop=*/true);
  MakeBranch();

  TString base_name(gSystem->BaseName(in_path));
  if (base_name.First("TT") != -1) {
    label_ = 1;
  } else if (base_name.First("QCD") != -1) {
    label_ = 0;
  } else {
    std::cout << base_name << std::endl;
    std::abort();
  }

  std::cout << "ctor" << std::endl;
}


ResolvedAnalyser::~ResolvedAnalyser() {
  out_tree_->Print();
  out_file_->Write();
  out_file_->Close();
  std::cout << "dtor" << std::endl;
}

void ResolvedAnalyser::MakeBranch() {
  Reset();
  gInterpreter->GenerateDictionary("vector<vector<Int_t> >", "vector"); 
  gInterpreter->GenerateDictionary("vector<vector<TLorentzVector> >", "vector"); 

  // per event
  out_tree_->Branch("label", &label_, "label/I");

  // per jet : vector of 
  out_tree_->Branch("jet_p4", &jet_p4_);

  out_tree_->Branch("jet_num_chad", &jet_num_chad_);
  out_tree_->Branch("jet_num_nhad", &jet_num_nhad_);
  out_tree_->Branch("jet_num_electron", &jet_num_electron_);
  out_tree_->Branch("jet_num_muon", &jet_num_muon_);
  out_tree_->Branch("jet_num_photon", &jet_num_photon_);

  out_tree_->Branch("jet_major_axis", &jet_major_axis_);
  out_tree_->Branch("jet_minor_axis", &jet_minor_axis_);
  out_tree_->Branch("jet_ptd", &jet_ptd_);

  out_tree_->Branch("jet_b_tag", &jet_b_tag_);
  out_tree_->Branch("jet_b_dr_matching", &jet_b_dr_matching_);
  // out_tree_->Branch("jet_b_tracking", &jets_b_tracking_);

  // per constituents : vector of vector of
  out_tree_->Branch("constituent_p4", &constituent_p4_);
  out_tree_->Branch("constituent_type", &constituent_type_);

  return ;
}


void ResolvedAnalyser::Reset() {
  jet_p4_.clear();

  jet_num_chad_.clear();
  jet_num_nhad_.clear();
  jet_num_electron_.clear();
  jet_num_muon_.clear();
  jet_num_photon_.clear();

  jet_major_axis_.clear();
  jet_minor_axis_.clear();
  jet_ptd_.clear();

  jet_b_tag_.clear();
  jet_b_dr_matching_.clear();
  // jet_b_tracking_.clear();

  constituent_p4_.clear();
  constituent_type_.clear();
}


std::vector<const Jet*> ResolvedAnalyser::SelectJet() {
  std::vector<const Jet*> selected_jets;

  for (Int_t i = 0; i < jets_->GetEntries(); i++) {
    auto jet = dynamic_cast<const Jet*>(jets_->At(i));

    if(jet->PT < 45.0f) continue;
    if(std::fabs(jet->Eta) > 2.4f) continue;

    selected_jets.push_back(jet);
  }
  return selected_jets;
}


Bool_t ResolvedAnalyser::SelectEvent() {
  if (jets_->GetEntries() < 6) {
    return false;
  }

  // NOTE Select objects
  selected_jets_ = SelectJet();
  if(selected_jets_.size() < 6) {
    return false;
  }

  return kTRUE;
}


void ResolvedAnalyser::AnalyseEvent() {

  std::vector<GenParticle*> bottom_quarks;
  for (Int_t p_idx = 0; p_idx < particles_->GetEntries(); p_idx++) {
    auto p = dynamic_cast<GenParticle*>(particles_->At(p_idx));
    if ((std::abs(p->PID) == 5) and (p->Status == 23)) {
      bottom_quarks.push_back(p);
    }
  } 

  for (const auto & jet : selected_jets_) {
    TLorentzVector j_p4 = jet->P4();
    Float_t j_eta = jet->Eta;
    Float_t j_phi = jet->Phi;

    constituent_p4_.push_back(std::vector<TLorentzVector>());
    constituent_type_.push_back(std::vector<Int_t>());

    Int_t num_chad = 0, num_nhad = 0, num_electron = 0, num_muon = 0, num_photon = 0;
    std::vector<Double_t> cons_deta, cons_dphi, cons_pt;

    for (Int_t con_idx = 0; con_idx < jet->Constituents.GetEntries(); con_idx++) {
      auto constituent = jet->Constituents.At(con_idx);

      TLorentzVector c_p4;
      Int_t c_type = -1;

      if (auto track = dynamic_cast<Track*>(constituent)) {
        c_p4 = track->P4();

        if (std::abs(track->PID) == kElectronPID_) { 
          c_type = kElectronType_;
          num_electron++;
        }
        else if (std::abs(track->PID) == kMuonPID_) {
          c_type = kMuonType_;
          num_muon++;
        }
        else { 
          c_type = kChargedHadronType_;
          num_chad++;
        }

      } else if (auto tower = dynamic_cast<Tower*>(constituent)) {
        c_p4 = tower->P4();

        if (tower->Eem == 0.0) {
          c_type = kNeutralHadronType_;
          num_nhad++;
        }
        else if (tower->Ehad == 0.0) {
          c_type = kPhotonType_;
          num_photon++;
        }
        else
          std::cout << ":p" << std::endl;

      } else {
        // FIXME cerr
        std::cout << ":p" << std::endl;
      }

      constituent_p4_.back().push_back(c_p4);
      constituent_type_.back().push_back(c_type); // photon

      cons_deta.push_back(c_p4.Eta() - j_eta);
      cons_dphi.push_back(TVector2::Phi_mpi_pi(c_p4.Phi() - j_phi));
      cons_pt.push_back(c_p4.Pt());

    } // end loop over constituents

    jet_num_chad_.push_back(num_chad);
    jet_num_nhad_.push_back(num_nhad);
    jet_num_electron_.push_back(num_electron);
    jet_num_muon_.push_back(num_muon);
    jet_num_photon_.push_back(num_photon);

    Float_t major_axis, minor_axis;
    std::tie(major_axis, minor_axis) = delphys::ComputeAxes(cons_deta, cons_dphi, cons_pt);

    jet_major_axis_.push_back(major_axis);
    jet_minor_axis_.push_back(minor_axis);

    // jet energy sharing variable
    Float_t pt_sum = std::accumulate(cons_pt.begin(), cons_pt.end(), 0.0f);
    Float_t pt_square_sum = std::inner_product(cons_pt.begin(), cons_pt.end(),
                                               cons_pt.begin(), 0.0f); 
    Float_t ptd = std::sqrt(pt_square_sum) / pt_sum;
    jet_ptd_.push_back(ptd);


    // NOTE b-jet properties
    jet_b_tag_.push_back(jet->BTag);

    // b-matching
    Bool_t found_matched_b = false;
    for (const auto & b : bottom_quarks) {
      Double_t delta_r = b->P4().DeltaR(jet->P4());
      if (delta_r < kBMatchingDeltaR_) {
        found_matched_b = true;
        break;
      }
    }
    jet_b_dr_matching_.push_back(found_matched_b);


    // is_from_b = IsFromB(jet);

    // jet_b_tracking_.push_back(true);


  } // end loop over jets

  out_tree_->Fill();
}

void ResolvedAnalyser::Loop() {
  Int_t kNumEntries = in_tree_->GetEntries();
  Int_t kPrintFreq = kNumEntries / 100;
  TString kMessageTemplate = "[%d/%d (%.2f %)] # of passed events: %d (%.2f %)\n";

  Int_t num_passed = 0;

  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);
    Reset();

    if (entry % kPrintFreq == 0) {
      TString msg = TString::Format(
          kMessageTemplate,
          entry + 1,
          kNumEntries,
          100 * static_cast<Float_t>(entry + 1) / kNumEntries,
          num_passed,
          100 * static_cast<Float_t>(num_passed) / (entry + 1));

      std::cout << msg;
    }

    // Bool_t is_selected = SelectEvent()
    if (SelectEvent()) {
      num_passed++;
      AnalyseEvent();
    }
  }
}

int main(int argc, char* argv[]) {
  TString in_path(argv[1]);
  TString out_path(argv[2]);

  ResolvedAnalyser analyser(in_path, out_path, "test");
  analyser.Loop();

  return 0;
}
