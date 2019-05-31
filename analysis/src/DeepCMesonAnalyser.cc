#include "delphys/analysis/interface/DeepCMesonAnalyser.h"
#include "TLorentzVector.h"

DeepCMesonAnalyser::DeepCMesonAnalyser(const TString & in_path,
                                       const TString & out_path,
                                       const TString & out_tree_name):
    BaseAnalyser(in_path, out_path, out_tree_name),
    jet_label_(-1), 
    jet_num_d0dau_(0),
    jet_num_track_(0),
    num_pion_cand_(0),
    num_kaon_cand_(0){
    setBranchAddress({"Vertex"},true);
    makeBranch();
}    

DeepCMesonAnalyser::~DeepCMesonAnalyser() {
    out_tree_->Print();
    out_file_->Write();
    out_file_->Close();
}


void DeepCMesonAnalyser::makeBranch() {
    std::cout << "Make Branch Begin" << std::endl;
    
    resetOnEachJet();

    out_tree_->Branch("track_deta", "std::vector<float>", &track_deta_);
    out_tree_->Branch("track_dphi", "std::vector<float>", &track_dphi_);
    out_tree_->Branch("track_d0", "std::vector<float>", &track_d0_);
    out_tree_->Branch("track_dz", "std::vector<float>", &track_dz_);
    out_tree_->Branch("track_xd", "std::vector<float>", &track_xd_);
    out_tree_->Branch("track_yd", "std::vector<float>", &track_yd_);
    out_tree_->Branch("track_zd", "std::vector<float>", &track_zd_);
    out_tree_->Branch("track_errd0", "std::vector<float>", &track_errd0_);
    
    out_tree_->Branch("track_pt", "std::vector<float>", &track_pt_);
    out_tree_->Branch("track_eta", "std::vector<float>", &track_eta_);
    out_tree_->Branch("track_phi", "std::vector<float>", &track_phi_);
    out_tree_->Branch("track_mass", "std::vector<float>", &track_mass_);
    out_tree_->Branch("track_P4", "std::vector<TLorentzVector>", &track_P4_);
    
    out_tree_->Branch("track_charge", "std::vector<int>", &track_charge_);
    out_tree_->Branch("track_pId", "std::vector<int>", &track_pId_);
    out_tree_->Branch("track_l", "std::vector<float>", &track_l_);
    
    out_tree_->Branch("track_costompId", "std::vector<int>", &track_costompId_);
    out_tree_->Branch("mother_pId", "std::vector<int>", &mother_pId_);
    
    out_tree_->Branch("pticle_is_d0dau", "std::vector<int>", &pticle_is_d0dau_);
    out_tree_->Branch("pticle_is_pion", "std::vector<int>", &pticle_is_pion_);
    out_tree_->Branch("pticle_is_kaon", "std::vector<int>", &pticle_is_kaon_);

    out_tree_->Branch("charge_pion_cand", "std::vector<int>", &charge_pion_cand_);
    out_tree_->Branch("charge_kaon_cand", "std::vector<int>", &charge_kaon_cand_);
    out_tree_->Branch("num_pion_cand", &num_pion_cand_, "num_pion_cand/I");
    out_tree_->Branch("num_kaon_cand", &num_kaon_cand_, "num_kaon_cand/I");
   
    out_tree_->Branch("pticle_label", "std::vector<int>", &pticle_label_);
    out_tree_->Branch("jet_label", &jet_label_, "jet_label/I");
    
    out_tree_->Branch("jet_num_d0dau", &jet_num_d0dau_, "jet_num_d0dau/I");
    out_tree_->Branch("jet_num_track", &jet_num_track_, "jet_num_track/I");
    out_tree_->Branch("jet_track_comb", "TLorentzVector", &jet_track_comb_);

    std::cout << "Make Branch End" << std::endl;
}

void DeepCMesonAnalyser::resetOnEachJet() {
    track_P4_.clear();
    
    track_deta_.clear();
    track_dphi_.clear();
    track_charge_.clear();
    track_pId_.clear();
    
    track_pt_.clear();
    track_eta_.clear();
    track_phi_.clear();
    track_mass_.clear();
    
    track_d0_.clear();
    track_dz_.clear();
    track_l_.clear();
    track_xd_.clear();
    track_yd_.clear();
    track_zd_.clear();    
    track_errd0_.clear();
   
    track_costompId_.clear();
    mother_pId_.clear();

    pticle_is_d0dau_.clear();
    pticle_is_pion_.clear();
    pticle_is_kaon_.clear();
    
    charge_pion_cand_.clear();
    charge_kaon_cand_.clear();

    pticle_label_.clear();    
    
    jet_label_ = -1; 
    jet_num_d0dau_ = 0;
    jet_num_track_ = 0;
    num_pion_cand_ = 0;
    num_kaon_cand_ = 0;

    jet_track_comb_.SetPtEtaPhiM(0,0,0,0);

    Idx.clear();
}

void DeepCMesonAnalyser::analyse(Int_t entry) {
    for (Int_t idx_jet = 0; idx_jet < jets_->GetEntries(); idx_jet++) {
        resetOnEachJet();
        auto jet = dynamic_cast<Jet*>(jets_->At(idx_jet));
       
        Float_t deta, dphi;
        for (Int_t idx_dau = 0; idx_dau < jet -> Constituents.GetEntries(); idx_dau++) {
            TObject* daughter = jet->Constituents.At(idx_dau);
            TLorentzVector jet_p4 = jet->P4();
            
            // Jet selection  
            if (abs(jet->Eta) > 2.4) continue;
            if (jet->PT < 20) continue;

            // We need only track particles
            if (auto track = dynamic_cast<Track*>(daughter)) {
                // Save all track informations
                jet_num_track_ = jet_num_track_ + 1;
                
                deta = track->Eta - jet->Eta;
                dphi = track->P4().DeltaPhi(jet_p4);
                
                track_deta_.push_back(deta);
                track_dphi_.push_back(dphi);
                
                // for momentum
                track_pt_.push_back(track->PT);
                track_eta_.push_back(track->Eta);
                track_phi_.push_back(track->Phi);
                track_mass_.push_back(track->P4().M());
                
                track_P4_.push_back(track->P4());

                track_d0_.push_back(track->D0);
                track_dz_.push_back(track->DZ);
                track_xd_.push_back(track->Xd);
                track_yd_.push_back(track->Yd);
                track_zd_.push_back(track->Zd);
                track_errd0_.push_back(track->ErrorD0);

                track_charge_.push_back(track->Charge);
                track_pId_.push_back(track->PID);
                track_l_.push_back(track->L);
                
                if (abs(track->PID) == 11) {
                    track_costompId_.push_back(1);
                } else if (abs(track->PID) == 13) {
                    track_costompId_.push_back(2);
                } else {
                    if (track->PID == 0) {
                        track_costompId_.push_back(0);
                    } else {
                        track_costompId_.push_back(3);
                    }
                }            

                // Check mother particle 
                auto gen = dynamic_cast<const GenParticle*>(track->Particle.GetObject());
                if (gen->M1 != -1) {
                    auto mother1 = dynamic_cast<const GenParticle*>(particles_->At(gen->M1));
                    mother_pId_.push_back(mother1->PID);

                    // Check if particle is from D0
                    if (abs(mother1->PID) == d0_pId_ and gen->M2 == -1) {
                        pticle_is_d0dau_.push_back(1);
                        jet_num_d0dau_ = jet_num_d0dau_ + 1; 
                        // Check if particle pId is kaon or pion 
                        if (abs(track->PID) == pion_pId_) {
                            pticle_label_.push_back(1);
                            pticle_is_pion_.push_back(1);
                            pticle_is_kaon_.push_back(0);
                            charge_pion_cand_.push_back(track->Charge);
                            num_pion_cand_ = num_pion_cand_+1;
                        } else if (abs(track->PID) == kaon_pId_) {
                            pticle_label_.push_back(1);
                            pticle_is_pion_.push_back(0);
                            pticle_is_kaon_.push_back(1);
                            charge_kaon_cand_.push_back(track->Charge);
                            num_kaon_cand_ = num_kaon_cand_+1;
                        } else {
                            pticle_label_.push_back(0);
                            pticle_is_pion_.push_back(0);
                            pticle_is_kaon_.push_back(0);
                        } 
                    } else { 
                        pticle_is_d0dau_.push_back(0);
                        pticle_label_.push_back(0);
                        if (abs(track->PID) == pion_pId_) {
                            pticle_is_pion_.push_back(1);
                            pticle_is_kaon_.push_back(0);
                        } else if (abs(track->PID) == kaon_pId_) {
                            pticle_is_pion_.push_back(0);
                            pticle_is_kaon_.push_back(1); 
                        } else {
                            pticle_is_pion_.push_back(0);
                            pticle_is_kaon_.push_back(0); 
                        }   
                    }  
                } else {
                    pticle_is_d0dau_.push_back(0);
                    pticle_is_pion_.push_back(0);
                    pticle_is_kaon_.push_back(0); 
                    pticle_label_.push_back(0);
                }//mother1
            }else continue; //track
        }//daughter 
    
        // Jet labelling    
        jet_label_ = 0;
        if (jet_num_d0dau_ >= 2) {
            jet_label_ = 1;
            if (num_pion_cand_ == 1 and num_kaon_cand_ ==1 ) {
                jet_label_ = 2;
                if (charge_pion_cand_[0] * charge_kaon_cand_[0] == -1) {
                    jet_label_ = 3;
                }
            }
        }
        
        // Particle relabelling
        if (jet_label_ != 3){
            replace(pticle_label_.begin(), pticle_label_.end(), 1,0);
        }
        
        // Get Particle Candidate
        Int_t size = pticle_label_.size(); 
        for ( Int_t i=0; i < size; i++ ){
            if (pticle_label_[i] == 1){
                Idx.push_back(i);
            }
        }
        
        // Particle four momentum
        if (Idx.size() == 2){
            jet_track_comb_ = track_P4_[Idx[0]] + track_P4_[Idx[1]];
        }

        // Jet selection last 
        if (jet_num_track_ > 1) {
        out_tree_->Fill();
        }
    }//jet
}
void DeepCMesonAnalyser::loop() {
  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);
    analyse(entry);
  }
} 
