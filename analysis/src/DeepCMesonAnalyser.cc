#include "delphys/analysis/interface/DeepCMesonAnalyser.h"
#include "TLorentzVector.h"

DeepCMesonAnalyser::DeepCMesonAnalyser(const TString & in_path,
                                       const TString & out_path,
                                       const TString & out_tree_name):
    BaseAnalyser(in_path, out_path, out_tree_name),
    jet_label_(-1), 
    jet_count_d0dau_(0),
    jet_count_track_(0),
    jet_count_pion_(0),
    jet_count_kaon_(0){
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

    out_tree_->Branch("pion_P4", "TLorentzVector", &pion_p4_);
    out_tree_->Branch("kaon_P4", "TLorentzVector", &kaon_p4_);
    out_tree_->Branch("d0cand_p4", "TLorentzVector", &d0cand_p4_);
    
    out_tree_->Branch("track_deta", "std::vector<float>", &track_deta_);
    out_tree_->Branch("track_dphi", "std::vector<float>", &track_dphi_);
   
    out_tree_->Branch("track_pt", "std::vector<float>", &track_pt_);
    out_tree_->Branch("track_eta", "std::vector<float>", &track_eta_);
    out_tree_->Branch("track_phi", "std::vector<float>", &track_phi_);
    out_tree_->Branch("track_mass", "std::vector<float>", &track_mass_);
    
    out_tree_->Branch("track_d0", "std::vector<float>", &track_d0_);
    out_tree_->Branch("track_dz", "std::vector<float>", &track_dz_);
    out_tree_->Branch("track_xd", "std::vector<float>", &track_xd_);
    out_tree_->Branch("track_yd", "std::vector<float>", &track_yd_);
    out_tree_->Branch("track_zd", "std::vector<float>", &track_zd_);
    out_tree_->Branch("track_errd0", "std::vector<float>", &track_errd0_);
    out_tree_->Branch("track_l", "std::vector<float>", &track_l_);
    
    out_tree_->Branch("track_pId", "std::vector<int>", &track_pId_);
    out_tree_->Branch("track_costompId", "std::vector<int>", &track_costompId_);
    out_tree_->Branch("track_charge", "std::vector<int>", &track_charge_);
    
    out_tree_->Branch("pticle_is_d0dau", "std::vector<int>", &pticle_is_d0dau_);
    out_tree_->Branch("pticle_is_pion", "std::vector<int>", &pticle_is_pion_);
    out_tree_->Branch("pticle_is_kaon", "std::vector<int>", &pticle_is_kaon_);
    out_tree_->Branch("pticle_label", "std::vector<int>", &pticle_label_);

    out_tree_->Branch("pion_charge", "std::vector<int>", &pion_charge_);
    out_tree_->Branch("kaon_charge", "std::vector<int>", &kaon_charge_);
    
    out_tree_->Branch("jet_label", &jet_label_, "jet_label/I");
    out_tree_->Branch("jet_count_d0dau", &jet_count_d0dau_, "jet_count_d0dau/I");
    out_tree_->Branch("jet_count_track", &jet_count_track_, "jet_count_track/I");
    out_tree_->Branch("jet_count_pion", &jet_count_pion_, "jet_count_pion/I");
    out_tree_->Branch("jet_count_kaon", &jet_count_kaon_, "jet_count_kaon/I");

    std::cout << "Make Branch End" << std::endl;
}

void DeepCMesonAnalyser::resetOnEachJet() {
    pion_p4_.SetPtEtaPhiM(0,0,0,0);
    kaon_p4_.SetPtEtaPhiM(0,0,0,0);
    d0cand_p4_.SetPtEtaPhiM(0,0,0,0);
    
    track_deta_.clear();
    track_dphi_.clear();
    
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
   
    track_charge_.clear();
    track_pId_.clear();
    track_costompId_.clear();
    
    mother_pId_.clear();
    mother_.clear();

    pticle_is_d0dau_.clear();
    pticle_is_pion_.clear();
    pticle_is_kaon_.clear();
    pticle_label_.clear();    
    
    pion_charge_.clear();
    pion_pt_.clear();
    pion_eta_.clear();
    pion_phi_.clear();
    pion_mass_.clear();
    
    kaon_charge_.clear();
    kaon_pt_.clear();
    kaon_eta_.clear();
    kaon_phi_.clear();
    kaon_mass_.clear();
    
    jet_label_ = -1; 
    jet_count_d0dau_ = 0;
    jet_count_track_ = 0;
    jet_count_pion_ = 0;
    jet_count_kaon_ = 0;

    Idx.clear();
}

void DeepCMesonAnalyser::analyse(Int_t entry) {
    for (Int_t idx_jet = 0; idx_jet < jets_->GetEntries(); idx_jet++) {
        resetOnEachJet();
        auto jet = dynamic_cast<Jet*>(jets_->At(idx_jet));
       
        TLorentzVector jet_p4 = jet->P4();
        Float_t deta, dphi;
        for (Int_t idx_dau = 0; idx_dau < jet -> Constituents.GetEntries(); idx_dau++) {
            TObject* daughter = jet->Constituents.At(idx_dau);
            
            // Jet selection  
            if (abs(jet->Eta) > 2.4) continue;
            if (jet->PT < 20) continue;

            // We need only track particles
            if (auto track = dynamic_cast<Track*>(daughter)) {
                // Save all track informations
                jet_count_track_ = jet_count_track_ + 1;
                
                deta = track->Eta - jet->Eta;
                dphi = track->P4().DeltaPhi(jet_p4);
                
                track_deta_.push_back(deta);
                track_dphi_.push_back(dphi);
                
                //std::cout << track->P4().M() << std::endl;
                
                track_d0_.push_back(track->D0);
                track_dz_.push_back(track->DZ);
                track_xd_.push_back(track->Xd);
                track_yd_.push_back(track->Yd);
                track_zd_.push_back(track->Zd);
                track_errd0_.push_back(track->ErrorD0);

                track_charge_.push_back(track->Charge);
                track_pId_.push_back(track->PID);
                track_pId_.push_back(track->PID);
                track_l_.push_back(track->L);
                
                if (abs(track->PID) == 11) {
                    track_costompId_.push_back(1);
                } else if (abs(track->PID) == 13) {
                    track_costompId_.push_back(2);
                } else if (track->PID == 0) {
                    track_costompId_.push_back(0);
                } else {track_costompId_.push_back(3);}
                          
                // Check mother particle 
                auto gen = dynamic_cast<const GenParticle*>(track->Particle.GetObject());
                if (gen->M1 != -1) {
                    auto mother1 = dynamic_cast<const GenParticle*>(particles_->At(gen->M1));
                    mother_pId_.push_back(mother1->PID);
                    //mother_.push_back(mother1->Mass);
                    mother_.push_back(mother1->M1);

                    // Check if particle is from D0
                    if (abs(mother1->PID) == d0_pId_ and gen->M2 == -1) {
                        pticle_is_d0dau_.push_back(1);
                        jet_count_d0dau_ = jet_count_d0dau_ + 1; 
                        // Check if particle pId is kaon or pion 
                        if (abs(track->PID) == pion_pId_) {
                            pticle_label_.push_back(1);
                            pticle_is_pion_.push_back(1);
                            pticle_is_kaon_.push_back(0);
                            
                            // for momentum
                            pion_pt_.push_back(track->PT);
                            pion_eta_.push_back(track->Eta);
                            pion_phi_.push_back(track->Phi);
                            pion_mass_.push_back(0.1396);
                            pion_charge_.push_back(track->Charge);
                            jet_count_pion_ = jet_count_pion_+1;
                        
                        } else if (abs(track->PID) == kaon_pId_) {
                            pticle_label_.push_back(1);
                            pticle_is_pion_.push_back(0);
                            pticle_is_kaon_.push_back(1);
                            
                            // for momentum
                            kaon_pt_.push_back(track->PT);
                            kaon_eta_.push_back(track->Eta);
                            kaon_phi_.push_back(track->Phi);
                            kaon_mass_.push_back(0.4937);
                            kaon_charge_.push_back(track->Charge);
                            jet_count_kaon_ = jet_count_kaon_+1;
                        
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
        if (jet_count_d0dau_ >= 2) {
            jet_label_ = 1;
            if (jet_count_pion_ == 1 and jet_count_kaon_ ==1){
                jet_label_ = 2;
                if ( pion_charge_[0] * kaon_charge_[0] == -1) {
                    jet_label_ = 3;
                }
            }
        }
        
        // Particle relabelling
        if (jet_label_ != 3) {
            replace(pticle_label_.begin(), pticle_label_.end(), 1, 0);
        }
        
        // Get Particle Candidate
        Int_t size = pticle_label_.size(); 
        for ( Int_t i=0; i < size; i++ ) {
            if (pticle_label_[i] == 1) {
                Idx.push_back(i);
            }
        }

        // D0 candidate four momentum
        if (Idx.size() == 2) {
            
            pion_p4_.SetPtEtaPhiM(pion_pt_[0], pion_eta_[0], pion_phi_[0], pion_mass_[0]); 
            kaon_p4_.SetPtEtaPhiM(kaon_pt_[0], kaon_eta_[0], kaon_phi_[0], kaon_mass_[0]);
            d0cand_p4_ = kaon_p4_ + pion_p4_;
            
            // Set mass window
            if (abs(d0cand_p4_.M() - d0_m_) > 0.2) {
                d0cand_p4_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
            }
            // check vars 
            //else {
            //    std::cout << "-----------------" << std::endl;
            //    std::cout << mother_pId_[Idx[0]] <<","<< mother_pId_[Idx[1]] << std::endl;
            //    std::cout << mother_[Idx[0]] <<","<< mother_[Idx[1]] << std::endl;
            //    std::cout << "comb mass: " << d0cand_p4_.M()<< ",d0 mass:" << d0_m_ << std::endl;
            //    std::cout << kaon_pt_[0]<<"," << kaon_eta_[0]<< "," <<kaon_phi_[0]<< "," <<kaon_mass_[0]  <<std::endl;
            //    std::cout << pion_pt_[0]<< "," <<pion_eta_[0]<< "," <<pion_phi_[0]<< "," <<pion_mass_[0]  <<std::endl;
            //}
        }
            
        // Last check
        if (jet_count_track_ < 2) continue;
        out_tree_->Fill();

    }//jet
}

void DeepCMesonAnalyser::loop() {
  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);
    analyse(entry);
  }
} 
