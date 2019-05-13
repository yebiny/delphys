#include "delphys/analysis/interface/DeepCMesonAnalyser.h"
#include "TLorentzVector.h"

DeepCMesonAnalyser::DeepCMesonAnalyser(const TString & in_path,
                                       const TString & out_path,
                                       const TString & out_tree_name):
    BaseAnalyser(in_path, out_path, out_tree_name) {
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
    out_tree_->Branch("jet_label_jpsi", &jet_label_jpsi_, "jet_label_jpsi/I");
    out_tree_->Branch("jet_label_d0", &jet_label_d0_, "jet_label_d0/I");
    
    out_tree_->Branch("num_d0Dau", &num_d0Dau_, "num_d0Dau/I");
    out_tree_->Branch("num_jpsiDau", &num_jpsiDau_, "num_jpsiDau/I");

    out_tree_->Branch("track_pt", "std::vector<float>", &track_pt_);
    out_tree_->Branch("track_deta", "std::vector<float>", &track_deta_);
    out_tree_->Branch("track_dphi", "std::vector<float>", &track_dphi_);
    out_tree_->Branch("track_charge", "std::vector<int>", &track_charge_);
    out_tree_->Branch("track_pid", "std::vector<int>", &track_pid_);

    out_tree_->Branch("track_d0", "std::vector<float>", &track_d0_);
    out_tree_->Branch("track_dz", "std::vector<float>", &track_dz_);
    out_tree_->Branch("track_l", "std::vector<float>", &track_l_);
    out_tree_->Branch("track_xd", "std::vector<float>", &track_xd_);
    out_tree_->Branch("track_yd", "std::vector<float>", &track_yd_);
    out_tree_->Branch("track_zd", "std::vector<float>", &track_zd_);

    out_tree_->Branch("tower_pt", "std::vector<float>", &tower_pt_);
    out_tree_->Branch("tower_deta", "std::vector<float>", &tower_deta_);
    out_tree_->Branch("tower_dphi", "std::vector<float>", &tower_dphi_);
    out_tree_->Branch("tower_charge", "std::vector<int>", &tower_charge_);

    out_tree_->Branch("Dau_jpsi", "std::vector<int>", &Dau_jpsi_);
    out_tree_->Branch("Dau_d0", "std::vector<int>", &Dau_d0_);
    out_tree_->Branch("Dau_M1", "std::vector<int>", &Dau_M1_);
    
    std::cout << "Make Branch End" << std::endl;
}

void DeepCMesonAnalyser::resetOnEachJet() {
    jet_label_jpsi_ = 100; 
    jet_label_d0_ = 100; 
    num_jpsiDau_ = 0; 
    num_d0Dau_ = 0; 

    track_pt_.clear();
    track_deta_.clear();
    track_dphi_.clear();
    track_charge_.clear();
    track_pid_.clear();
    
    track_d0_.clear();
    track_dz_.clear();
    track_l_.clear();
    track_xd_.clear();
    track_yd_.clear();
    track_zd_.clear();    
    
    tower_pt_.clear();
    tower_deta_.clear();
    tower_dphi_.clear();
    tower_charge_.clear();
    
    Dau_jpsi_.clear();
    Dau_d0_.clear();
    Dau_M1_.clear();

}


void DeepCMesonAnalyser::analyse(Int_t entry) {
    for (Int_t idx_jet = 0; idx_jet < jets_->GetEntries(); idx_jet++) {
        resetOnEachJet();
        auto jet = dynamic_cast<Jet*>(jets_->At(idx_jet));
       
        Float_t deta, dphi;
        for (Int_t idx_dau = 0; idx_dau < jet -> Constituents.GetEntries(); idx_dau++) {
            TObject* daughter = jet->Constituents.At(idx_dau);
            TLorentzVector jet_p4 = jet->P4();

            if (auto tower = dynamic_cast<Tower*>(daughter)) {
                deta = tower->Eta - jet->Eta;
                dphi = tower->P4().DeltaPhi(jet_p4);
          
                auto gen = dynamic_cast<const GenParticle*>(tower->Particles.At(0));
                if ( gen->M1 != -1 ) {
                    auto mother1 = dynamic_cast<const GenParticle*>(particles_->At(gen->M1));
                    Dau_M1_.push_back(mother1->PID);
                    if (mother1->PID == 443 and gen->M2 == -1) {
                        Dau_jpsi_.push_back(1);
                        Dau_d0_.push_back(0);
                        num_jpsiDau_ = num_jpsiDau_ + 1;
                    }
                    else if (mother1->PID == 421 and gen->M2 == -1) {
                        Dau_jpsi_.push_back(0);
                        Dau_d0_.push_back(1);
                        num_d0Dau_ = num_d0Dau_ + 1;
                    } 
                    else {
                        Dau_jpsi_.push_back(0);
                        Dau_d0_.push_back(0);
                    }
                }
                
                tower_pt_.push_back(tower->ET);
                tower_deta_.push_back(deta);
                tower_dphi_.push_back(dphi);
                tower_charge_.push_back(0);
            
            }
            else if (auto track = dynamic_cast<Track*>(daughter)) {
                deta = track->Eta - jet->Eta;
                dphi = track->P4().DeltaPhi(jet_p4);

                auto gen = dynamic_cast<const GenParticle*>(track->Particle.GetObject());
                if ( gen->M1 != -1 ) {
                    auto mother1 = dynamic_cast<const GenParticle*>(particles_->At(gen->M1));
                    Dau_M1_.push_back(mother1->PID);
                    if (mother1->PID == 443 and gen->M2 == -1) {
                        Dau_jpsi_.push_back(1);
                        Dau_d0_.push_back(0);
                        num_jpsiDau_ = num_jpsiDau_ + 1;
                    }
                    else if (mother1->PID == 421 and gen->M2 == -1) {
                        Dau_jpsi_.push_back(0);
                        Dau_d0_.push_back(1);
                        num_d0Dau_ = num_d0Dau_ + 1;
                    } 
                    else {
                        Dau_jpsi_.push_back(0);
                        Dau_d0_.push_back(0);
                    }
                }
               

                track_pt_.push_back(track->PT);
                track_deta_.push_back(deta);
                track_dphi_.push_back(dphi);
                track_charge_.push_back(track->Charge);
                track_pid_.push_back(track->PID);
                
                track_d0_.push_back(track->D0);
                track_dz_.push_back(track->DZ);
                track_l_.push_back(track->L);
                track_xd_.push_back(track->Xd);
                track_yd_.push_back(track->Yd);
                track_zd_.push_back(track->Zd);
            }
            
        } // daughter 
   
    if ( num_jpsiDau_ == 2) { jet_label_jpsi_ = 1; }  
    else if ( num_jpsiDau_ == 0) { jet_label_jpsi_ = 0; }  
    else { jet_label_jpsi_ = -1; }  
    if ( num_d0Dau_ == 2) { jet_label_d0_ = 1; }  
    else if ( num_d0Dau_ == 0) { jet_label_d0_ = 0; }  
    else { jet_label_d0_ = -1; }  

    out_tree_->Fill();

    } // jet
}


void DeepCMesonAnalyser::Loop() {
  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);
    //if (not selectEvent()) continue;
    analyse(entry);
  }
}
