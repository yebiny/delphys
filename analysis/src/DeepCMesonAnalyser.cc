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

    out_tree_->Branch("d0_gen_p4", "TLorentzVector", &d0_gen_p4_);
    out_tree_->Branch("d0_rec_p4", "TLorentzVector", &d0_rec_p4_);
    
    out_tree_->Branch("track_deta", "std::vector<float>", &track_deta_);
    out_tree_->Branch("track_dphi", "std::vector<float>", &track_dphi_);
   
    out_tree_->Branch("track_pt", "std::vector<float>", &track_pt_);
    out_tree_->Branch("track_eta", "std::vector<float>", &track_eta_);
    out_tree_->Branch("track_phi", "std::vector<float>", &track_phi_);
    
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
    
    out_tree_->Branch("dau_label", "std::vector<int>", &dau_label_);

    out_tree_->Branch("gen_charge", "std::vector<int>", &gen_charge_);
    
    out_tree_->Branch("jet_label", &jet_label_, "jet_label/I");
    out_tree_->Branch("jet_count_d0dau", &jet_count_d0dau_, "jet_count_d0dau/I");
    out_tree_->Branch("jet_count_track", &jet_count_track_, "jet_count_track/I");
    out_tree_->Branch("jet_count_pion", &jet_count_pion_, "jet_count_pion/I");
    out_tree_->Branch("jet_count_kaon", &jet_count_kaon_, "jet_count_kaon/I");

    std::cout << "Make Branch End" << std::endl;
}

void DeepCMesonAnalyser::resetOnEachJet() {
    pion_gen_p4_.SetPtEtaPhiM(0,0,0,0);
    kaon_gen_p4_.SetPtEtaPhiM(0,0,0,0);
    pion_rec_p4_.SetPtEtaPhiM(0,0,0,0);
    kaon_rec_p4_.SetPtEtaPhiM(0,0,0,0);
    d0_gen_p4_.SetPtEtaPhiM(0,0,0,0);
    d0_rec_p4_.SetPtEtaPhiM(0,0,0,0);
    
    track_deta_.clear();
    track_dphi_.clear();
    
    track_pt_.clear();
    track_eta_.clear();
    track_phi_.clear();
    
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
    mother_num_.clear();

    dau_label_.clear();    
    
    gen_charge_.clear();
    gen_pt_.clear();
    gen_eta_.clear();
    gen_phi_.clear();
    gen_mass_.clear();
    
    jet_label_ = -1; 
    jet_count_d0dau_ = 0;
    jet_count_track_ = 0;
    jet_count_pion_ = 0;
    jet_count_kaon_ = 0;

    kaon_Idx.clear();
    pion_Idx.clear();
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

            // Check if track particle
            if (auto track = dynamic_cast<Track*>(daughter)) {
                
                // Check has mother particle
                auto gen = dynamic_cast<const GenParticle*>(track->Particle.GetObject());
                if (gen->M1 == -1) continue;
                
                // Now, save all track informations
                jet_count_track_ = jet_count_track_ + 1;
                
                deta = track->Eta - jet->Eta;
                dphi = track->P4().DeltaPhi(jet_p4);
                
                track_deta_.push_back(deta);
                track_dphi_.push_back(dphi);
                
                track_pt_.push_back(track->PT);
                track_eta_.push_back(track->Eta);
                track_phi_.push_back(track->Phi);
                
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
               
                // check mother particle info
                auto mother1 = dynamic_cast<const GenParticle*>(particles_->At(gen->M1));
                mother_pId_.push_back(mother1->PID);
                mother_num_.push_back(gen->M1);
                
                gen_charge_.push_back(gen->Charge);
                gen_pt_.push_back(gen->PT);
                gen_eta_.push_back(gen->Eta);
                gen_phi_.push_back(gen->Phi);
                gen_mass_.push_back(gen->Mass);
                
                // Check if particle is from D0
                if (abs(mother1->PID) == d0_pId_ and gen->M2 == -1) {
                    jet_count_d0dau_ = jet_count_d0dau_ + 1; 
                    
                    // Check if particle pId is kaon or pion 
                    if (abs(gen->PID) == kaon_pId_) {
                        
                        dau_label_.push_back(1);
                        jet_count_kaon_ = jet_count_kaon_+1;
                    
                    } else if (abs(gen->PID) == pion_pId_) {
                       
                        dau_label_.push_back(2);
                        jet_count_pion_ = jet_count_pion_+1;

                    } else{
                        dau_label_.push_back(-1);
                    }
                } else{
                    dau_label_.push_back(0);
                }//d0 daughter
            } else continue;//track
        }//constituents 
    
        Int_t size = dau_label_.size(); 
        // Get Particle Candidate
        for ( Int_t i=0; i < size; i++ ) {
            if (dau_label_[i] == 1) {
                kaon_Idx.push_back(i);
            }
            else if (dau_label_[i] == 2) {
                pion_Idx.push_back(i);
            }
        }
        Int_t kaon_size = kaon_Idx.size(); 
        Int_t pion_size = pion_Idx.size(); 
        
        // Jet labelling    
        jet_label_ = 0;
        if ( (jet_count_d0dau_ >=2) and (jet_count_pion_>=1) and (jet_count_kaon_>=1) ){
            jet_label_ = 1;
            for (Int_t i=0; i<kaon_size; i++){
                Int_t knum = kaon_Idx[i]; 

                for (Int_t j=0; j<pion_size; j++){
                    Int_t pnum = pion_Idx[j]; 
                    
                    if (mother_num_[knum] != mother_num_[pnum]) continue;
                    if (mother_pId_[knum] != mother_pId_[pnum]) continue;
                    if (gen_charge_[knum]*gen_charge_[pnum] >0) continue;
                        jet_label_=2;
                        
                        Int_t mnum = mother_num_[knum];
                        Int_t n_mnum = 0;
                        
                        for (Int_t k=0; k<size; k++){
                            if ((mnum == mother_num_[k]) and (abs(mother_pId_[k]) == 421)){
                                n_mnum = n_mnum+ 1;
                            }
                        }

                        if ( n_mnum == 2 ){
                            jet_label_=3;

                            kaon_rec_p4_.SetPtEtaPhiM(track_pt_[knum],track_eta_[knum],track_phi_[knum],gen_mass_[knum]);
                            pion_rec_p4_.SetPtEtaPhiM(track_pt_[pnum],track_eta_[pnum],track_phi_[pnum],gen_mass_[pnum]);
                            d0_rec_p4_ = kaon_rec_p4_+pion_rec_p4_;
                            kaon_gen_p4_.SetPtEtaPhiM(gen_pt_[knum],gen_eta_[knum],gen_phi_[knum],gen_mass_[knum]);
                            pion_gen_p4_.SetPtEtaPhiM(gen_pt_[pnum],gen_eta_[pnum],gen_phi_[pnum],gen_mass_[pnum]);
                            d0_gen_p4_ = kaon_gen_p4_+pion_gen_p4_;

                            if ( abs(d0_gen_p4_.M() - d0_m_) < 0.05 ) {
                                jet_label_=4;
                                
                                //Particle relabelling for reconstruction
                                if (pnum <5 and knum <5 ){
                                    jet_label_ =5;
                                    replace(dau_label_.begin(), dau_label_.end(), 1, 0);
                                    replace(dau_label_.begin(), dau_label_.end(), 2, 0);
                                    dau_label_[pnum]=1;
                                    dau_label_[knum]=1;
                                }

                            }else {
                                d0_gen_p4_.SetPtEtaPhiM(0,0,0,0);
                                d0_rec_p4_.SetPtEtaPhiM(0,0,0,0);
                            } 
                        }//label 3
                }//label 2: pion
            }//label 2: kaon  
        }//label 1
        
        //if (jet_label_ != 4) {
        //    replace(dau_label_.begin(), dau_label_.end(), 1, 0);
        //    replace(dau_label_.begin(), dau_label_.end(), 2, 0);
        //}
        
        // Fill only jet which has least 2 track particles
        if (jet_count_track_ < 2) continue;
        out_tree_->Fill();

        //CHECK     
        if (jet_label_ > 3){
            std::cout << "----------------------" << std::endl;
            std::cout << "Jet tag: "<< jet_label_ << std::endl;
            std::cout << "----------------------" << std::endl;
            for ( Int_t i=0; i < size; i++ ) {
                std::cout << dau_label_[i]<<"|"<< mother_num_[i]<<"|"<< mother_pId_[i] <<"|"<< gen_charge_[i]<<std::endl;
            }
        }    
    }//jet
}

void DeepCMesonAnalyser::loop() {
  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);
    analyse(entry);
  }
} 
