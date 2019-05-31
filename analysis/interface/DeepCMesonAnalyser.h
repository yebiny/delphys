#ifndef DELPHYS_ANALYSIS_DEEPCMESONANALYSER_H_
#define DELPHYS_ANALYSIS_DEEPCMESONANALYSER_H_

#include "delphys/analysis/interface/BaseAnalyser.h"

class DeepCMesonAnalyser : private BaseAnalyser  {
    public:
        DeepCMesonAnalyser( const TString & in_path,
                            const TString & out_path,
                            const TString & out_tree_name);
        ~DeepCMesonAnalyser();
        void loop();

    private:
        void makeBranch();
        void analyse(Int_t entry);
        void resetOnEachJet();

        std::vector<TLorentzVector> track_P4_;
        TLorentzVector jet_track_comb_;
        std::vector<Int_t>          Idx;
        
        std::vector<Float_t>        track_deta_;
        std::vector<Float_t>        track_dphi_;
        
        std::vector<Float_t>        track_pt_;
        std::vector<Float_t>        track_eta_;
        std::vector<Float_t>        track_phi_;
        std::vector<Float_t>        track_mass_;
        
        std::vector<Float_t>        track_d0_;
        std::vector<Float_t>        track_dz_;
        std::vector<Float_t>        track_xd_;
        std::vector<Float_t>        track_yd_;
        std::vector<Float_t>        track_zd_;
        std::vector<Float_t>        track_errd0_;
       
        std::vector<Float_t>        track_l_;
        std::vector<Int_t>          track_charge_;
        std::vector<Int_t>          track_pId_;

        std::vector<Int_t>          track_costompId_;
        std::vector<Int_t>          mother_pId_;
        
        std::vector<Int_t>          pticle_is_d0dau_; 
        std::vector<Int_t>          pticle_is_pion_;
        std::vector<Int_t>          pticle_is_kaon_;
        
        std::vector<Int_t>          pticle_label_;
        Int_t                       jet_label_;
        Int_t                       jet_num_d0dau_;
        Int_t                       jet_num_track_;
        
        Int_t                       num_pion_cand_;
        Int_t                       num_kaon_cand_;
        std::vector<Int_t>          charge_pion_cand_;
        std::vector<Int_t>          charge_kaon_cand_;
       
        static const int pion_pId_ = 211, kaon_pId_ = 321;
        static const int d0_pId_ = 421;   
};

#endif
