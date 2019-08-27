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

        TLorentzVector              pion_gen_p4_;
        TLorentzVector              kaon_gen_p4_;
        TLorentzVector              pion_rec_p4_;
        TLorentzVector              kaon_rec_p4_;
        
        TLorentzVector              d0_gen_p4_;
        TLorentzVector              d0_rec_p4_;
        TLorentzVector              d0_gen_;
        TLorentzVector              d0_rec_;
        
        std::vector<Float_t>        track_deta_;
        std::vector<Float_t>        track_dphi_;
        
        std::vector<Float_t>        track_pt_;
        std::vector<Float_t>        track_eta_;
        std::vector<Float_t>        track_phi_;
        
        std::vector<Float_t>        track_d0_;
        std::vector<Float_t>        track_dz_;
        std::vector<Float_t>        track_xd_;
        std::vector<Float_t>        track_yd_;
        std::vector<Float_t>        track_zd_;
        std::vector<Float_t>        track_errd0_;
        std::vector<Float_t>        track_l_;
        
        std::vector<Int_t>          track_pId_;
        std::vector<Int_t>          track_costompId_;
        std::vector<Int_t>          track_charge_;
        
        std::vector<Int_t>          mother_pId_;
        std::vector<Int_t>          mother_num_;
        
        std::vector<Int_t>          dau_label_;
        
        std::vector<Int_t>          gen_charge_;
        std::vector<Float_t>        gen_pt_;
        std::vector<Float_t>        gen_eta_;
        std::vector<Float_t>        gen_phi_;
        std::vector<Float_t>        gen_mass_;
        
        Int_t                       jet_label_;
        Int_t                       jet_count_d0dau_;
        Int_t                       jet_count_track_;
        Int_t                       jet_count_pion_;
        Int_t                       jet_count_kaon_;
       
        std::vector<Int_t>          kaon_Idx;
        std::vector<Int_t>          pion_Idx;
       
        static const int pion_pId_ = 211, kaon_pId_ = 321, d0_pId_ = 421;
        static constexpr float pion_m_ = 0.1396, kaon_m_ = 0.4937, d0_m_ = 1.865;   
};

#endif
