#ifndef DELPHYS_ANALYSIS_DEEPCMESONANALYSER_H_
#define DELPHYS_ANALYSIS_DEEPCMESONANALYSER_H_

#include "delphys/analysis/interface/BaseAnalyser.h"


class DeepCMesonAnalyser : private BaseAnalyser  {
    public:
        DeepCMesonAnalyser( const TString & in_path,
                            const TString & out_path,
                            const TString & out_tree_name);
        ~DeepCMesonAnalyser();
        void Loop();


    private:
        void makeBranch();
        //Bool_t selectEvent();
        void analyse(Int_t entry);
        void resetOnEachJet();
        //void resetOnEachEvent();

        std::vector<Float_t>        track_pt_;
        std::vector<Float_t>        track_deta_;
        std::vector<Float_t>        track_dphi_;
        std::vector<Int_t>          track_charge_;
        std::vector<Int_t>          track_pid_;
       
        std::vector<Float_t>        track_d0_;
        std::vector<Float_t>        track_dz_;
        std::vector<Float_t>        track_l_;
        std::vector<Float_t>        track_xd_;
        std::vector<Float_t>        track_yd_;
        std::vector<Float_t>        track_zd_;



        std::vector<Float_t>        tower_pt_;
        std::vector<Float_t>        tower_deta_;
        std::vector<Float_t>        tower_dphi_;
        std::vector<Int_t>          tower_charge_;
        
        std::vector<Int_t>          Dau_jpsi_;
        std::vector<Int_t>          Dau_d0_;
        std::vector<Int_t>          Dau_M1_;
        
        
        Int_t jet_label_d0_;
        Int_t jet_label_jpsi_;
        Int_t num_d0Dau_;
        Int_t num_jpsiDau_;
        

};

#endif
