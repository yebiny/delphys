#ifndef DELPHYS_ANALYSIS_DEEPJETSELECTOR_H_
#define DELPHYS_ANALYSIS_DEEPJETSELECTOR_H_

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

class DeepJetSelector{

    public:
        DeepJetSelector(
                         const TString & in_path,   
                         const TString & out_path   
                       );
        ~DeepJetSelector();
        void loop();

    private:
        void analyse(Int_t entry);
        void setBranchAddress();
        TFile* in_file_;
        TFile* out_file_;
        TTree* in_tree_;
        TTree* out_tree_;
            
        // branches
        Int_t jet_label_d0_;

        // 
        Int_t num_bkg_d0_;
        Int_t num_sig_d0_;
};

#endif

    

