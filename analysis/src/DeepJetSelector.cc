#include "delphys/analysis/interface/DeepJetSelector.h"
#include "TLorentzVector.h"

#include <algorithm> // std::fill
#include <cmath> // std::abs, std::fabs, std::sqrt
#include <tuple> // tie, ignore
#include <iterator>
#include <numeric> // std::accumulate, std::inner_product

#include <iostream>

DeepJetSelector::DeepJetSelector(const TString & in_path,
                                 const TString & out_path
                                 ) :
  num_bkg_d0_(0),
  num_sig_d0_(0){

    std::cout << "ctor begin" << std::endl;
    in_file_ = TFile::Open(in_path, "READ");
    in_tree_ = dynamic_cast<TTree*>(in_file_->Get("delphys"));
    
    out_file_ = TFile::Open(out_path, "RECREATE");
    out_tree_ = in_tree_->CloneTree(0);
    out_tree_->SetDirectory(out_file_);

    setBranchAddress();
    std::cout << "ctor end" << std::endl;
}


DeepJetSelector::~DeepJetSelector() {
    std::cout << "dtor begin" << std::endl;
    out_tree_->Print();
    out_file_->Write();
    out_file_->Close();
    std::cout << "dtor end" << std::endl;

}


void DeepJetSelector::setBranchAddress() {
    std::cout << "setbranch begin" << std::endl;
    in_tree_->SetBranchAddress("jet_label_d0", &jet_label_d0_);
    std::cout << "setbranch end" << std::endl;
}


void DeepJetSelector::analyse(Int_t entry) {
    
    if (jet_label_d0_ == 1) { 
        out_tree_->Fill();
        num_sig_d0_++;
            if (num_bkg_d0_%5000 == 0)
            std::cout << "!!SIG!!"<<num_sig_d0_ << std::endl;
    }
    else {
        num_bkg_d0_++;
        if (num_bkg_d0_% 23 == 0){
            out_tree_->Fill();
            //std::cout << "!!BKG!!"<<num_bkg_d0_ << std::endl;
        }
       
    }
}

void DeepJetSelector::Loop() {

  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);
    analyse(entry);
  }
}
