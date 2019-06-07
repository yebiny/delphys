#include "delphys/external/interface/TreeShuffler.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

#include <iostream>
#include <cstdlib> // std::rand, std::abort
#include <memory>
#include <algorithm> // find


TreeShuffler::TreeShuffler(std::vector<TString> paths,
                           TString tree_name)
    : tree_name_(tree_name) {

  num_total_ = paths.size();

  num_total_entries_ = 0;
  for (UInt_t idx=0; idx < num_total_; idx++) {
    TString & path = paths.at(idx);

    TFile* root_file = TFile::Open(path, "READ");
    TTree* tree = dynamic_cast<TTree*>(root_file->Get(tree_name_));

    Int_t num_entries = tree->GetEntries();
    num_total_entries_ += num_entries;

    candidates_[idx] = {path, root_file, tree, 0, num_entries};

    weights_.push_back(num_entries);
  }

  engine_ = std::mt19937(device_());
  dist_ = std::discrete_distribution<UInt_t>(weights_.begin(), weights_.end());
  progress_ = 0;
}


TreeShuffler::~TreeShuffler() {
  for (auto each : candidates_) {
    each.second.file->Close();
  }
}


UInt_t TreeShuffler::RandomChoice() {
  UInt_t index = dist_(engine_);

  // check if this index correspond
  auto iter = std::find(exhausted_.begin(), exhausted_.end(), index);

  if (iter != exhausted_.end()) {
    return RandomChoice();
  }
  return index;
}


void TreeShuffler::GetEntry() {
  progress_++;

  UInt_t idx = RandomChoice();

  candidates_[idx].tree->GetEntry(candidates_[idx].entry);
  candidates_[idx].entry++;

  if (candidates_[idx].entry == candidates_[idx].num_entries) {
    exhausted_.push_back(idx);

    std::cout << candidates_[idx].path << "'s tree is exhausted" << std::endl;

    std::cout << progress_ << " / " << num_total_entries_ << std::endl;

    Int_t num_alives = num_total_ - exhausted_.size();
    std::cout << TString::Format("Total: %d | Remaining: %d\n",
                                 num_total_, num_alives)
              << std::endl;
  }
}


TTree* TreeShuffler::CloneTree() {
  return candidates_[0].tree->CloneTree(0);
}

void TreeShuffler::CopyAddress(TTree* tree) {
  for (const auto & each : candidates_)
    tree->CopyAddresses(each.second.tree);
}


Bool_t TreeShuffler::IsExhausted() {
  return (exhausted_.size() == num_total_);
}

Int_t TreeShuffler::GetEntries() {
  return num_total_entries_;
}
