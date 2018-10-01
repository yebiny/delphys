#include "delphys/external/interface/TreeShuffler.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

#include <iostream>
#include <cstdlib> // std::rand, std::abort
#include <memory>


TreeShuffler::TreeShuffler(std::vector<TString> paths,
                           TString tree_name)
    : tree_name_(tree_name) {

  num_total_ = static_cast<Int_t>(paths.size());

  num_total_entries_ = 0;
  for (Int_t idx=0; idx < num_total_; idx++) {
    TString & path = paths.at(idx);

    // FIXME unique_ptr
    TFile* root_file = TFile::Open(path, "READ");
    TTree* tree = dynamic_cast<TTree*>(root_file->Get(tree_name_));

    Int_t num_entries = tree->GetEntries();
    num_total_entries_ += num_entries;

    candidates_[idx] = {path, root_file, tree, 0, num_entries};

    alives_.push_back(idx);

  }
  num_alives_ = static_cast<Int_t>(alives_.size());
}


TreeShuffler::~TreeShuffler() {
  for (auto each : candidates_) {
    each.second.file->Close();
  }
}


void TreeShuffler::GetEntry() {
  // random choice
  Int_t alives_idx = std::rand() % alives_.size();
  Int_t idx = alives_[alives_idx];

  candidates_[idx].tree->GetEntry(candidates_[idx].entry);
  candidates_[idx].entry++;

  if (candidates_[idx].entry == candidates_[idx].num_entries) {
    alives_.erase(alives_.begin() + alives_idx);
    num_alives_ = static_cast<Int_t>(alives_.size());
    std::cout << candidates_[idx].path << "'s tree is exhausted" << std::endl;
    std::cout << TString::Format("Total: %d | Remaining: %d\n", num_total_, num_alives_) << std::endl;
  }
}


TTree* TreeShuffler::CloneTree() {
  return candidates_[0].tree->CloneTree(0);
}

void TreeShuffler::CopyAddress(TTree* tree) {
  for(const auto & each : candidates_)
    tree->CopyAddresses(each.second.tree);
}


Bool_t TreeShuffler::IsExhausted() {
  Bool_t is_exhausted;
  if (num_alives_ > 0) {
    is_exhausted = true;
  } else if (num_alives_ == 0) {
    is_exhausted = false;
  } else {
    std::abort();
  }

  return is_exhausted;
}


Int_t TreeShuffler::GetEntries() {
  return num_total_entries_;
}
