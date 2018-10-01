#ifndef DELPHYS_EXTERNAL_TREESHUFFLER_H_
#define DELPHYS_EXTERNAL_TREESHUFFLER_H_

#include "TString.h"

#include <vector>
#include <map>

class TTree;
class TFile;

class TreeShuffler {
 public:
  TreeShuffler(std::vector<TString> paths, TString tree_name);
  ~TreeShuffler();

  TTree* CloneTree();
  void CopyAddress(TTree* tree);
  void GetEntry();
  Int_t GetEntries();
  Bool_t IsExhausted();

 private:
  struct TreeData {
    TString path;
    TFile*  file;
    TTree*  tree;
    Int_t   entry;
    Int_t   num_entries;
  };


  TString tree_name_;
  std::map<Int_t, TreeData> candidates_;

  std::vector<UInt_t> alives_;
  Int_t num_alives_;

  Int_t num_total_;
  Int_t num_total_entries_;
};



#endif // DELPHYS_EXTERNAL_TREESHUFFLER_H_
