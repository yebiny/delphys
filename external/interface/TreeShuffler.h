#ifndef DELPHYS_EXTERNAL_TREESHUFFLER_H_
#define DELPHYS_EXTERNAL_TREESHUFFLER_H_

#include <vector>
#include <map>
#include <random>

#include "TString.h"

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
  UInt_t RandomChoice();

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


  UInt_t num_total_;
  Int_t num_total_entries_;

  std::random_device device_;
  std::mt19937 engine_;
  std::discrete_distribution<UInt_t> dist_;
  std::vector<Int_t> weights_;
  std::vector<UInt_t> exhausted_;
  UInt_t progress_;
};



#endif // DELPHYS_EXTERNAL_TREESHUFFLER_H_
