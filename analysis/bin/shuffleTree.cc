#include "delphys/external/interface/TreeShuffler.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>


int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "shufflerTrees [tree_name] [out_path] [input_paths ...]" << std::endl;
    return 1;
  }

  TString tree_name = argv[1];
  TString out_path = argv[2];

  std::vector<TString> in_paths;
  for (Int_t arg_idx = 3; arg_idx < argc; arg_idx++) {
    in_paths.push_back(TString(argv[arg_idx]));
  }

  auto shuffler = new TreeShuffler(in_paths, tree_name);
  const Int_t kTotalInputEntries = shuffler->GetEntries();

  auto out_file = TFile::Open(out_path, "RECREATE");
  TTree* out_tree = shuffler->CloneTree();
  shuffler->CopyAddress(out_tree);
  out_tree->SetDirectory(out_file);

  Int_t kPrintFreq = kTotalInputEntries / 10;

  for (Int_t entry = 0; entry < kTotalInputEntries; entry++) {
    if (entry % kPrintFreq == 0) {
      TString msg = TString::Format(
          "[%d/%d] %.2f",
          entry,
          kTotalInputEntries,
          static_cast<Float_t>(entry) / kTotalInputEntries * 100);
      std::cout << msg << " %" << std::endl;
    }


    shuffler->GetEntry();
    out_tree->Fill();
  }

  out_tree->Print();
  out_file->Write();
  out_file->Close();

  return 0;
}
