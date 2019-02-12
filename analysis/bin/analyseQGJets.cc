#include "delphys/analysis/interface/QGJetsAnalyser.h"

#include "TSystem.h"


int main(int argc, char* argv[]) {
  TString in_path(argv[1]);
  TString out_path(argv[2]);

  TString in_name = gSystem->BaseName(in_path);

  // quark-enriched sample has even label
  Int_t label;
  Bool_t is_dijet;
  if (in_name.Contains("qq")) {
    label = 0;
    is_dijet = true;
  } else if (in_name.Contains("gg")) {
    label = 1;
    is_dijet = true;
  } else if (in_name.Contains("zq")) {
    label = 2;
    is_dijet = false;
  } else if (in_name.Contains("zg")) {
    label = 3;
    is_dijet = false;
  } else if (in_name.Contains("zj")) {
    label = 4;
    is_dijet = false;
  } else if (in_name.Contains("jj")) {
    label = 5;
    is_dijet = true;
  } else {
    std::cout << in_path << std::endl;
    return 1;
  }


  QGJetsAnalyser analyser(in_path, out_path, "jetAnalyser", is_dijet, label);
  analyser.Loop();
 
  return 0;
}
