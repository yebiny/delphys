#include "delphys/analysis/interface/TTAllJetsAnalyser.h"

#include "TSystem.h"

int main(int argc, char* argv[]) {
  TString in_path(argv[1]);
  TString out_path(argv[2]);

  TString in_name = gSystem->BaseName(in_path);

  Int_t label = -1;
  if (in_name.Contains("TT")) {
    label = 1;
  } else if (in_name.Contains("QCD")) {
    label = 0;
  } else {
    return 1;
  } 

  TTAllJetsAnalyser analyser(in_path, out_path, label);
  analyser.Loop();

  return 0;
}
