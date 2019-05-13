#include "delphys/analysis/interface/DeepCMesonAnalyser.h"

#include "TSystem.h"


int main(int argc, char* argv[]) {
  TString in_path(argv[1]);
  TString out_path(argv[2]);

  DeepCMesonAnalyser analyser(in_path, out_path, "delphys");
  analyser.loop();
 
  return 0;
}
