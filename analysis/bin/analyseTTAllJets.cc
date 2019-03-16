#include "delphys/analysis/interface/TTAllJetsAnalyser.h"

int main(int argc, char* argv[]) {
  TString in_path(argv[1]);
  TString out_path(argv[2]);

  TTAllJetsAnalyser analyser(in_path, out_path, "delphys");
  analyser.Loop();

  return 0;
}
