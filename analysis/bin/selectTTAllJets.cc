#include "delphys/analysis/interface/TTAllJetsSelector.h"

int main(int argc, char* argv[]) {
  TString in_path(argv[1]);
  TString out_path(argv[2]);

  TTAllJetsSelector analyser(in_path, out_path);
  analyser.Loop();

  return 0;
}
