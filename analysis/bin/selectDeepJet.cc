#include "delphys/analysis/interface/DeepJetSelector.h"
#include "TSystem.h"
#include <cstdlib>

int main(int argc, char* argv[]) {


    TString in_path(argv[1]);
    TString out_path(argv[2]);
    //Int_t out_number = std::atoi(argv[3]);


    DeepJetSelector analyser(in_path, out_path);
    analyser.loop();

    return 0;

}
                
