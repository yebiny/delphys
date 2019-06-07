#ifndef DELPHYS_ANALYSIS_TTALLJETSSELECTOR_H_
#define DELPHYS_ANALYSIS_TTALLJETSSELECTOR_H_

#include "delphys/analysis/interface/BaseAnalyser.h"

namespace pdgid {
  static const Int_t kBottom = 5;
  static const Int_t kTop = 6;
  static const Int_t kElectron = 11;
  static const Int_t kMuon = 13;
  static const Int_t kTau = 15;
  static const Int_t kPhoton = 22;
  static const Int_t kWBoson = 24;
} // pdgid

static const std::set<Int_t> kSkipPID = {
    1, 2, 3, 4, 5, // quarks
    12, 14, 16 // neutrinos
};


namespace p8status {
  static const Int_t kTop = 62;
} // p8status


class TTAllJetsSelector : private BaseAnalyser {
 public:
  TTAllJetsSelector(const TString & in_path,
                    const TString & out_path);
  ~TTAllJetsSelector();
  void Loop();

 private:
  Bool_t selectEvent();
  void analyse();

  Bool_t isAllJetsChannel();

  std::map<Int_t, Int_t> decay_channel_count_;

};

#endif //  DELPHYS_ANALYSIS_TTALLJETSSELECTOR_H_
