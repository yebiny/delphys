#ifndef DELPHYS_ANALYSIS_TTALLJETSSELECTOR_H_
#define DELPHYS_ANALYSIS_TTALLJETSSELECTOR_H_

#include "delphys/analysis/interface/BaseAnalyser.h"


class TTAllJetsSelector : private BaseAnalyser {
 public:
  TTAllJetsSelector(const TString & in_path,
                    const TString & out_path);
  ~TTAllJetsSelector();
  void Loop();

 private:
  void MakeBranch();
  void Reset();
  Bool_t SelectEvent();
  void Analyse();

  Bool_t IsAllJetsChannel();

  std::map<Int_t, Int_t> decay_channel_count_;
};

#endif //  DELPHYS_ANALYSIS_TTALLJETSSELECTOR_H_
