#ifndef DELPHYS_ANALYSIS_TTALLJETSFILTER_H_
#define DELPHYS_ANALYSIS_TTALLJETSFILTER_H_

#include "delphys/analysis/interface/BaseAnalyser.h"


class TTAllJetsFilter : private BaseAnalyser {
 public:
  TTAllJetsFilter(const TString & in_path,
                  const TString & out_path);
  ~TTAllJetsFilter();
  void Loop();

 private:
  void MakeBranch();
  void Reset();
  Bool_t SelectEvent();
  void AnalyseEvent();

  Bool_t IsAllJetsChannel();

  std::map<Int_t, Int_t> decay_channel_count_;
};

#endif //  DELPHYS_ANALYSIS_TTALLJETSFILTER_H_
