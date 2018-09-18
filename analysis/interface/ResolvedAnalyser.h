#ifndef DELPHYS_ANALYSIS_RESOLVEDALYSER_H_
#define DELPHYS_ANALYSIS_RESOLVEDALYSER_H_

#include "delphys/analysis/interface/BaseAnalyser.h"





class ResolvedAnalyser : private BaseAnalyser {
 public:
  ResolvedAnalyser(const TString & in_path,
                   const TString & out_path,
                   const TString & out_tree_name);
  ~ResolvedAnalyser();
  void Loop();

 private:
  // inherited
  void MakeBranch() override;
  void Reset() override;
  Bool_t SelectEvent() override;
  void AnalyseEvent() override;

  // Bool_t VetoElectron();
  // Bool_t VetoMuon();
  std::vector<Jet> SelectJet();
  // TClonesArray SelectJet();
  
  std::vector<TLorentzVector> jet_p4_;
  std::vector<TLorentzVector> j0_p4_;
};
#endif //  DELPHYS_ANALYSIS_RESOLVEDALYSER_H_
