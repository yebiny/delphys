import ROOT, math
hlist = { "bfc": {"tt":dict(), "hh":dict(), "hh_B6":dict(), "hh_B11":dict()}, "afc" : {"tt":dict(), "hh":dict(), "hh_B6":dict(), "hh_B11":dict()} }
limit_list = { 'll_pt': { "bfc":[60,0,400], "afc":[60,0,400] },
               'bb_pt': { "bfc":[60,0,500], "afc":[60,0,500] }, 
               'missing_et': { "bfc":[60,0,400], "afc":[60,0,400]}, 
               'missing_et_phi': { "bfc":[60,-4,4], "afc":[60,-4,4] }, 
               'bbll_mass': { "bfc":[100,0,700], "afc":[100,0,700] }, 
               "bb_mass": { "bfc":[100,0,400], "afc":[100,50,400] }, 
               "ll_mass": { "bfc":[100,0,400], "afc":[100,0,150] },
               "basic_MT2_332_bbll" : { "bfc":[80,0,500], "afc":[80,0,400] },
               "basic_MT2_332_blbl" : { "bfc":[80,0,500], "afc":[80,0,400] },
               "basic_MT2_332_b" : { "bfc":[80,0,500], "afc":[80,0,400] },
               "basic_MT2_332_l" : { "bfc":[80,0,500], "afc":[80,0,200] },
               "ch_bisect_MT2_332" : { "bfc":[80,0,500], "afc":[80,0,400] },
               "ch_bisect_MT2_332_b" : { "bfc":[80,0,500], "afc":[80,0,400] },
               "ch_bisect_MT2_332_l" : { "bfc":[80,0,500], "afc":[80,0,400] },
               "bb_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "ll_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "bbll_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "bl_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "bl_min_deltaR" : { "bfc":[50,0,5], "afc":[50,0,5] },
               "mT" : { "bfc":[75,0,300], "afc":[75,0,300] },
              }
def getHistogram(title,binning,lower_limit,upper_limit):
    return ROOT.TH1D(title,title,binning,lower_limit,upper_limit)

for key in ["bfc","afc"]:
    for key2 in ["hh","tt","hh_B6","hh_B11"]:
        for key3 in limit_list.keys():
            hlist[key][key2][key3] = getHistogram(key2+" "+key+" "+key3.replace("_"," "),limit_list[key3][key][0],limit_list[key3][key][1],limit_list[key3][key][2])

tt_c = ROOT.TChain("events")
tt_c.Add("TT_200.root")
hh_c = ROOT.TChain("events")
hh_c.Add("HH_SM.root")
hh_B6_c = ROOT.TChain("events")
hh_B6_c.Add("HH_B6.root")
hh_B11_c = ROOT.TChain("events")
hh_B11_c.Add("HH_B11.root")
sample_list = {"tt":tt_c,"hh":hh_c,"hh_B6":hh_B6_c,"hh_B11":hh_B11_c}

f = ROOT.TFile("plots.root","RECREATE")

for sample in ["tt", "hh", "hh_B6", "hh_B11"]:
    c = sample_list[sample]
    for ie, e in enumerate(c):
        cut = ""
        if e.step<4:
            cut = "bfc"
        else: cut = "afc"
        hlist[cut][sample]['ll_pt'].Fill(e.ll.Pt())
        hlist[cut][sample]['bb_pt'].Fill(e.bb.Pt())
        hlist[cut][sample]["missing_et_phi"].Fill(e.MET.Phi())
        hlist[cut][sample]["missing_et"].Fill(e.MET.E())
        hlist[cut][sample]["bbll_mass"].Fill(e.bbll.M())
        hlist[cut][sample]["bb_mass"].Fill(e.bb.M())
        hlist[cut][sample]["ll_mass"].Fill(e.ll.M())
        hlist[cut][sample]["basic_MT2_332_bbll"].Fill(e.basic_MT2_332_bbll)
        hlist[cut][sample]["basic_MT2_332_blbl"].Fill(e.basic_MT2_332_blbl)
        hlist[cut][sample]["basic_MT2_332_b"].Fill(e.basic_MT2_332_b)
        hlist[cut][sample]["basic_MT2_332_l"].Fill(e.basic_MT2_332_l)
        hlist[cut][sample]["ch_bisect_MT2_332"].Fill(e.ch_bisect_MT2_332)
        hlist[cut][sample]["ch_bisect_MT2_332_b"].Fill(e.ch_bisect_MT2_332_b)
        hlist[cut][sample]["ch_bisect_MT2_332_l"].Fill(e.ch_bisect_MT2_332_l)
        hlist[cut][sample]["bb_deltaR"].Fill(e.bb_deltaR)
        hlist[cut][sample]["ll_deltaR"].Fill(e.ll_deltaR)
        hlist[cut][sample]["bbll_deltaR"].Fill(e.bbll_deltaR)
        hlist[cut][sample]["bl_min_deltaR"].Fill(e.bl_min_deltaR)
        hlist[cut][sample]["mT"].Fill(e.mT)
        for bld in e.bl_deltaR:
            hlist[cut][sample]["bl_deltaR"].Fill(bld)

for key in ["bfc","afc"]:
    for key2 in ["tt", "hh", "hh_B6", "hh_B11"]:
        for key3 in hlist[key][key2].keys():
            #integ = hlist[key][key2][key3].Integral()
            #if integ == 0:
            #    print "%s %s %s has 0 integral" % (key, key2, key3)
            #else:
            #    hlist[key][key2][key3].Scale(1/integ)
            hlist[key][key2][key3].Write()

f.Close()
