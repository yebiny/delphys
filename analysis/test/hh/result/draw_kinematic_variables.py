import ROOT, math
from ROOT import gStyle, gROOT

hlist = { "bfc": {"tt":dict(), "hh":dict(), "hh_B6":dict(), "hh_B11":dict()}, "afc" : {"tt":dict(), "hh":dict(), "hh_B6":dict(), "hh_B11":dict()} }
sample_list = ["tt", "hh", "hh_B6", "hh_B11" ]
order_list = { 'll_pt': { "bfc":[0,3,2,1], "afc":[0,2,3,1] },
               'bb_pt': { "bfc":[0,2,3,1], "afc":[0,2,3,1] }, 
               'missing_et': { "bfc":[0,1,2,3], "afc":[0,1,2,3]}, 
               'missing_et_phi': { "bfc":[1,2,3,0], "afc":[0,2,3,1] }, 
               'bbll_mass': { "bfc":[1,2,3,0], "afc":[0,2,3,1] }, 
               "bb_mass": { "bfc":[1,2,3,0], "afc":[3,2,1,0] }, 
               "ll_mass": { "bfc":[3,2,1,0], "afc":[2,3,1,0] },
               "basic_MT2_332_bbll": { "bfc":[3,2,1,0], "afc":[2,3,1,0] },
               "basic_MT2_332_b": { "bfc":[3,2,1,0], "afc":[0,1,2,3] },
               "basic_MT2_332_l": { "bfc":[3,2,1,0], "afc":[2,3,1,0] },
               "basic_MT2_332_blbl": { "bfc":[3,2,1,0], "afc":[0,1,2,3] },
               "ch_bisect_MT2_332" : { "bfc":[3,2,1,0], "afc":[2,3,1,0] },
               "ch_bisect_MT2_332_b" : { "bfc":[3,2,1,0], "afc":[0,1,2,3] },
               "ch_bisect_MT2_332_l" : { "bfc":[3,2,1,0], "afc":[2,3,1,0] },
               "bb_deltaR": { "bfc":[3,2,1,0], "afc":[2,3,1,0] },
               "ll_deltaR": { "bfc":[3,2,1,0], "afc":[2,3,1,0] },
               "bbll_deltaR": { "bfc":[1,2,3,0], "afc":[2,3,1,0] },
               "bl_deltaR": { "bfc":[3,2,1,0], "afc":[3,2,1,0] },
               "bl_min_deltaR": { "bfc":[0,1,2,3], "afc":[0,2,3,1] },
               "mT": { "bfc":[3,2,1,0], "afc":[2,3,1,0] }
              }
color_list = {"tt":ROOT.kGreen+2, "hh":ROOT.kRed, "hh_B6":ROOT.kBlue+1, "hh_B11":ROOT.kViolet}

f = ROOT.TFile("plots.root")
cvs = ROOT.TCanvas()
cvs.SetGrid()
cvs.SetLeftMargin(0.15)
cvs.SetBottomMargin(0.12)

for key in ["bfc","afc"]:
    for key2 in ["hh","hh_B6","hh_B11","tt"]:
        for key3 in order_list.keys():
            hlist[key][key2][key3] = f.Get(key2+" "+key+" "+key3.replace("_"," "))

for key in ["bfc","afc"]:
    for key3 in order_list.keys():
        legend = ROOT.TLegend(0.70,0.55,0.90,0.90)
        for i in [0,1,2,3]:
            #legend = ROOT.TLegend(0.78,0.5,0.98,0.75)
            h = hlist[key][sample_list[i]][key3]
            h.SetLineColor(color_list[sample_list[i]])
            legend.SetTextSize(0.06)
            legend.AddEntry(h,sample_list[i])
            h.SetTitle("")
            h.SetStats(False)
            h.SetLabelSize(0.05,"X")
            h.SetLabelSize(0.05,"Y")
            h.SetLineWidth(3)
            h.SetTickLength(0.03)
            if 'deltaR' in key3:
                h.GetXaxis().SetTitle(key3.replace("_"," "))
            else:
                h.GetXaxis().SetTitle(key3.replace("_"," ")+" [GeV]")
            h.GetXaxis().SetTickSize(0.03)
            h.GetXaxis().SetTitleSize(0.06)
            h.GetXaxis().SetAxisColor(ROOT.kBlack)
            cvs.RedrawAxis()
            h.GetYaxis().SetAxisColor(ROOT.kBlack)
            h.GetYaxis().SetTitle("Normalized to unity")
            h.GetYaxis().SetTitleSize(0.05)
        for ik, ikey2 in enumerate(order_list[key3][key]):
            key2 = sample_list[ikey2]
            h = hlist[key][key2][key3]
            integ = h.Integral()
            if integ != 0:
                h.Scale(1/integ)
            if ik==0:
                h.Draw("hist")
            else:
                h.Draw("hist same")
        legend.Draw("same")
        cvs.SetTitle(key3.replace("_"," "))
        cvs.SaveAs(key3+"_"+key+".png")

f.Close()
