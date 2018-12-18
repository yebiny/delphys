import ROOT, math
from ROOT import gStyle, gROOT

color_list = {"tt":ROOT.kGreen+2, "hh_SM":ROOT.kRed, "hh_B6":ROOT.kBlue+1, "hh_B11":ROOT.kViolet}

f = ROOT.TFile("plots.root")
cvs = ROOT.TCanvas()
cvs.SetGrid()
cvs.SetLeftMargin(0.15)
cvs.SetBottomMargin(0.12)
gStyle.SetPadTickX(True)
gStyle.SetPadTickY(True)
gStyle.SetLineWidth(2)

h_SM = f.Get("hh afc basic MT2 332 bbll")
h_B6 = f.Get("hh_B6 afc basic MT2 332 bbll")
h_B11 = f.Get("hh_B11 afc basic MT2 332 bbll")
h_tt = f.Get("tt afc basic MT2 332 blbl")
hlist = {"hh_SM":h_SM, "tt":h_tt, "hh_B6":h_B6, "hh_B11":h_B11}
legend = ROOT.TLegend(0.70,0.55,0.90,0.90)
for sample in ["hh_SM","tt","hh_B6","hh_B11"]:
    h = hlist[sample]
    h.SetLineColor(color_list[sample])
    h.SetTitle("")
    h.SetStats(False)
    h.SetLabelSize(0.05,"X")
    h.SetLabelSize(0.04,"Y")
    h.SetLineWidth(3)
    h.GetXaxis().SetTitle("MT2 [GeV]")
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetAxisColor(ROOT.kBlack)
    cvs.RedrawAxis()
    h.GetYaxis().SetAxisColor(ROOT.kBlack)
    h.GetYaxis().SetTitle("Normalized to unity")
    h.GetYaxis().SetTitleSize(0.05)
    integ = h.Integral()
    if integ != 0:
        h.Scale(1/integ)
    if sample=="hh":
        h.Draw("hist")
    else:
        h.Draw("hist same")
for sample in ["tt","hh_SM","hh_B6","hh_B11"]:
    h = hlist[sample]
    legend.AddEntry(h,sample)
legend.SetTextSize(0.06)
legend.Draw("same")
cvs.SetTitle("MT2")
cvs.SaveAs("MT2.png")

f.Close()
