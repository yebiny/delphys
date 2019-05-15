from __future__ import division
from collections import OrderedDict
import numpy as np

from ROOT import TH1F
import os
import matplotlib.pyplot as plt
import numpy
import ROOT, sys

folder = './'+sys.argv[1:][0]

info = np.load(folder+"/model_info.npz")
train_sig_response = info['train_sig_response']
train_bkg_response = info['train_bkg_response']
test_sig_response = info['test_sig_response']
test_bkg_response = info['test_bkg_response']

canvas = ROOT.TCanvas("c", "c", 1200, 800)
canvas.cd()
print(  "make canvas")

h0 = TH1F("untitled", "Model resoponse" , 50, 0, 1.)
h0.SetXTitle("Model response")

h_test_bkg = TH1F("test_bkg", "test_bkg" , 50, 0, 1.)
h_train_sig = TH1F("train_sig", "train_sig" , 50, 0, 1.)
h_train_bkg = TH1F("train_bkg", "train_bkg" , 50, 0, 1.)
hists = [
    h0,
    h_test_bkg,
    h_train_sig,
    h_train_bkg
]
print(  "save hists")
 
for each in train_sig_response:
    h_train_sig.SetFillColor(3)
    h_train_sig.SetFillStyle(3004)
    h_train_sig.Fill(each)
print( "fill train sig"   ) 
for each in train_bkg_response:
    h_train_bkg.SetFillColor(5)
    h_train_bkg.SetFillStyle(3004)
    h_train_bkg.Fill(each)
print( "fill train bkg")
for each in test_sig_response: 
    h0.SetLineColor(7)
    h0.SetLineWidth(5)
    h0.Fill(each)
print( "fill test sig"  ) 
for each in test_bkg_response: 
    h_test_bkg.SetLineColor(6)
    #h_test_bkg.SetLineColorAlpha(6,0.7)
    h_test_bkg.SetLineWidth(4)
    h_test_bkg.Fill(each)
print( "fill test bkg"  ) 

h0.Scale(1.0/h0.Integral())
h_test_bkg.Scale(1.0/h_test_bkg.Integral())
h_train_sig.Scale(1.0/h_train_sig.Integral())
h_train_bkg.Scale(1.0/h_train_bkg.Integral())
print( "scale")

max_value = max(
    h0.GetMaximum(),
    h_test_bkg.GetMaximum(),
    h_train_sig.GetMaximum(),
    h_train_bkg.GetMaximum(),
)

h0.SetMaximum(1.4 * max_value)
h0.Draw("HIST")
h_test_bkg.Draw("HIST same")
h_train_sig.Draw("HIST same")
h_train_bkg.Draw("HIST same")

leg = ROOT.TLegend()
leg.AddEntry(h0, "Test sig", 'l')
leg.AddEntry(h_test_bkg, "Test bkg", 'l')
leg.AddEntry(h_train_sig, "Train sig", 'f')

leg.AddEntry(h_train_bkg, "Train bkg", 'f')
leg.Draw("same")

canvas.Draw()
canvas.SaveAs(folder+"/responce.png")

