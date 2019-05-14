import ROOT

root_file = ROOT.TFile("reout8.root")
tree = root_file.delphys

d0 = tree.AsMatrix(["jet_label_d0"])

print(d0[:100])
print(d0[-100:])
