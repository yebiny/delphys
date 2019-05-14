import ROOT

trainset = "train_set9.root"
rootFile = ROOT.TFile("../3-Dataset/"+trainset)
#rootFile = ROOT.TFile("out.root")
tree = rootFile.Get("delphys")

label = [entry.jet_label_d0 for entry in tree]

num_sig = tree.Project("sig_jet_label_d0", 'jet_label_d0', 'jet_label_d0 == 1')
num_bkg = tree.Project("bkg_jet_label_d0", 'jet_label_d0', 'jet_label_d0 != 1')


print trainset
print "SIG:",num_sig
print "BKG:", num_bkg
