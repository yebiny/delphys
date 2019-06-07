import ROOT,sys

# Set Input Folder 
data_name = 'test'
if len(sys.argv) == 2:
    data_name = sys.argv[1:][0]
elif len(sys.argv) > 2 : 
    print "Wrong Input"
    sys.exit()

file_name = "../2-Analyser/result_"+data_name+"/deepc_test.root" 
rootFile = ROOT.TFile(file_name)
tree = rootFile.Get("delphys")

label = [entry.jet_label for entry in tree]
num_sig = tree.Project("sig_jet_label", 'jet_label', 'jet_label == 3')
num_bkg = tree.Project("bkg_jet_label", 'jet_label', 'jet_label != 3')

print file_name
print "SIG:",num_sig
print "BKG:", num_bkg
