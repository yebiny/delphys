import os
from collections import OrderedDict
import ROOT
from ROOT import gRandom
from ROOT import gROOT


def get_name(path):
    basename = os.path.basename(path)
    root, _ = os.path.splitext(basename)
    return root


def get_random_color():
    r, g, b = [gRandom.Uniform(0, 1) for _ in range(3)]
    return ROOT.TColor.GetColor(r, g, b)


class Factory(object):
    def __init__(self, paths, tree_name):
        self._files = OrderedDict({get_name(each): ROOT.TFile.Open(each, "READ") for each in paths})
        self._names = self._files.keys()
        self._trees = OrderedDict({name: self._files[name].Get(tree_name) for name in self._names})
        self._colors = {name: get_random_color() for name in self._names}

        self._memory = OrderedDict()

    def __getitem__(self, key):
        return self._memory[key]

    def __iter__(self):
        return self._trees.iteritems()


    def book(self, hist_name, title, nbinsx, xup, xlow, **kwargs):
        canvas = ROOT.TCanvas("c_" + hist_name, hist_name, 1200, 800)
        bkg = ROOT.TH1F(hist_name, title, nbinsx, xup, xlow)

        hists = {}
        for tree_name in self._trees.keys():
            name = tree_name + "_" + hist_name
            hists[tree_name] = ROOT.TH1F(name, title, nbinsx, xup, xlow)

        setattr(self, hist_name, hists)

        self._memory[hist_name] = {
            "canvas": canvas,
            "bkg": bkg,
            "histograms": hists
        }

    def draw(self):
        print("Drawing")
        for hist_name in self._memory.keys():
            bkg = self._memory[hist_name]["bkg"]
            histograms = self._memory[hist_name]["histograms"]
            canvas = self._memory[hist_name]["canvas"]
            canvas.cd()

            legend = ROOT.TLegend(400, 200)
            for tree_name in self._trees.keys():
                hist = self._memory[hist_name]["histograms"][tree_name]
                hist.SetLineColor(self._colors[tree_name])
                hist.SetLineWidth(4)
                hist.Scale(1.0 / hist.GetEntries())
                legend.AddEntry(hist, tree_name, "l")

            max_value = max(each.GetMaximum() for each in self._memory[hist_name]["histograms"].values()) 
            bkg.SetMaximum(max_value * 1.1)

            bkg.Draw("hist")
            legend.Draw()

            for t_name in self._trees.keys():
                hist = histograms[t_name]
                hist.Draw("hist same")

            canvas.SaveAs("{}.png".format(hist_name))
            self._memory.pop(hist_name)


if __name__ == "__main__":
    from utils import load_delphes
    load_delphes("/cms/ldap_home/slowmoyang/Install/Delphes-3.4.1")

    gROOT.SetBatch(True)
    seed = ROOT.TDatime().Get()
    gRandom.SetSeed(seed)

    datasets = [
        "QCD_HT300to500",
        "QCD_HT500to700",
        "QCD_HT2000toInf",
        "TT-alljets"
    ]
    path_template = "/xrootd/store/user/seyang/Data/TopTagging/{0}/{0}_1.root"
    paths = [path_template.format(each) for each in datasets]

    factory = Factory(paths, "Delphes")

    factory.book("scalar_ht", "Scalar sum of transverse momenta", 100, 0, 4000)
    factory.book("leading_pt", "Leading p_{T}", 100, 0, 2000)
    factory.book("six_jets_pt_sum", "\sum_{j=1}^{6} p_{T}^j", 100, 0, 4000)
    factory.book("num_b_tagged_jets", "# of b-tagged jets", 10, -0.5, 9.5)
    factory.book("num_electrons", "# of electrons", 10, -0.5, 9.5)
    factory.book("num_muons", "# of muons", 10, -0.5, 9.5)

    for name, tree in factory:
        print(name)
        for entry in tree:
            if entry.Jet_size < 6:
                continue
            factory.scalar_ht[name].Fill(entry.ScalarHT.At(0).HT)
            factory.leading_pt[name].Fill(entry.Jet.At(0).PT)

            six_jets_pt_sum = sum(entry.Jet.At(i).PT for i in range(6))
            factory.six_jets_pt_sum[name].Fill(six_jets_pt_sum)

            factory.num_b_tagged_jets[name].Fill(sum(each.BTag for each in entry.Jet))

            factory.num_electrons[name].Fill(entry.Electron.GetEntries())
            factory.num_muons[name].Fill(entry.Muon.GetEntries())

    factory.draw()
