import os
import ROOT
from ROOT import TLorentzVector
from array import array
import numpy as np

from tensorflow import keras
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Input, Bidirectional, Dropout
from tensorflow.keras.layers import LSTM
from tensorflow.keras.utils import Sequence 
from tensorflow.keras.callbacks import ModelCheckpoint

from sklearn.utils.class_weight import compute_class_weight
from pprint import pprint

class CMesonDataset(Sequence):

    def __init__(self, path, batch_size):
        self.root_file = ROOT.TFile(path)
        self.tree = self.root_file.delphys
        self.num_entries = self.tree.GetEntries()
        self.batch_size = batch_size
        
    def __len__(self):
        return int(self.num_entries / float(self.batch_size))

    def __getitem__(self, index):
        start = index * self.batch_size
        end = (index + 1) * self.batch_size
        
        x = []
        y = []

        for entry in range(start, end):
                       
            self.tree.GetEntry(entry)
            dau_pt = np.array(list(self.tree.track_pt), dtype=np.float32)
            dau_deta = np.array(list(self.tree.track_deta), dtype=np.float32)
            dau_dphi = np.array(list(self.tree.track_dphi), dtype=np.float32)
            dau_d0 = np.array(list(self.tree.track_d0), dtype=np.float32)
            dau_dz = np.array(list(self.tree.track_dz), dtype=np.float32)
            dau_xd = np.array(list(self.tree.track_xd), dtype=np.float32)
            dau_yd = np.array(list(self.tree.track_yd), dtype=np.float32)
            dau_zd = np.array(list(self.tree.track_zd), dtype=np.float32)
            
            
            #Sorting
            order_pt = np.argsort(dau_pt)
                                 
            dau_pt = dau_pt[order_pt][::-1]
            dau_deta = dau_deta[order_pt][::-1]
            dau_dphi = dau_dphi[order_pt][::-1]
            dau_d0 = dau_d0[order_pt][::-1]
            dau_dz = dau_dz[order_pt][::-1]
            dau_xd = dau_xd[order_pt][::-1]  
            dau_yd = dau_yd[order_pt][::-1]    
            dau_zd = dau_zd[order_pt][::-1]
            
            #print jet_num, "th Jet's Dauthers are : ", len(dau_pt)
            #print " dau_pt : "
            #print  dau_pt
                  
            dau_set = []
                
            for i in range(0, len(dau_pt)):
                    dau_set.append([])
           
                    dau_set[-1].append(dau_pt[i])
                    dau_set[-1].append(dau_deta[i])
                    dau_set[-1].append(dau_dphi[i])
                    dau_set[-1].append(dau_d0[i])
                    dau_set[-1].append(dau_dz[i])
                    dau_set[-1].append(dau_xd[i])
                    dau_set[-1].append(dau_yd[i])
                    dau_set[-1].append(dau_zd[i])
         
            x.append(dau_set)
            
            if (self.tree.jet_label_d0 ==1): label = 1
            else : label = 0
            
            y.append(label)

        x = keras.preprocessing.sequence.pad_sequences(x, maxlen=15, padding='post', truncating='post', dtype=np.float32)
        y = np.array(y)
        
        return x, y
    
def get_datasets():        
    datasets = [
        CMesonDataset("../2-Selector/reout1.root", batch_size=256),
        CMesonDataset("../2-Selector/reout8.root", batch_size=256),
        CMesonDataset("../2-Selector/reout9.root", batch_size=256),
    ]
    
    trainset, valset, testset = sorted(datasets, key=lambda dset: len(dset), reverse=True)
    return trainset, valset, testset

def main():
    train_set, val_set, test_set = get_datasets()

    print("Train Set : ",trainset, len(trainset) )
    print("Test Set : ",testset, len(testset) )
    print("Val Set : ",valset, len(valset) )
    #print( train_set[0])

if __name__ == '__main__':
    main()
