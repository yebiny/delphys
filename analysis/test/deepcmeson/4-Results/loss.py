from __future__ import division
from collections import OrderedDict
import numpy as np
from ROOT import TH1F
import matplotlib.pyplot as plt
import os,sys

folder = './'+sys.argv[1:][0]

info = np.load(folder+"/model_info.npz")
vloss =  info['y_vloss']
loss = info['y_loss']

x_len = np.arange(len(loss))
plt.plot(x_len, vloss, marker='.', c = 'red', label = 'Valdation-set Loss')
plt.plot(x_len, loss, marker='.', c = 'blue', label = 'Train-Set Loss')
plt.legend(fontsize=15)
plt.grid()
plt.savefig(folder+"/loss.png")
