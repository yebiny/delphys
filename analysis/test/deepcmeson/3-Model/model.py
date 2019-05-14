import os
os.environ['CUDA_VISIBLE_DEVICES'] = "0"

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
from dataset import get_datasets


def build_model():
    model = Sequential()
    model.add(Bidirectional(LSTM(64, return_sequences=True), input_shape=(15, 8)) )
    model.add(Bidirectional(LSTM(128)))
    model.add(Dropout(0.5))
    model.add(Dense(1, activation='sigmoid'))
    # try using different optimizers and different optimizer configs
    model.compile('adam', 'binary_crossentropy', metrics=['accuracy'])
    return model

    
def getyinfo(model, xset):
    y_true = []
    y_score = []
    
    for idx in range(len(xset)):
        x, y = xset[idx]
        y_predict = model.predict_on_batch(x)
        
        y_true.append(y)
        y_score.append(y_predict)
   
    y_true = np.concatenate(y_true)
    y_score = np.concatenate(y_score)
    
    return y_true, y_score

def evaluate(model, y_vloss, y_loss,  train_set, test_set, path):
    train_y_true, train_y_score = getyinfo(model, train_set)
    test_y_true, test_y_score = getyinfo(model, test_set)

    
    train_is_sig = train_y_true.astype(np.bool)
    train_is_bkg = np.logical_not(train_is_sig)
    test_is_sig = test_y_true.astype(np.bool)
    test_is_bkg = np.logical_not(test_is_sig)

    train_sig_response = train_y_score[train_is_sig]
    train_bkg_response = train_y_score[train_is_bkg]
    test_sig_response = test_y_score[test_is_sig]
    test_bkg_response = test_y_score[test_is_bkg] 

    np.savez(path,
        y_vloss=y_vloss,
        y_loss=y_loss,
        train_sig_response=train_sig_response,
        train_bkg_response=train_bkg_response,
        test_sig_response=test_sig_response,
        test_bkg_response=test_bkg_response,
        )
def main():
    train_set, val_set, test_set = get_datasets()
    model = build_model()
    
    checkpointer = ModelCheckpoint(filepath='./weights.hdf5', verbose=1, save_best_only=True)

    # training
    history = model.fit_generator(
    generator=train_set,
    validation_data=val_set,
    steps_per_epoch=len(train_set), 
    epochs=20,
    callbacks=[checkpointer]
    )
    y_vloss = history.history['val_loss']
    y_loss = history.history['loss']


    # evaluation
    evaluate(model, y_vloss, y_loss, train_set, test_set, './model3_Info.npz')
    keras.utils.plot_model(model, to_file='./model3.png', show_shapes=True, show_layer_names=True)

if __name__ == '__main__':
    main()



