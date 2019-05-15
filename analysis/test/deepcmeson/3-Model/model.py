import os
os.environ['CUDA_VISIBLE_DEVICES'] = "0"

import ROOT, sys
from ROOT import TLorentzVector
from array import array
import numpy as np

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Input, Bidirectional, Dropout

if tf.test.is_gpu_available(cuda_only=True):
    from tensorflow.keras.layers import CuDNNLSTM as LSTM
else:
    from tensorflow.keras.layers import LSTM
from tensorflow.keras.utils import Sequence 
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.callbacks import ModelCheckpoint


from sklearn.utils.class_weight import compute_class_weight
from pprint import pprint
from dataset import get_datasets


''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''                                                '''
'''   python3 model.py [title] [number of epochs]  '''
'''                                                '''
''''''''''''''''''''''''''''''''''''''''''''''''''''''

def build_model(x_shape):
    
    model = Sequential()
    model.add(Bidirectional(LSTM(64, return_sequences=True), input_shape=(x_shape[1],x_shape[2])))
    model.add(Bidirectional(LSTM(128)))
    model.add(Dropout(0.5))
    model.add(Dense(1, activation='sigmoid'))
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

def evaluate(model,train_set, test_set):
    
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

    return train_sig_response, train_bkg_response, test_sig_response, test_bkg_response

def main():

    # set name
    modelname = 'test'
    epochs = 1
    
    modelname = sys.argv[1:][0]
    epochs = int(sys.argv[1:][1])

    # set save path
    savepath = '../4-Results/'+modelname+'/'
    if os.path.isdir(savepath):
        print("Already exist. Continue?")    
        a = input()
        if a != 'y': 
            sys.exit()
    if not os.path.isdir(savepath):
        os.mkdir(savepath)

    # set datasets
    train_set, val_set, test_set = get_datasets()
    tmp_x, tmp_y = train_set[0]
    x_shape = tmp_x.shape

    # set model
    model = build_model(x_shape)
    
    # set checkpointer
    checkpointer = ModelCheckpoint(filepath=savepath+'weights.hdf5', verbose=1, save_best_only=True)

    # training
    history = model.fit_generator(
        generator = train_set,
        validation_data = val_set,
        steps_per_epoch = len(train_set), 
        epochs = epochs,
        callbacks = [checkpointer]
    )
    
    # evaluation
    train_s_res, train_b_res, test_s_res, test_b_res = evaluate(model, train_set, test_set)
    y_vloss = history.history['val_loss']
    y_loss = history.history['loss']
   
    # save model and result informations
    keras.utils.plot_model(model, to_file=savepath+'model_plot.png', show_shapes=True, show_layer_names=True)
    np.savez(
        savepath+'model_info.npz',
        y_vloss = y_vloss,
        y_loss = y_loss,
        train_sig_response = train_s_res,
        train_bkg_response = train_b_res,
        test_sig_response = test_s_res,
        test_bkg_response = test_b_res,
        train_len = len(train_set),
        val_len = len(val_set),
        test_len = len(test_set)
    )

if __name__ == '__main__':
    main()
