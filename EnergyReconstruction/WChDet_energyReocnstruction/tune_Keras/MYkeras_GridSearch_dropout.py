from keras.models import Sequential
from keras.layers.core import Dense, Activation
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=80)
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from keras.wrappers.scikit_learn import KerasRegressor
from keras.layers import Dropout
from keras.constraints import maxnorm

import ROOT
from root_numpy import root2array, tree2array, fill_hist
import os

def create_model(dropout_rate=0.0, weight_constraint=0):
    model = Sequential()
    model.add(Dense(10, input_dim=5, init='he_uniform', activation='relu'))
    model.add(Dropout(dropout_rate))
    model.add(Dense(1, init='he_uniform', activation='relu'))
    model.compile(loss='mse', optimizer='Nadam', metrics=['accuracy'])
    #model.compile(loss='mse',optimizer=optimizer, metrics=['accuracy'])
    return model

##history = model.fit(data_input_n, data_truth_n, validation_split = 0.1, nb_epoch=100)
#history = model.fit(x_train, y_train, nb_epoch=50, validation_data=(x_test, y_test)) 
chain = ROOT.TChain('nu_eneNEW')
for i in range(1040,1099):
    input_file = '/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_'+str(i)+'_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root'
    if os.path.exists(input_file):
        chain.Add(input_file)
E_threshold_lo = 200
E_threshold_hi = 600
select='trueKE>'+str(E_threshold_lo)+'&&'+'trueKE<'+str(E_threshold_hi)
data = tree2array(chain,selection=select)
#data = tree2array(chain)
data_reduced   = data[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]
data_reduced_n = data_reduced.view(data_reduced.dtype[0]).reshape(data_reduced.shape + (-1,))
data_target    = data['trueKE']/1e3

model = KerasRegressor(build_fn=create_model, nb_epoch=100, batch_size=10, verbose=0)
#batch_size = [10, 20, 40, 60, 80, 100]
#nb_epoch = [10, 50, 100]
#optimizer = ['SGD', 'RMSprop', 'Adagrad', 'Adadelta', 'Adam', 'Adamax', 'Nadam']
#learn_rate = [0.001, 0.01, 0.1, 0.2, 0.3]
#momentum = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9]
#init_mode = ['uniform', 'lecun_uniform', 'normal', 'zero', 'glorot_normal', 'glorot_uniform', 'he_normal', 'he_uniform']
#activation = ['softmax', 'softplus', 'softsign', 'relu', 'tanh', 'sigmoid', 'hard_sigmoid', 'linear']
dropout_rate = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
weight_constraint = [1, 2, 3, 4, 5]
#neurons = [1, 5, 10, 15, 20, 25, 30]

param_grid = dict(dropout_rate=dropout_rate, weight_constraint=weight_constraint)

grid = GridSearchCV(model, param_grid, scoring='neg_mean_squared_error', n_jobs=-1)
grid_result = grid.fit(data_reduced_n, data_target)
# summarize results
print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
means = grid_result.cv_results_['mean_test_score']
stds = grid_result.cv_results_['std_test_score']
params = grid_result.cv_results_['params']
for mean, stdev, param in zip(means, stds, params):
    print("%f (%f) with: %r" % (mean, stdev, param))

#######
#Best: -0.002352 using {'nb_epoch': 100, 'batch_size': 10}
#Best: -0.002464 using {'learn_rate': 0.1, 'momentum': 0.8}
#Best: -0.001723 using {'dropout_rate': 0.0, 'weight_constraint': 5}
#######


