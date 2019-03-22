import sys
sys.path.insert(0,'/usr/local/lib/python2.7/site-packages/')
import numpy as np
import matplotlib.pyplot as plt
from keras.models import Sequential
from keras.layers.core import Dense, Activation
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn import linear_model, ensemble, grid_search, metrics
#from hyperas import optim
#from hyperas.distributions import choice, uniform, conditional
import ROOT
ROOT.gROOT.SetBatch(True)
from root_numpy import root2array, tree2array, fill_hist
import os

chain = ROOT.TChain('nu_eneNEW')
chain.Add('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
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

#E_threshold_lo = 200
#E_threshold_hi = 600
#select='trueKE>'+str(E_threshold_lo)+'&&'+'trueKE<'+str(E_threshold_hi)
#data = tree2array(chain,selection=select)
#data_reduced   = data[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
#data_reduced_n = data_reduced.view(data_reduced.dtype[0]).reshape(data_reduced.shape + (-1,))
#data_target    = data['trueKE']/1e3

#params = {'n_estimators': 1000, 'max_depth': 10, 'learning_rate': 0.01, 'loss': 'lad'}
#net_hi_E = ensemble.GradientBoostingRegressor(**params)

#scoresBDTG = cross_val_score(net_hi_E, data_reduced_n, data_target, cv=5, scoring='neg_mean_squared_error')
#print scoresBDTG
#print("BDTG MSE: %0.3f (+/- %0.3f)" % (scoresBDTG.mean(), scoresBDTG.std() * 2))

# Training the ML algo
#clf = linear_model.SGDRegressor()
#def base_model0():
#    model0 = Sequential()
#    model0.add(Dense(1000, input_dim=6, init='uniform', activation='tanh'))
#    model0.add(Dense(1000, activation='tanh'))
#    model0.add(Dense(1000, activation='relu'))
#    model0.add(Dense(1000, activation='relu'))
#    model0.add(Dense(100, activation='relu'))
#    model0.add(Dense(100, activation='relu'))
#    model0.add(Dense(100, activation='relu'))
#    model0.add(Dense(10, activation='relu'))
#    model0.add(Dense(10, activation='relu'))
#    model0.add(Dense(10, activation='relu'))
#    model0.add(Dense(1))
#    model0.compile(loss='mse',optimizer='adam')
#    return model0

def base_model():
    model = Sequential()
    model.add(Dense(1000, input_dim=5, init='uniform', activation='tanh'))
    model.add(Dense(1000, activation='tanh'))
    model.add(Dense(1000, activation='relu'))
    model.add(Dense(1000, activation='relu'))
    model.add(Dense(100, activation='relu'))
    model.add(Dense(100, activation='relu'))
    model.add(Dense(100, activation='relu'))
#    model.add(Dense(15, activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(1, activation='relu'))
    model.add(Dense(1))
    model.compile(loss='mse',optimizer='adam')
    return model

#history = base_model().fit(data_reduced_n, data_target, validation_split = 0.1, nb_epoch=200)

#model0 = KerasRegressor(build_fn=base_model, nb_epoch=100, batch_size=2000, verbose=0)
model0 = KerasRegressor(build_fn=base_model, nb_epoch=100, verbose=0)
scores0 = cross_val_score(model0, data_reduced_n, data_target, cv=5, scoring='neg_mean_squared_error')
print scores0
print("MSE: %0.3f (+/- %0.3f)" % (scores0.mean(), scores0.std() * 2))
print '----------------'
model = KerasRegressor(build_fn=base_model, nb_epoch=200, verbose=0)
scores = cross_val_score(model, data_reduced_n, data_target, cv=5, scoring='neg_mean_squared_error')
print scores
print("MSE: %0.3f (+/- %0.3f)" % (scores.mean(), scores.std() * 2))

model1 = KerasRegressor(build_fn=base_model, nb_epoch=500, verbose=0)
scores1 = cross_val_score(model1, data_reduced_n, data_target, cv=5, scoring='neg_mean_squared_error')
print scores1
print("MSE: %0.3f (+/- %0.3f)" % (scores1.mean(), scores1.std() * 2))

model2 = KerasRegressor(build_fn=base_model, nb_epoch=1000, verbose=0)
scores2 = cross_val_score(model2, data_reduced_n, data_target, cv=5, scoring='neg_mean_squared_error')
print scores2
print("MSE: %0.3f (+/- %0.3f)" % (scores2.mean(), scores2.std() * 2))

print '________________'




