from keras.models import Sequential
from keras.layers.core import Dense, Activation
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=80)
from sklearn.model_selection import train_test_split
from keras.layers import Dropout
from keras.constraints import maxnorm
from keras.optimizers import Nadam

import ROOT
from root_numpy import root2array, tree2array, fill_hist

model = Sequential()
#model.add(Dense(200, input_dim=5, init='he_uniform', activation='relu'))
model.add(Dense(190, input_dim=5, init='he_uniform', activation='relu'))#best..
model.add(Dropout(0.0))
model.add(Dense(1, init='he_uniform', activation='relu'))
#model.add(Dense(1, activation='relu'))
#model.add(Dense(1))

nadam = Nadam(lr=0.1)#, momentum=0.8)
model.compile(loss='mse',optimizer=nadam)  #'Nadam')

# fix random seed for reproducibility
#seed = 500 #100 #G 100 epochs
seed = 700 #G
#seed = 5000 #F
#seed = 10000 #E
#seed = 1000 #D
#seed = 110 #C
#seed = 50 #B
np.random.seed(seed)

columns = ['total_hits2', 'total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2','trueKE']
data = root2array('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root', branches=columns, treename='nu_eneNEW')

data_input = data[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]
data_truth = data[['trueKE']]
data_input_n = data_input.view(data_input.dtype[0]).reshape(data_input.shape + (-1,))
data_truth_n = data_truth.view(data_truth.dtype[0]).reshape(data_truth.shape + (-1,))#/1e3
print type(data_input_n), data_input_n.shape
print data_input_n[:2]
print data_truth_n[:2]

print ' KE: ' , data_truth_n.shape
x_train, x_test, y_train, y_test = train_test_split(data_input_n, data_truth_n, test_size=0.10, random_state=42)

#history = model.fit(data_input_n, data_truth_n, validation_split = 0.1, nb_epoch=100)
history = model.fit(x_train, y_train, nb_epoch=100, batch_size=10, validation_data=(x_test, y_test)) 
### visualise the model:
#from keras.utils.visualize_util import plot
#plot(model, to_file='model.png')

# list all data in history
print(history.history.keys())
# summarize history for accuracy
f, ax2 = plt.subplots(1,1)
# summarize history for loss
ax2.plot(history.history['loss'])
ax2.plot(history.history['val_loss'])
ax2.set_title('model loss')
ax2.set_ylabel('Performance')
ax2.set_xlabel('Epochs')
ax2.legend(['train', 'test'], loc='upper left')
plt.savefig("keras_train_test.pdf")
#######

predicted_energy = model.predict(data_input_n)
print 'predicted energy: ', predicted_energy[:2]
chain = ROOT.TChain('nu_eneNEW')
for i in range(1040,1099):
    chain.Add('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_'+str(i)+'_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
test_data = tree2array(chain, branches=columns)

print 'len(test_data) ',len(test_data)

test_data_input = test_data[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]
test_data_truth = test_data[['trueKE']]
test_data_input_n = test_data_input.view(test_data_input.dtype[0]).reshape(test_data_input.shape + (-1,))
test_data_truth_n = test_data_truth.view(test_data_truth.dtype[0]).reshape(test_data_truth.shape + (-1,))#/1e3
test_predicted_energy = model.predict(test_data_input_n)
test_predicted_energy.shape

test_data_hi_E = tree2array(chain, selection='trueKE')#>'+str(E_threshold))
test_data_hi_Eneu = tree2array(chain, selection='neutrinoE')
recoDwall = tree2array(chain, selection='recoDWall_2')
recoToWall = tree2array(chain, selection='recoToWall_2')
#float recoDWall_2=recoDWall/550.;   float recoToWall_2=recoToWall/2500.;
test_data_trueKE_hi_E = test_data_hi_E['trueKE']
test_data_neutrinoE_hi_E = test_data_hi_Eneu['neutrinoE']
recoDwall_data = recoDwall['recoDWall_2']
recoToWall_data = recoToWall['recoToWall_2']
test_data_recoKE_hi_E = model.predict(test_data_input_n)

##### store in file:
#outputFile = ROOT.TFile.Open("keras_reco_TrainTestSample100ep.root", "RECREATE")
outputFile = ROOT.TFile.Open("keras_reco_TrainTestSample.root", "RECREATE")
outputTuple = ROOT.TNtuple("tuple", "tuple", "trueKE:neutrinoE:recoKE:recoDwall:recoToWall")
for i in range(len(test_data_recoKE_hi_E)):
     outputTuple.Fill(test_data_trueKE_hi_E[i], test_data_neutrinoE_hi_E[i],test_data_recoKE_hi_E[i], recoDwall_data[i], recoToWall_data[i])
outputTuple.Write()
outputFile.Close()

