from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import pandas as pd
import tensorflow as tf
import tempfile
import random
import matplotlib.pyplot as plt
#get_ipython().magic('matplotlib inline')
import csv
import seaborn as sns
from array import array
from sklearn import datasets
from sklearn import metrics
from sklearn import model_selection
from sklearn import preprocessing

import ROOT
ROOT.gROOT.SetBatch(True)
from root_numpy import root2array, tree2array, fill_hist

# Define the algorithm flags
#tf.app.flags.DEFINE_integer("seed", 1, "RNG seed")
#FLAGS = tf.app.flags.FLAGS

#tf.app.flags.DEFINE_float("augment_hue", 0.5, "hue augment factor")
#tf.app.flags.DEFINE_boolean("no_augment", False, "whether to disable augmentation entirely")

##########################
def main(unused_argv):
  # Set TF random seed to improve reproducibility
  seed = 150
  np.random.seed(seed)
  #tf.set_random_seed(FLAGS.seed)
  #print('FLAGS.seed ',FLAGS.seed)

  #with tf.Session(config=tf.ConfigProto(allow_soft_placement = True)) as sess:
  #  if FLAGS.broken:
  #    with tf.device('/cpu:0'):
  #      input = mk_input()
  #  else:
  #    input = mk_input()

  #tf.global_variables_initializer().run()
  #tf.train.start_queue_runners(sess = sess)
  #print('nn: TF session initialized')
  #res = sess.run()
  #print('res: ', res)

  # Load dataset
  #filein = open("data_forRecoLength_0.csv","r")
  #filein = open("data_forRecoLength.csv","r")
  filein = open("data_forRecoLength_9_10MRD.csv","r")
  Dataset = np.array(pd.read_csv(filein))
  #features, labels = np.split(Dataset,[-1],axis=1)
  features, lambdamax, labels = np.split(Dataset,[6003,6004],axis=1)
  #lambdamax = np.zeros(shape=(811,6000))
  #print('labels: ',labels[:2])
  #print('lambdamax: ',lambdamax[:2])
  #print('features: ',features[:2])
 
  #filein2 = open("data_forRecoLength_0.csv","r")
  filein2 = open("data_forRecoLength_0_8MRD.csv","r")
  Dataset2 = np.array(pd.read_csv(filein2))
  #features2, labels2 = np.split(Dataset2,[-1],axis=1) 
  features2, lambdamax2, labels2 = np.split(Dataset2,[6003,6004],axis=1)
  print('labels: ',labels2[0]," last length ", labels2[1678])
  #print('lambdamax: ',lambdamax2[0]," lambdamax2[1]: ",lambdamax2[1]," last lambdamax2[1679] ",lambdamax2[1679])
  print('lambdamax[:2]: ',lambdamax2[:2])
  #print('features: ',features2[:2])  
  print("len(lambdamax2) ", len(lambdamax2))

 # with open('data_forRecoLength.csv') as inf:
 #      reader = csv.reader(inf, delimiter=",")  
 #      lambdamax0 = list(zip(*reader))[6000]
 # with open('data_forRecoLength.csv') as inf0:
 #      reader0 = csv.reader(inf0, delimiter=",")
 #      y0 = list(zip(*reader0))[6003]
  #print('y0: ',y0, ' y0[0] ',y0[1] )  
  #print('lambdamax0: ',lambdamax0)

  num_events, num_pixels = features.shape
  print(num_events, num_pixels)

  #boston = datasets.load_boston()
  #x, y = boston.data, boston.target

  #Split dataset into train / test
  np.random.seed(0)
  #rnd_indices = np.random.rand(len(features)) < 0.80
  train_x = features#[rnd_indices]
  train_y = labels#[rnd_indices]
  #test_x = features[~rnd_indices]
  #test_y = labels[~rnd_indices]
  #lambdamax_test = lambdamax[~rnd_indices]
  test_x = features2
  test_y = labels2
  lambdamax_test = lambdamax2
  print('train: ',len(train_y),' test: ', len(test_y), ' & ', len(lambdamax_test))
  #train_x, test_x, train_y, test_y = model_selection.train_test_split(features, labels, test_size=0.2, random_state=42)
  print("max: ",np.amax(train_x)," shape- ", train_x.shape)
  ## Split dataset into train / test
  #train_x, test_x, train_y, test_y = model_selection.train_test_split(
  #x, y, test_size=0.2, random_state=42)

  #Use dataset from different files
  #train_x =features
  #train_y = labels
  #test_x = features2
  #test_y = labels2
  #lambdamax_test = lambdamax2

  # Scale data (training set) to 0 mean and unit standard deviation.
  scaler = preprocessing.StandardScaler()
  train_x = scaler.fit_transform(train_x)
 
  # Build 2 layer fully connected DNN with 10, 10 units respectively.
  feature_columns = [
      tf.feature_column.numeric_column('x', shape=np.array(train_x).shape[1:])]
  regressor = tf.estimator.DNNRegressor(
     feature_columns=feature_columns, hidden_units=[70, 20])
     #feature_columns=feature_columns, hidden_units=[200, 100, 20])
     #feature_columns=feature_columns, hidden_units=[10, 10])

  # Train.
  print('traiining....')
  batch_size = 1#2
  epochs_no= 2000
  n_batches = int(np.ceil(num_events / batch_size))
  train_input_fn = tf.estimator.inputs.numpy_input_fn(
      x={'x': train_x}, y=train_y, batch_size=batch_size, num_epochs=epochs_no, shuffle=False,num_threads=1)
  regressor.train(input_fn=train_input_fn,steps=1000) #1000)
 # print('before scores..')
 # train2_fn = tf.estimator.inputs.numpy_input_fn( 
 #       x={'x': train_x}, y=train_y, num_epochs=epochs_no, shuffle=False)
 # scoresT = regressor.evaluate(input_fn=train2_fn)
 # for epoch in range(epochs_no):
 #     if epoch % 100 == 0:
 #        print("Epoch number", epoch, "----> Loss:", 'MSE (tensorflow): {0:f}'.format(scoresT['average_loss']))

  # Predict.
  print('predicting...')
  x_transformed = scaler.transform(test_x)
  test_input_fn = tf.estimator.inputs.numpy_input_fn(
      x={'x': x_transformed}, y=test_y, shuffle=False)
  predictions = regressor.predict(input_fn=test_input_fn)
  y_predicted = np.array(list(p['predictions'] for p in predictions))
  y_predicted = y_predicted.reshape(np.array(test_y).shape)

  # Score with sklearn.
  score_sklearn = metrics.mean_squared_error(y_predicted, test_y)
  print('MSE (sklearn): {0:f}'.format(score_sklearn))

  # Score with tensorflow.
  scores = regressor.evaluate(input_fn=test_input_fn)
  print('MSE (tensorflow): {0:f}'.format(scores['average_loss']))

  # pred_y = sess.run(prediction,feed_dict={x:test_x})
  #  mse = tf.reduce_mean(tf.square(pred_y - test_y))
  #  print("MSE: %.4f" % sess.run(mse))
  #  pred_y
        
  # testing
  #correct = tf.equal(tf.argmax(prediction, 1), tf.argmax(y, 1))
  #accuracy = tf.reduce_mean(tf.cast(correct, 'float'))
  #print('Accuracy:', accuracy.eval({x: mnist.test.images, y: mnist.test.labels}))

  fig, ax = plt.subplots()
  ax.scatter(test_y,y_predicted)
  ax.plot([test_y.min(),test_y.max()],[test_y.min(),test_y.max()],'k--',lw=3)
  ax.set_xlabel('Measured')
  ax.set_ylabel('Predicted')
  #plt.show()
  plt.savefig("plotsRecolengthB/test_recolength.png")

  fig, ax = plt.subplots()
  ax.scatter(test_y,y_predicted)
  ax.plot([test_y.min(),test_y.max()],[test_y.min(),test_y.max()],'k--',lw=3)
  ax.set_xlabel('Measured')
  ax.set_ylabel('Predicted')
  ax.set_xlim(-50.,400.)
  ax.set_ylim(-50.,400.)
  #plt.show()
  plt.savefig("plotsRecolengthB/test_recolength_ZOOM.png")

  fig, ax = plt.subplots()
  ax.scatter(test_y,lambdamax_test)
  ax.plot([test_y.min(),test_y.max()],[test_y.min(),test_y.max()],'k--',lw=3)
  ax.set_xlabel('Measured')
  ax.set_ylabel('Predicted')
  ax.set_xlim(-50.,400.)
  ax.set_ylim(-50.,400.)
  #plt.show()
  plt.savefig("plotsRecolengthB/MYrecolength.png")

 # #plt.figure(figsize=(10, 6))
 # #plt.plot(range(len(cost_history)),cost_history, 'b-')
 # #plt.axis([0,epochs_no,0,30000]) #np.max(cost_history)])
 # #plt.plot(range(len(cost_history2)),cost_history2, 'r--')
  #plt.axis([0,epochs_no,0,30000]) #np.max(cost_history)])
  #plt.savefig("plotsRecolengthB/deviation_recolength_WCSim.png")

  #plt.legend(loc='upper right')
  #plt.xlabel('Number of Estimators')
  #plt.ylabel('Least Absolute Deviation [MeV]')
  #plt.savefig("deviation_train_test.png")

  #----------------------------------------
  data = abs(y_predicted-test_y)
  dataprev = abs(lambdamax_test-test_y)
  nbins=np.arange(0,200,10)
  ##n, bins, patches = plt.hist(data, 100, alpha=1,normed='true')
  fig,ax=plt.subplots(ncols=1, sharey=False)#, figsize=(8, 6))
  f0=ax.hist(data, nbins, histtype='step', fill=False, color='blue',alpha=0.75) 
  f1=ax.hist(dataprev, nbins, histtype='step', fill=False, color='red',alpha=0.75)
  ax.set_xlim(0.,200.)
  ax.set_xlabel('$\Delta R$ [cm]')
  ##ax.set_ylabel('Number of Entries [%]')
  ##ax.xaxis.set_label_coords(0.95, -0.08)
  ##ax.yaxis.set_label_coords(-0.1, 0.71)

  #from scipy.stats import norm
  ## Fit a normal distribution to the data:
  ##mu, std = norm.fit(data)

  xmin, xmax = plt.xlim()
  ##x = np.linspace(xmin, xmax, 100)
  ##p = norm.pdf(x, mu, std)
  ##plt.plot(x, p, 'k', linewidth=2)
  title = "mean = %.2f, std = %.2f, Prev: mean = %.2f, std = %.2f " % (data.mean(), data.std(),dataprev.mean(), dataprev.std())
  plt.title(title)
  plt.savefig("plotsRecolengthB/resol_distr2l_WCSim.png")
  ##plt.show()

  #wirte output in .root file
  outputFile = ROOT.TFile.Open("OUTMRD.root", "RECREATE")
  outputTuple = ROOT.TNtuple("tuple", "tuple", "newrecolength")
  #print('type: ', type(y_predicted)) 
  for i in range(len(y_predicted)):
      outputTuple.Fill( y_predicted[i] )
  #    print('y_predicted: ', y_predicted[i])

  outputTuple.Write()
  outputFile.Close()

  rfile = ROOT.TFile('vars_forEreco26_NEW_0_8MRD.root')
  intree = rfile.Get('nu_eneNEW')
  arr_hi_E0 = tree2array(intree)
  trueKE = arr_hi_E0['trueKE']
  neutrinoE = arr_hi_E0['neutrinoE']
  lambda_max_2=arr_hi_E0['lambda_max_2']
  TrueTrackLengthInWater=arr_hi_E0['TrueTrackLengthInWater']
  TrueTrackLengthInMrd=arr_hi_E0['TrueTrackLengthInMrd']
  diffDirAbs2=arr_hi_E0['diffDirAbs2']
  z1b=arr_hi_E0['z1b']
  z2b=arr_hi_E0['z2b']
  z3b=arr_hi_E0['z3b']
  z4b=arr_hi_E0['z4b']
  z5b=arr_hi_E0['z5b']
  z6b=arr_hi_E0['z6b']
  recoDWallR2=arr_hi_E0['recoDWallR2']
  recoDWallZ2=arr_hi_E0['recoDWallZ2']
  totalLAPPDs=arr_hi_E0['totalLAPPDs']
  totalPMTs=arr_hi_E0['totalPMTs']
  vtxX=arr_hi_E0['vtxX']
  vtxY=arr_hi_E0['vtxY']
  vtxZ=arr_hi_E0['vtxZ']
  dirX=arr_hi_E0['dirX']
  dirY=arr_hi_E0['dirY']
  dirZ=arr_hi_E0['dirZ']
  truedirX=arr_hi_E0['truedirX']
  truedirY=arr_hi_E0['truedirY']
  truedirZ=arr_hi_E0['truedirZ']
  truevtxX=arr_hi_E0['truevtxX']
  truevtxY=arr_hi_E0['truevtxY']
  truevtxZ=arr_hi_E0['truevtxZ']  
  trueMomentumTransfer=arr_hi_E0['TrueMomentumTransfer']
  trueMuonAngle=arr_hi_E0['TrueMuonAngle']
 
  #print('lambda_max_2: ',500.*lambda_max_2[0]," lambda_max_2[1] ",500.*lambda_max_2[1]," last lambda_max_2 ",500.*lambda_max_2[1679])
  #print('lambda_max_2[:2]: ',500.*lambda_max_2[:2])
  print("checking..."," len(trueKE): ",len(trueKE)," len(y_predicted): ", len(y_predicted))

  canvas = ROOT.TCanvas()
  canvas.cd(1)
  th2f = ROOT.TH2F("Emu_newrecolength", " ;E_{#mu} [MeV];Reco Track Length [cm]", 100, 0, 2000., 40, 0., 400.)
  for i in range(len(trueKE)):
      th2f.Fill(1000.*trueKE[i], y_predicted[i])  
  th2f.Draw("ColZ") 
  th2f.SetStats(0)
  canvas.Draw()  
  canvas.SaveAs("plotsRecolengthB/Emu_newrecolength.png")

  canvas.cd(1)
  th2f = ROOT.TH2F("True_RecoLength", " ;MC Track Length [cm];Reco Track Length [cm]", 40, 0, 400., 40, 0., 400.)
  for i in range(len(trueKE)):
      th2f.Fill(test_y[i], y_predicted[i])
  line = ROOT.TLine(0.,400.,0.,400.)
  th2f.SetStats(0)
  th2f.Draw("ColZ")
  line.SetLineColor(2)
  line.Draw("same")
  canvas.Draw()
  canvas.SaveAs("plotsRecolengthB/MClength_newrecolength.png")

  canvas.cd(1)
  th2f = ROOT.TH2F("Emu_newrecolength", " ;E_{#mu} [MeV];Reco Track Length [cm]", 100, 0, 2000., 40, 0., 400.)
  for i in range(len(trueKE)):
      th2f.Fill(1000.*trueKE[i], 500.*lambda_max_2[i])
  th2f.Draw("ColZ")
  canvas.Draw()
  canvas.SaveAs("plotsRecolengthB/Emu_lambdamax.png")

  canvas.cd(1)
  th2f = ROOT.TH2F("Emu_newrecolength", " ;E_{#mu} [MeV];MC Track Length [cm]", 100, 0, 2000., 40, 0., 400.)
  for i in range(len(trueKE)):
      th2f.Fill(1000.*trueKE[i], 500.*TrueTrackLengthInWater[i])
  th2f.Draw("ColZ")
  th2f.SetStats(0)
  canvas.Draw()
  canvas.SaveAs("plotsRecolengthB/Emu_MClength.png")

  rfile2 = ROOT.TFile('OUTMRD.root')
  intree2 = rfile2.Get('tuple')
  arr = tree2array(intree2)
  newrecolength=arr['newrecolength']
  print("check: ", len(newrecolength)," | ",len(trueKE))
  outputFile2 = ROOT.TFile.Open("TreeforEnergyRecoB.root", "RECREATE")
  #outputTuple2 = ROOT.TNtuple("tuple", "tuple", "trueKE:neutrinoE:TrueTrackLengthInMrd:diffDirAbs2:z1b:z2b:z3b:z4b:z5b:z6b:recoDWallR2:recoDWallZ2:totalLAPPDs:totalPMTs:newrecolength")#:vtxX:vtxY:vtxZ:dirX:dirY:dirZ:TrueMomentumTransfer:TrueMuonAngle:truedirX:truedirY:truedirZ:truevtxX:truevtxY:truevtxZ")
  outputTuple2 = ROOT.TNtuple("tuple", "tuple", "trueKE:neutrinoE:TrueTrackLengthInMrd:diffDirAbs2:recoDWallR2:recoDWallZ2:totalLAPPDs:totalPMTs:newrecolength:vtxX:vtxY:vtxZ")
  for i in range(len(trueKE)):
     #outputTuple2.Fill(trueKE[i],neutrinoE[i],TrueTrackLengthInMrd[i],diffDirAbs2[i],z1b[i],z2b[i],z3b[i],z4b[i],z5b[i],z6b[i],recoDWallR2[i],recoDWallZ2[i],totalLAPPDs[i],totalPMTs[i],(newrecolength[i]/600.))
    outputTuple2.Fill(trueKE[i],neutrinoE[i],TrueTrackLengthInMrd[i],diffDirAbs2[i],recoDWallR2[i],recoDWallZ2[i],totalLAPPDs[i],totalPMTs[i],(newrecolength[i]/600.),vtxX[i]/150.,vtxY[i]/200.,vtxZ[i]/150.)
  outputTuple2.Write()
  outputFile2.Close()
 
  outputFile3 = ROOT.TFile.Open("TreeforMCInfoB.root", "RECREATE")
  outputTuple3 = ROOT.TNtuple("tuple", "tuple", "trueKE:neutrinoE:vtxX:vtxY:vtxZ:dirX:dirY:dirZ:TrueMomentumTransfer:TrueMuonAngle:truedirX:truedirY:truedirZ:truevtxX:truevtxY")#:truevtxZ")
  for i in range(len(trueKE)):
      outputTuple3.Fill(trueKE[i],neutrinoE[i],vtxX[i],vtxY[i],vtxZ[i],dirX[i],dirY[i],dirZ[i],trueMomentumTransfer[i],trueMuonAngle[i],truedirX[i],truedirY[i],truedirZ[i],truevtxX[i],truevtxY[i])
  outputTuple3.Write()
  outputFile3.Close()

if __name__ == '__main__':
   #flags = tf.app.flags
   #flags.DEFINE_integer("random_seed", 1, "Value of random seed")
   #tf.flags.DEFINE_integer("tf_random_seed", 1,
   #                     """Random seed for TensorFlow initializers. Setting
   #                     this value allows consistency between reruns.""")
   #tf_random_seed=1
   tf.logging.set_verbosity(tf.logging.INFO)
   tf.app.run(main=main)


