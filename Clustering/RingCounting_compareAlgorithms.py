#https://github.com/Rockyzsu/base_function/blob/master/sklearn_basic.py
from __future__ import unicode_literals
import datetime
from collections import Counter

import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from sklearn import cluster
from sklearn import datasets


#iris = datasets.load_iris()
## X_iris = iris.data[50: 100]
#X_iris = iris.data
#Y_iris = iris.target
geo = 231
#print("X_iris: ",X_iris)
#print("Y_iris: ",Y_iris)
# X_iris = np.delete(X_iris, 3, axis=1)
# X_iris /= 10.

#Read data
infile = "data/EventTableRings.csv"
filein = open(str(infile))
print("data in: ",filein)
df0=pd.read_csv(filein)
print(df0.head())
booldf = df0["Radius in m"].isnull()
null_columns=df0.columns[df0.isnull().any()]
nulldf = df0[df0["Radius in m"].isnull()][null_columns]
#print(nulldf)
print("nulldf.index.values: ",nulldf.index.values)
#print(df0[140:145])
count = 0
for g, dfs in nulldf.groupby(np.arange(len(nulldf)) // 3):
    print("----- g: ",g," dfs: ",dfs)
    a = dfs.index.values
    #print("null indices : ",a)
    if g==0: 
       df_evt = df0.loc[:a[0]-1,:]
       count = a[0]+3
    else:
       #print("count: ",count)
       #print("from index: ", count , " to index: ", a[0]-1)
       df_evt = df0.loc[count:a[0]-1,:]
       count = a[0]+3
    print("df_evt.head(): ",df_evt.head())

    #drop LAPPDs: 
    #indexNames = df_evt[df_evt['Detector type']=='LAPPD'].index
    #df_evt.drop(indexNames , inplace=True)
    #print("PMTS only-df_evt.head(): ",df_evt.head())
 
    #drop some columns from the data sample
    X  = df_evt.drop(["Detector type", "Ring number", "Radius in m","Angle in rad","Height in m"], axis=1).values
    #print(X)
    print("X.shape: ",X.shape," type(X)",type(X)," len(X) ", len(X))
    #store true labels
    labels_true = df_evt["Ring number"].values
    print("labels_true.shape: ",labels_true.shape," type(labels_true): ",type(labels_true))
    print("Number of Rings: ",labels_true[0])
    Y = labels_true[0]

    # #############################################################################
    def timeit(name=None):
        """
        @decorate
        """
        def wrapper2(func):
            def wrapper1(*args, **kargs):
                start = datetime.datetime.now()
                r = func(*args, **kargs)
                end = datetime.datetime.now()
                print('---------------')
                print('project name: %s' % name)
                print('start at: %s' % start)
                print('end at:   %s' % end)
                print('cost:     %s' % (end - start))
                print('res:      %s' % r)
#            print('err1:     %s' \)
            #      % (50 - Counter(r[: 50]).most_common()[0][1])
#            print('err2:     %s' \)
            #      % (50 - Counter(r[50: 100]).most_common()[0][1])
#            print('err3:     %s' \)
            #      % (50 - Counter(r[100: 150]).most_common()[0][1])
                print('---------------')
                return r
            return wrapper1
        return wrapper2


    def randrange(n, vmin, vmax):
        return (vmax - vmin) * np.random.rand(n) + vmin


    @timeit('target')
    def target(fig):
        global X, Y, geo
        ax = fig.add_subplot(geo + 0, projection='3d', title='target')
        for n, i in enumerate(X):
            #print("n: ", n," i: ",i)
#            ax.scatter(*i[: 3], c=['r', 'y', 'g'][Y[n]], marker='o')
             ax.scatter(*i[: 3], marker='o')

        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.title('True number of clusters: %d' % labels_true[0])
        return Y


    # kmeans
    @timeit('kmeans')
    def kmeans(fig):
        global X, geo
        ax = fig.add_subplot(geo + 1, projection='3d', title='k-means')
        k_means = cluster.KMeans(init='random', n_clusters=3)
        k_means.fit(X)
        res = k_means.labels_
        n_clusters_ = len(set(res)) - (1 if -1 in res else 0)
        n_noise_ = list(res).count(-1)
        print('in kmenas - Estimated number of clusters: %d' % n_clusters_)
        print('in kmeans - Estimated number of noise points: %d' % n_noise_)
        for n, i in enumerate(X):
            ax.scatter(*i[: 3], c='bgrcmyk'[res[n] % 7], marker='o')

        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.title('k-means: Number of clusters: %d' % n_clusters_)
        return res


    @timeit('mini_batch_kmeans')
    def mini_batch(fig):
        global X, geo
        ax = fig.add_subplot(geo + 2, projection='3d', title='mini-batch')
        mini_batch = cluster.MiniBatchKMeans(init='random', n_clusters=3)
        mini_batch.fit(X)
        res = mini_batch.labels_
        n_clusters_ = len(set(res)) - (1 if -1 in res else 0)
        n_noise_ = list(res).count(-1)
        print('in mini_batch_kmeans - Estimated number of clusters: %d' % n_clusters_)
        print('in mini_batch_kmeans - Estimated number of noise points: %d' % n_noise_)
        for n, i in enumerate(X):
            ax.scatter(*i[: 3], c='bgrcmyk'[res[n] % 7], marker='o')

        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.title('mini-batch-kmeans: Number of clusters: %d' % n_clusters_)
        return res


    @timeit('affinity')
    def affinity(fig):
        global X, geo
        ax = fig.add_subplot(geo + 3, projection='3d', title='affinity')
        affinity = cluster.AffinityPropagation(preference=-50)
        affinity.fit(X)
        res = affinity.labels_
        n_clusters_ = len(set(res)) - (1 if -1 in res else 0)
        n_noise_ = list(res).count(-1)
        print('in affinity - Estimated number of clusters: %d' % n_clusters_)
        print('in affinity - Estimated number of noise points: %d' % n_noise_)
        for n, i in enumerate(X):
            ax.scatter(*i[: 3], c='bgrcmyk'[res[n] % 7], marker='o')

        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.title('affinity: Number of clusters: %d' % n_clusters_)
        return res


    @timeit('mean_shift')
    def mean_shift(fig):
        global X, geo
        ax = fig.add_subplot(geo + 4, projection='3d', title='mean_shift')
        bandwidth = cluster.estimate_bandwidth(X, quantile=0.2, n_samples=50)
        mean_shift = cluster.MeanShift(bandwidth=bandwidth, bin_seeding=True)
        mean_shift.fit(X)
        res = mean_shift.labels_
        n_clusters_ = len(set(res)) - (1 if -1 in res else 0)
        n_noise_ = list(res).count(-1)
        print('in mean_shift - Estimated number of clusters: %d' % n_clusters_)
        print('in mean_shift - Estimated number of noise points: %d' % n_noise_)
        for n, i in enumerate(X):
            ax.scatter(*i[: 3], c='bgrcmyk'[res[n] % 7], marker='o')

        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.title('mean_shift: Number of clusters: %d' % n_clusters_)
        return res


    @timeit('dbscan')
    def dbscan(fig):
        global X, geo
        ax = fig.add_subplot(geo + 5, projection='3d', title='dbscan')
        dbscan = cluster.DBSCAN()
        dbscan.fit(X)
        res = dbscan.labels_
        core = dbscan.core_sample_indices_
        print(repr(core))
        size = [5 if i not in core else 40 for i in range(len(X))]
        print(repr(size))
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(res)) - (1 if -1 in res else 0)
        n_noise_ = list(res).count(-1)
        print('in dbscan - Estimated number of clusters: %d' % n_clusters_)
        print('in dbscan - Estimated number of noise points: %d' % n_noise_)
        for n, i in enumerate(X):
            #print ("in dbscan -n: ",n," i: ",i)
            ax.scatter(*i[: 3], s=size[n], c='bgrcmyk'[res[n] % 7],
                       alpha=0.8, marker='o')

        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.title('dbscan: Number of clusters: %d' % n_clusters_)
        return res


    def main():
        print("I'm in g: ", g)
        if g<5:
            fig = plt.figure(figsize=(20,10))
            target(fig)
            kmeans(fig)
            mini_batch(fig)
            affinity(fig)
            mean_shift(fig)
            dbscan(fig)

#            plt.show()
            plt.savefig("plots/VisRings"+str(g)+".png")

    if __name__ == '__main__':
       main()
       #print(" in main part ")
