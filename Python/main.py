### path/to/sklearn-env/Scripts/python main.py  <--- Windows
### path/to/sklearn-env/bin/python3 main.py  <---- Ubuntu

from simData import simData
from sklearn.ensemble import IsolationForest 

def mainIsolationForest():
    sd = simData("../Dat/GTIVTimeSim.dat", n=3, socAri=True, slidingAri=False, nSuteru=2)
    
    print("type(mat)=", type(sd.mat), "mat.shape=", sd.mat.shape)
    print("type(timeStr)=", type(sd.timeStr), "timeStr.shape=", sd.timeStr.shape)
    # for i,vec in enumerate(sd.mat):
    #     print(sd.timeStr[i], vec)
    # exit()####################
    
    iff = IsolationForest(random_state=0).fit_predict(sd.mat)
    print(iff)
    # iff = IsolationForest(random_state=0).fit(sd.mat)
    # df = iff.decision_function(sd.mat)
    # print(df)
    
    for i in range(0, sd.mat.shape[0]):
        print(sd.timeStr[i], iff[i])
        # print(sd.timeStr[i], df[i])

import sklearn.cluster as cl 

def mainKMeans():
    sd = simData("../Dat/GTIVTimeSim.dat", n=5, socAri=False, slidingAri=True)
    
    print("type(mat)=", type(sd.mat), "mat.shape=", sd.mat.shape)
    print("type(timeStr)=", type(sd.timeStr), "timeStr.shape=", sd.timeStr.shape)
    
    kmeans = cl.KMeans(2).fit(sd.mat)
    print(kmeans.labels_)
    
    for i in range(0, sd.mat.shape[0]):
        # print(sd.timeStr[i], sd.mat[i], kmeans.labels_[i])
        print(sd.timeStr[i], kmeans.labels_[i])

# mainKMeans()
mainIsolationForest()
