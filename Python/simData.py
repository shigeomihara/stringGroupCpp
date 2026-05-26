### path/to/sklearn-env/Scripts/python main.py  <--- Windows
### path/to/sklearn-env/bin/python3 main.py  <---- Ubuntu

import numpy as np
import csv

class simData:
    def __init__(self, fileName, n=5, socAri=True, slidingAri=True, nSuteru=0):
        self.fileName = fileName
        self.n = n   # number of times of a sliding window
        # self.n = 1   # number of times of a sliding window
        # self.socAri = False
        self.socAri = socAri
        self.slidingAri = slidingAri

        if self.socAri == False:
            self.matDim1Size = self.n*4
        else:
            # self.matDim1Size = self.n*5  # for SOC
            self.matDim1Size = self.n  # for SOC#########################
            
        self.mat = np.zeros((0, self.matDim1Size), dtype=np.float64)
        self.timeStr = np.zeros((0,))
        
        np.set_printoptions(linewidth=75)
        np.set_printoptions(threshold=np.inf)

        if self.slidingAri == True:
            self.readCsvFile(self.fileName)
        else:
            self.readCsvFileNoSliding(self.fileName, nSuteru)
            
        self.mat, self.timeStr = self.preprocessing()

    # return system output coefficient, G W/m^2, T degC, P W
    def getSOC(self, G, T, P):
        return P/(G/1000.0*(1.0-0.0044*(T-25.0)))

    def saveSOCs(self, csvReader):
        f = open("../Dat/soc-TimeSim.dat", "w")
        for row in csvReader:
            soc = self.getSOC(G=float(row[0]), T=float(row[1]),\
                                P=float(row[5])*float(row[6]))
            f.write(row[4]+", "+str(soc)+'\n')
        f.close()
        
    def readCsvFile(self, fileName):
        f = open(fileName, "r", newline='')
        csvReader = csv.reader(f, delimiter=',')

        # self.saveSOCs(csvReader) ####################
        # f.seek(0)
        # exit()#########################
        
        if self.n == 1:
            for row in csvReader:
                addData = []
                for i in (0, 1, 5, 6):
                    addData.append(float(row[i]))
                if self.socAri == True:
                    soc = self.getSOC(G=float(row[0]), T=float(row[1]),\
                                      P=float(row[5])*float(row[6]))
                    addData.append(soc)
                self.mat = np.append(self.mat, [addData], axis=0)
                self.timeStr = np.append(self.timeStr, row[4])
        else:    
            # read n-1 lines
            preLines = []
            preTimes = []
            preSOCs = []
            count = 0
            for row in csvReader:
                if count == self.n-1:
                    break
                preLines.append(row)
                preTimes.append(row[4])
                preSOCs.append(self.getSOC(G=float(row[0]), T=float(row[1]),\
                                           P=float(row[5])*float(row[6])))
                # print("soc["+str(count)+"]="+str(preSOCs[count]))##################
                count += 1

            #count = 0 #################
            for row in csvReader:
                addData = []

                for k in range(self.n-1):
                    for i in (0, 1, 5, 6):
                        addData.append(float(preLines[k][i]))
                    if self.socAri == True:
                        addData.append(preSOCs[k])
                for i in (0, 1, 5, 6):
                    addData.append(float(row[i]))
                    
                if self.socAri == True:
                    soc = self.getSOC(G=float(row[0]), T=float(row[1]),\
                                               P=float(row[5])*float(row[6]))
                    addData.append(soc)

                    preSOCs.pop(0)
                    preSOCs.append(soc)

                self.mat = np.append(self.mat, [addData], axis=0)

                preLines.pop(0)
                preLines.append(row)

                preTimes.append(row[4])
                self.timeStr = np.append(self.timeStr, preTimes[len(preTimes)-self.n])
                preTimes.pop(0)

    def readCsvFileNoSliding(self, fileName, nSuteru):
        f = open(fileName, "r", newline='')
        csvReader = csv.reader(f, delimiter=',')

        # self.saveSOCs(csvReader) ####################
        # f.seek(0)

        if nSuteru > 0:
            ns = 0
            for row in csvReader:
                ns += 1
                if ns == nSuteru:
                    break
            
        count = 0
        addData = []
        timeStr = ""
        for row in csvReader:
            # for i in (0, 1, 5, 6):######################
            #     addData.append(float(row[i]))

            if self.socAri == True:
                soc = self.getSOC(G=float(row[0]), T=float(row[1]),\
                                           P=float(row[5])*float(row[6]))
                addData.append(soc)

            if count == 0:
                timeStr = row[4]
                
            count += 1

            if count == self.n:
                self.mat = np.append(self.mat, [addData], axis=0)
                addData.clear()
                self.timeStr = np.append(self.timeStr, timeStr)
                count = 0

    def preprocessing(self):
        mean = np.mean(self.mat, axis=0)
        mat0 = self.mat-mean
        # print("mat0.shape=", mat0.shape)
        # print(mat0)
        sum = 0
        for i in range(0, mat0.shape[0]):
            sum += np.inner(mat0[i][:], mat0[i][:])
        sigma = np.sqrt(sum/mat0.shape[0])
        # print("sigma=", sigma)

        matRet = np.ndarray((0, self.matDim1Size))
        timeStrRet = np.ndarray((0,))
        for i in range(0, mat0.shape[0]):
            if np.linalg.norm(mat0[i][:]) <= 3*sigma:
                matRet = np.append(matRet, [mat0[i]/sigma], axis=0)
                timeStrRet = np.append(timeStrRet, self.timeStr[i])

        return matRet, timeStrRet
