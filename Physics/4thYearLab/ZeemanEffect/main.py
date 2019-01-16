import matplotlib.pyplot as plt
import numpy as np
import sys, os, re
from scipy import optimize, signal
from pathlib import Path
import cv2, pandas

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def viewPlots(fringeData):     

    figNo = 1
    for key, item in fringeData.items():

        figure = plt.figure(figNo)
        axis = figure.gca()

        title = key
        array = item 
        xValues = np.linspace(0, len(array)-1, len(array))

        n = len(xValues)                          #the number of data
        mean = sum(xValues*array)/n                   #note this correction
        sigma = sum(array*(xValues-mean)**2)/n   
        popt,pcov = optimize.curve_fit(gaus,xValues,array,p0=[1,46,6])
        if figNo == 1: 
            peaks = signal.find_peaks(array,height=25)
            print("length:", len(peaks[0]) )
            print(peaks)


        axis.set_yticks(np.arange(0, 300, 25))
        axis.plot(xValues, array, label='data')
        axis.grid(linewidth=0.5)

        axis.plot(xValues, gaus(xValues,*popt),'r',  label='fit')
        plt.legend()
        
        figNo += 1
    
def fourFringes(fringeData, currentRun, writeCSV = True):


    for key, item in fringeData.items():

        if not Path("../Logs/Data/{0}/{1}".format(currentRun, key)).is_dir(): os.makedirs("../Logs/Data/{0}/{1}".format(currentRun, key))
        else: 
            try: os.remove("../Logs/Data/{0}/{1}/fringes.csv".format(currentRun, key))
            except FileNotFoundError: pass
        

        array = item

        peaks = signal.find_peaks(array,height=25)
        lastPeakIndex = peaks[0][16]
        newArray = array[:lastPeakIndex]

        fringeData[key] = newArray
        xValues = np.linspace(0, len(newArray)-1, len(newArray))

        try: 
            fullFringeData = np.column_stack((xValues, newArray, peaks))
            np.save("../Logs/Data/{0}/{1}/fringeData".format(currentRun, key), fullFringeData)
        except ValueError: pass



        np.save("../Logs/Data/{0}/{1}/peaks".format(currentRun, key), peaks[0])
        np.save("../Logs/Data/{0}/{1}/xValues".format(currentRun, key), xValues)
        np.save("../Logs/Data/{0}/{1}/yValues".format(currentRun, key), newArray)


        if writeCSV:

            data = {"pixelNo":xValues, "Amplitude":newArray}
            df = pandas.DataFrame(data=data)
            df.name = key
            with open("../Logs/Data/{0}/{1}/fringes.csv".format(currentRun, key), 'a') as f:
                df.to_csv(f, index=False)





if __name__ == "__main__":

    logFile = open("../Logs/ZeemanLog.log", 'w')
    # sys.stdout = logFile

    currentRun = "Fringes_3"
    fringesFolder = r"D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Zeeman Effect\Data\{0}".format(currentRun)
    fringeData = dict()
    np.set_printoptions(threshold=np.inf)

    for filename in os.listdir(fringesFolder):
        if filename.endswith(('jpg','PNG', 'jpeg')):

            # print(filename)
            fringeFile = cv2.imread("{0}{1}{2}".format(fringesFolder,os.sep, filename))
            centerFringeMatrix = fringeFile[725:735, :]
            centerAverageArray = centerFringeMatrix.mean(axis=(2,0))

            fringeData[os.path.splitext(filename)[0]] = centerAverageArray
            # print(centerAverageArray, centerAverageArray.shape)
            # cv2.imwrite("../Logs/test.jpg", centerFringeMatrix)
            # print(fringeFile.shape)
            # sys.exit()
            # fringeData = fringeFile[]

    fourFringes(fringeData, currentRun)
    viewPlots(fringeData)


    plt.show()









