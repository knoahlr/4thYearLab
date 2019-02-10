import matplotlib.pyplot as plt
import numpy as np
import sys, os, re
from scipy import optimize, signal
from pathlib import Path
import cv2, pandas

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def diagonalCut(fringeFile, pixelWidth):


    j= pixelWidth
    newMatrix = np.zeros((pixelWidth, fringeFile.shape[1] - pixelWidth, 3) )

    matrixIndex = 0
    i = 0
    count = 0
    x=0

    while j < (fringeFile.shape[1] - pixelWidth) - 1:
        
        while True:
            try:
                newMatrix[x+count][matrixIndex] = fringeFile[i+ count][j]
                count+=1
                if not count < pixelWidth: break
            except IndexError:
                print(i,j,matrixIndex, count,x)
            
        i+=1
        j+=1
        x=0
        count = 0
        matrixIndex += 1
        if not i< (fringeFile.shape[0] - pixelWidth - 1): break


    return newMatrix


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
        # popt,pcov = optimize.curve_fit(gaus,xValues,array,p0=[1,46,6])
        # if figNo == 1: 
            # peaks = signal.find_peaks(array,height=25)
            # print("length:", len(peaks[0]) )
            # print(peaks)


        axis.set_yticks(np.arange(0, 300, 25))
        axis.plot(xValues, array, label='data')
        axis.grid(linewidth=0.5)

        # axis.plot(xValues, gaus(xValues,*popt),'r',  label='fit')
        plt.legend()
        
        figNo += 1
    
def fourFringes(fringeData, currentRun, writeCSV = True):


    for key, item in fringeData.items():

        if not Path("../Logs/Data/{0}/{1}".format(currentRun, key)).is_dir(): os.makedirs("../Logs/Data/{0}/{1}".format(currentRun, key))
        else: 
            try: os.remove("../Logs/Data/{0}/{1}/fringes.csv".format(currentRun, key))
            except FileNotFoundError: pass
        

        array = item
        # widths = np.full((1, len(array)), 20)
        # peaks = signal.find_peaks_cwt(array, widths)
        # peaks = signal.find_peaks(array, height=20, distance = 30)
        peaks, peakH = signal.find_peaks(array, height=20, distance = 10)

        peakHeights = peakH['peak_heights']
        peakProminence = signal.peak_prominences(array, peaks)[0]
        peakProminenceMean = np.mean(peakProminence)
        peakProminenceStd = np.std(peakProminence)

        newPeaks = peaks
        index = 0
        while True:
            if index < len(newPeaks):
                if peakProminence[index] < 10 and peakHeights[index] < 40: #abs((peakProminenceMean - peakProminenceStd))) < 0
                    print(index)
                    newPeaks = np.delete(newPeaks, index)
                    peakProminence = np.delete(peakProminence, index)
                    peakHeights = np.delete(peakHeights, index)
                else: 
                    index  += 1
            else: break
 
                

        # print(peaks[0])
 




        print("{0}\n\n{1}\n\n{2}\n{3} {4}".format(key, newPeaks, peakProminence, peakProminenceMean, peakProminenceStd))
  


        lastPeakIndex = newPeaks[16]
        newArray = array[:lastPeakIndex]
        fourRingPeaks = newPeaks[:17]
        
        print("\n fourRingPeaks {0} \n\n peaks Again {1}".format(fourRingPeaks, newPeaks))
        fringeData[key] = newArray
        xValues = np.linspace(0, len(newArray)-1, len(newArray))

        try: 
            fullFringeData = np.column_stack((xValues, newArray, peaks))
            np.save("../Logs/Data/{0}/{1}/fringeData".format(currentRun, key), fullFringeData)
        except ValueError: pass



        np.save("../Logs/Data/{0}/{1}/peaks".format(currentRun, key), newPeaks)
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
    sys.stdout = logFile

    currentRun = "Fringes_3"
    fringesFolder = r"D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Zeeman Effect\Data\{0}".format(currentRun)
    fringeData = dict()
    np.set_printoptions(threshold=np.inf)

    for filename in os.listdir(fringesFolder):
        if filename.endswith(('jpg','PNG', 'jpeg')):

            # print(filename)
            fringeFile = cv2.imread("{0}{1}{2}".format(fringesFolder,os.sep, filename))
            centerFringeMatrix = fringeFile[725:735,:]
            # centerFringeMatrix = np.matrix.transpose(fringeFile[:,850:900])
            # centerFringeMatrix = diagonalCut(fringeFile,100)
            centerAverageArray = centerFringeMatrix.mean(axis=(2,0))

            fringeData[os.path.splitext(filename)[0]] = centerAverageArray
            # print(centerAverageArray, centerAverageArray.shape)
            cv2.imwrite("../Logs/test.jpg", centerFringeMatrix)
            # print(fringeFile.shape)
            # sys.exit()
            # fringeData = fringeFile[]

    fourFringes(fringeData, currentRun)
    viewPlots(fringeData)


    plt.show()









