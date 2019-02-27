import numpy as np
import pandas as pd
import scipy as sp
from scipy.optimize import curve_fit
from scipy import signal
import cv2
import matplotlib.pyplot as plt


from pathlib import Path

import os, sys

def gauss2(x, a1, x1, sigma1, a2, x2, sigma2, c):

    return a1*np.exp(-(x-x1)**2/(2*sigma1**2)) + a2*np.exp(-(x-x2)**2/(2*sigma2**2)) + c

def tempFringeFunctionFit(x, k, a, c, phase):

    # return a*(1/2 + (1/2)*np.cos(2*np.pi*k*x + phase))
    return a + c*np.cos(2*np.pi*np.exp(-k)*x + phase)

def fringeFunctionFit(x, k, a, x0, sigma, phase):#

    # return a*(np.cos((k*x)/2 + phase )*np.cos((k*x)/2 + phase)) #+ offset

    # return b*np.sin(f*x) + a*(1/2 + (1/2)*np.cos(k*x + phase))

    return a*np.exp(-(x-x0)**2/(sigma**2))*(1/2 + (1/2)*np.cos(k*x + phase))

def findAverageTemp(tempList, range):

    fourFringeTemp = tempList[range[0]:range[1]]
    return np.mean(fourFringeTemp)

def findRatio(etalonStrip, meas):


    figure, (axis1, axis2) = plt.subplots(1, 2, sharey=False)


    yValues = etalonStrip
    xValues = np.linspace(0, len(yValues)-1, len(yValues))

    peaks, peakH = sp.signal.find_peaks(yValues, height=0.04, distance=100)
    print("etalon", peaks, peakH)

    peakHeights = peakH['peak_heights']
    
    oneDiv = yValues[0: int(np.round( (peaks[1] + peaks[2])/2, 0) ) ]
    x_oneDiv = xValues[0: int(np.round( (peaks[1] + peaks[2])/2, 0) )]
    print("Popt Before", [np.mean(peakHeights), peaks[0], 40, np.mean(peakHeights), peaks[1], 40])

    popt, pcov = curve_fit(gauss2, x_oneDiv, oneDiv, p0=[np.mean(peakHeights), peaks[0], 50, np.mean(peakHeights), peaks[1], 50, np.amin(oneDiv)] , maxfev=100000)

    print("popt After", popt)
    newY = gauss2(x_oneDiv, *popt)


    




    if meas == "etalon":
        print("0.1mm = {0} pixels".format(popt[4]-popt[2]) )
        axis1.set_title("Etalon")
        axis2.set_title("One Division curve fit - etalon")
        axis1.plot(xValues, yValues,'r' ,label='etalon, oneDiv = 0.1mm')
        axis2.plot(x_oneDiv, oneDiv, label='oneDiv')
        axis2.plot(x_oneDiv, newY, label='Curve_Fit')

    if meas == "mesh":
        axis1.set_title("Mesh Grid ")
        axis2.set_title("One Division curve fit - mesh")
        axis1.plot(xValues, yValues,'r' ,label='mesh, oneDiv = 0.1mm')
        axis2.plot(x_oneDiv, oneDiv, label='oneDiv')
        axis2.plot(x_oneDiv, newY, label='Curve_Fit')

    axis1.grid(linewidth=0.5)
    axis2.grid(linewidth=0.5)

    axis1.legend()
    axis2.legend()

    return popt[4]-popt[2]

def gratingRatio(gratingData):

    figure, (axis1, axis2) = plt.subplots(1, 2, sharey=False)

    # figure = plt.figure()
    # axis = figure.gca()

    yValues = gratingData
    xValues = np.linspace(0, len(yValues)-1, len(yValues))

    peaks, peakH = sp.signal.find_peaks(yValues, height=0.01, distance=40)
    print("grating", peaks, peakH)

    peakHeights = peakH['peak_heights']
    
    oneDiv = yValues[0: int(np.round( (peaks[1] + peaks[2])/2, 0) ) ]
    x_oneDiv = xValues[0: int(np.round( (peaks[1] + peaks[2])/2, 0) )]

    print("Popt Before", [np.mean(peakHeights), peaks[0], 40, np.mean(peakHeights), peaks[1], 40])


    popt, pcov = curve_fit(gauss2, x_oneDiv, oneDiv, p0=[np.mean(peakHeights), peaks[0], 3.3, np.mean(peakHeights), peaks[1], 3.3, np.amin(oneDiv)] , maxfev=1000)
    
    
    print("popt After", popt)


    newY = gauss2(x_oneDiv, *popt)
    

    # axis.grid(linewidth=0.5)
    # axis.plot(xValues, yValues, 'r')


    axis1.set_title("Mesh on Camera")
    axis2.set_title("One Division curve fit - mesh")
    axis1.plot(xValues, yValues,'r' ,label='Mesh, oneDiv')
    axis2.plot(x_oneDiv, oneDiv, label='oneDiv')
    axis2.plot(x_oneDiv, newY, label='Curve_Fit')

    return popt[4]-popt[2]


def viewPlots(fringeData):     

    ''' Multiple fringe pattern images'''
    figNo = 1
    for key, item in fringeData.items():

        figure = plt.figure(figNo)
        axis = figure.gca()

        title = key
        yValues = item 
        xValues = np.linspace(0, len(yValues)-1, len(yValues))

        n = len(xValues)                          #the number of data
        mean = sum(xValues*yValues)/n                   #note this correction
        sigma = sum(yValues*(xValues-mean)**2)/n  


        ''' Curve fit '''
        peaks, peakH= sp.signal.find_peaks(yValues, height=20, distance = 10)
        peakHeights = peakH['peak_heights']
        heightGuess = np.amax(peakHeights)

        freqGuess = 1/abs(peaks[1]-peaks[2])
        sinFreqGuess = 2/abs(peaks[0]-peaks[-1])
        print("sinFreq", sinFreqGuess)
        maxHeightIndex = np.where(yValues==heightGuess)[0][0]
        # maxHeighIndex = yValues.index(peakHeights)


        sigma= abs(len(xValues) - maxHeightIndex)    
        

        print("peaks and Freq", freqGuess, maxHeightIndex, sigma)
        popt, pcov = curve_fit(fringeFunctionFit, xValues, yValues, p0=[2*np.pi*freqGuess, heightGuess, maxHeightIndex,sigma/3, 120],  maxfev=100000)
        # popt, pcov = curve_fit(fringeFunctionFit, xValues, yValues, p0=[0.1335, 125, 201.3, 67.73, 118.52]) #freqGuess, heightGuess, heightGuess/4, maxHeightIndex,sigma, 120],  maxfev=100000)
        newY = fringeFunctionFit(xValues, *popt)


        print(popt)

        axis.set_yticks(np.arange(0, 300, 25))
        axis.plot(xValues, yValues,'r.' ,label='data')
        axis.plot(xValues, newY, label='Curve_Fit')
        axis.grid(linewidth=0.5)
        axis.set_title(title)

        # axis.plot(xValues, gaus(xValues,*popt),'r',  label='fit')
        plt.legend()
        
        figNo += 1


def fourFringePeaks(fringeWaterData):

    ''' Fringe data collected over time'''

    figNo = 8
    fringes =  fringeWaterData['fringes']
    temp = fringeWaterData['temp']
    fourFringe = dict()

    #  = plt.figure(figNo)
    figure, (axis1, axis2) = plt.subplots(1, 2, sharey=False)
    # axis = figure.gca()

    title = 'Four Fringes - heated Water'
    yValues = fringes#item 
    xValues = np.linspace(0, len(yValues)-1, len(yValues))

    n = len(xValues)                          #the number of data
    mean = sum(xValues*yValues)/n                   #note this correction
    sigma = sum(yValues*(xValues-mean)**2)/n  

    ''' Curve fit '''
    peaks, peakH = sp.signal.find_peaks(yValues, height=0.2, distance = 25)
    peakHeights = peakH['peak_heights']
    heightGuess = np.amax(peakHeights)

    index = 0
    fourFringeCount = 1
    while True:
        if (index + 4) < len(peaks):
            fourFringeTime = peaks[index + 4] - peaks[index]
            fourFringeAvgTemp = findAverageTemp(temp, (peaks[index], peaks[index+4]) )
            fourFringe["fringe{0}".format(fourFringeCount)] = [fourFringeTime, fourFringeAvgTemp]
            fourFringeCount += 1
            index += 4
        else: break


    fringeTimeValues = [item[0] for key, item in fourFringe.items()]
    tempValues = [item[1] for key, item in fourFringe.items()]
    linearX =  np.linspace(0, len(fringeTimeValues)-1, len(fringeTimeValues))

    # 

    # axis.set_yticks(np.arange(0, 300, 25))
    # axis.set_xticks(np.arange(0, 8000, 1))
    # plt.ylim(0, 0.4)s
    # axis.plot(xValues, yValues,'r.' ,label='data')
    axis1.set_ylim(75,500)
    # axis2.set_ylim(0.01, 0.2)
    axis1.plot(linearX, fringeTimeValues,'r.' ,label='Four_fringes')
    axis2.plot(linearX, tempValues,'bo', label='temperature')
    axis1.grid(linewidth=0.5)
    axis2.grid(linewidth=0.5)
    axis1.set_title(title)
    axis2.set_title("Temperature")


    axis1.legend()
    axis2.legend()

    # axis.plot(xValues, gaus(xValues,*popt),'r',  label='fit')
    # plt.legend()

    # print(peaks)
    # freqGuess = 1/abs(peaks[0]-peaks[1])

    

    # print("peaks and Freq four-Fringe", peaks, freqGuess)


    # popt, pcov = curve_fit(tempFringeFunctionFit, xValues, yValues, p0=[-np.log(freqGuess), heightGuess, 0.08, 0],  maxfev=100000)
    # newY = tempFringeFunctionFit(xValues, *popt)

    # print(popt)



if __name__ == "__main__":

    AIR_FRINGE_PATH  = Path(r"D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Mach Zender Experiment\Data\Feb 19\Fringe_NoWater")
    WATER_FRINGE_PATH = Path(r"D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Mach Zender Experiment\Data\Feb 19\Fringes_With_Water")

    logFile = open("../../Logs/MZ.log", 'w')
    sys.stdout = logFile

    fringeData = dict()
    fringeDataWithWater = dict()
    np.set_printoptions(threshold=np.inf)

    
    for fringeFilePath in os.listdir(AIR_FRINGE_PATH):

        print("Noah",fringeFilePath)
        fringeFile = cv2.imread("{0}\{1}".format(AIR_FRINGE_PATH ,fringeFilePath))
        centerFringeMatrix = fringeFile[95:115,:]
        centerAverageArray = centerFringeMatrix.mean(axis=(2,0))

        fringeData[os.path.splitext(fringeFilePath)[0]] = centerAverageArray

        cv2.imwrite("../../Logs/test_{0}.jpg".format(os.path.splitext(fringeFilePath)[0]), centerFringeMatrix)
    
    tempData = pd.read_csv(r"D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Mach Zender Experiment\Data\Feb 19\Fringes_With_Water\fringeTempData.csv")
    etalonData = cv2.imread(r"D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Mach Zender Experiment\Data\Feb_25\etalon_picture.jpg")
    meshData = cv2.imread(r"D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Mach Zender Experiment\Data\Feb_25\mesh_picture.jpg")
    meshDataCamera = cv2.imread(r"D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Mach Zender Experiment\Data\Feb_25\grating.bmp")

    fringeDataWithWater['fringes'] = tempData.values[2300:,0]
    fringeDataWithWater['temp'] = tempData.values[2300:,1]
    
    viewPlots(fringeData)
    fourFringePeaks(fringeDataWithWater)

    etalonStripMatrix = etalonData[1020:1040,:]
    meshStripMatrix = meshData[1020:1040,:]
    meshCameraMatrix = meshDataCamera[20:460,40:600]

    etalonStripAverage = etalonStripMatrix.mean(axis=(2,0))
    meshStripAverage = meshStripMatrix.mean(axis=(2,0))
    meshStripCameraAverage = meshCameraMatrix.mean(axis=(2,0))

    inverseEtalonStrip = 1/etalonStripAverage
    inverseMeshStrip = 1 / meshStripAverage
    inverseMeshCamera = 1/meshStripCameraAverage


    etalonOneDivPixels = findRatio(inverseEtalonStrip, meas="etalon")
    meshOneDivPixels = findRatio(inverseMeshStrip, meas="mesh")
    gratingCameraPixels = gratingRatio(inverseMeshCamera)
    meshOneDivMM = (meshOneDivPixels/etalonOneDivPixels) * 0.1


    print("Mesh --- one division = {0} pixels = {1} mm".format(meshOneDivPixels, meshOneDivMM))
    print("grating --- one division on camera = {0} pixels".format(gratingCameraPixels))

    plt.show()
