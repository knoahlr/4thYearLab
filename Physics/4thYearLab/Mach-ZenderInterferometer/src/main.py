import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import cv2
import matplotlib.pyplot as plt


from pathlib import Path

import os, sys


def fringeFunctionFit(x, k, a, phase):#

    return a*(np.cos((k*x)/2 + phase )*np.cos((k*x)/2 + phase)) #+ offset






def viewPlots(fringeData):     

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
        popt, pcov = curve_fit(fringeFunctionFit, xValues, yValues, p0=[0.05, 45, np.pi/2], maxfev=100000)
        newY = fringeFunctionFit(xValues, *popt)
        print(popt)

        axis.set_yticks(np.arange(0, 300, 25))
        axis.plot(xValues, yValues,'r.' ,label='data')
        axis.plot(xValues, newY, label='Curve_Fit')
        axis.grid(linewidth=0.5)

        # axis.plot(xValues, gaus(xValues,*popt),'r',  label='fit')
        plt.legend()
        
        figNo += 1




if __name__ == "__main__":

    AIR_FRINGE_PATH  = Path(r"D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Mach Zender Experiment\Data\Fringes_NoWater")
    WATER_FRINGE_PATH = None

    logFile = open("../../Logs/MZ.log", 'w')
    sys.stdout = logFile

    fringeData = dict()
    np.set_printoptions(threshold=np.inf)

    
    for fringeFilePath in os.listdir(AIR_FRINGE_PATH):

        print("Noah",fringeFilePath)
        fringeFile = cv2.imread("{0}\{1}".format(AIR_FRINGE_PATH ,fringeFilePath))
        centerFringeMatrix = fringeFile[220:260,:]
        centerAverageArray = centerFringeMatrix.mean(axis=(2,0))

        fringeData[os.path.splitext(fringeFilePath)[0]] = centerAverageArray

        cv2.imwrite("../../Logs/test.jpg", centerFringeMatrix)

    viewPlots(fringeData)
    plt.show()
