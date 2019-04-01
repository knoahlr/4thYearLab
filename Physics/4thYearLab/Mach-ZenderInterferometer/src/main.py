import numpy as np
import pandas as pd
import scipy as sp
from scipy.optimize import curve_fit
from scipy import signal
import cv2
import matplotlib.pyplot as plt
from sympy import sympify, latex, diff, symbols, sqrt, lambdify, exp

from pathlib import Path

import os, sys

def gauss2(x, a1, x1, sigma1, a2, x2, sigma2, c):

    return a1*np.exp(-(x-x1)**2/(2*sigma1**2)) + a2*np.exp(-(x-x2)**2/(2*sigma2**2)) + c

def tempFringeFunctionFit(x, k, a, c, phase):

    # return a*(1/2 + (1/2)*np.cos(2*np.pi*k*x + phase))
    return a + c*np.cos(2*np.pi*np.exp(-k)*x +  phase)

def fringeFunctionFit(x, k, a, x0, sigma, phase):#

    # return a*(np.cos((k*x)/2 + phase )*np.cos((k*x)/2 + phase)) #+ offset

    # return b*np.sin(f*x) + a*(1/2 + (1/2)*np.cos(k*x + phase))

    return a*np.exp(-(x-x0)**2/(sigma**2))*(1/2 + (1/2)*np.cos( 2*np.pi*k*x + phase))

def tempFit(x, a, k, c):

    return a*np.exp(-k*x) + c

    # return sp.special.expit(-k*x)+ c
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
    
    
    if meas == "etalon":

        # oneDiv = yValues[int(np.round( (peaks[0] + peaks[1])/2, 0) ):  int(np.round( (peaks[2] + peaks[3])/2, 0) )]
        # x_oneDiv = xValues[int(np.round( (peaks[0] + peaks[1])/2, 0) ) :int(np.round( (peaks[2] + peaks[3])/2, 0) )]
        oneDiv = yValues[0:  int(np.round( (peaks[1] + peaks[2])/2, 0) )]
        x_oneDiv = xValues[0:int(np.round( (peaks[1] + peaks[2])/2, 0) )]
        print("Popt Before", [np.mean(peakHeights), peaks[0], 40, np.mean(peakHeights), peaks[1], 40])
        popt, pcov = curve_fit(gauss2, x_oneDiv, oneDiv, p0=[np.mean(peakHeights), peaks[0], 50, np.mean(peakHeights), peaks[1], 50, np.amin(oneDiv)] , maxfev=100000)
        print("popt After", popt)
        newY = gauss2(x_oneDiv, *popt)

        print("0.1mm = {0} pixels".format(popt[4]-popt[1]) )
        axis1.set_title("Etalon")
        axis2.set_title("One Division curve fit - etalon")
        
        axis1.plot(xValues, yValues,'r' ,label='etalon, oneDiv = 0.1mm')
        axis2.plot(x_oneDiv, oneDiv, label='oneDiv')
        axis2.plot(x_oneDiv, newY, label='Gaussian Curve Fit')
        axis1.set_xlabel('distance in pixels')
        axis2.set_xlabel('distance in pixels')
        axis1.set_ylabel('Amplitude')
        axis2.set_ylabel('Amplitude')

        print('Etalon', popt[1], popt[4], popt[4]-popt[1], popt[2]/np.sqrt(len(oneDiv)), popt[5]/np.sqrt(len(oneDiv)),  np.sqrt( (popt[2]/np.sqrt(len(oneDiv)) )**2 + (popt[2]/np.sqrt(len(oneDiv)) )**2) )
    else:

        oneDiv = yValues[0: int(np.round( (peaks[1] + peaks[2])/2, 0) )]
        x_oneDiv = xValues[0: int(np.round( (peaks[1] + peaks[2])/2, 0) ) ]
    
        popt, pcov = curve_fit(gauss2, x_oneDiv, oneDiv, p0=[np.mean(peakHeights), peaks[0], 50, np.mean(peakHeights), peaks[1], 50, np.amin(oneDiv)] , maxfev=100000)

        newY = gauss2(x_oneDiv, *popt)

        axis1.set_title("Grating ")
        axis2.set_title("One Division curve fit - Grating")
        axis1.plot(xValues, yValues,'r' ,label='Grating')
        axis2.plot(x_oneDiv, oneDiv, label='oneDiv')
        axis2.plot(x_oneDiv, newY, label='Gaussian Curve Fit')
        axis1.set_xlabel('distance in pixels')
        axis2.set_xlabel('distance in pixels')

        axis1.set_ylabel('Amplitude')
        axis2.set_ylabel('Amplitude')

        print('Grating', popt[1], popt[4], popt[4]-popt[1], popt[2]/np.sqrt(len(oneDiv)), popt[5]/np.sqrt(len(oneDiv)),  np.sqrt( (popt[2]/np.sqrt(len(oneDiv)) )**2 + (popt[2]/np.sqrt(len(oneDiv)) )**2) )


    
    

    axis1.grid(linewidth=0.5)
    axis2.grid(linewidth=0.5)

    axis1.legend()
    axis2.legend()

    return popt[4]-popt[1]

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
    

    axis1.grid(linewidth=0.5)
    axis2.grid(linewidth=0.5)
    # axis.plot(xValues, yValues, 'r')


    # axis1.set_title("Mesh on Camera")
    # axis2.set_title("One Division curve fit - mesh")
    # axis1.plot(xValues, yValues,'r' ,label='Mesh, oneDiv')
    # axis2.plot(x_oneDiv, oneDiv, label='oneDiv')
    # axis2.plot(x_oneDiv, newY, label='Curve_Fit')

    axis1.set_title("Grating on Detector A ")
    axis2.set_title("One Division curve fit - Grating")
    axis1.plot(xValues, yValues,'r' ,label='Grating')
    axis2.plot(x_oneDiv, oneDiv, label='oneDiv')
    axis2.plot(x_oneDiv, newY, label='Gaussian Curve Fit')
    axis1.set_xlabel('distance in pixels')
    axis2.set_xlabel('distance in pixels')

    print('Grating camera', popt[1], popt[4], popt[4]-popt[1], popt[2]/np.sqrt(len(oneDiv)), popt[5]/np.sqrt(len(oneDiv)),  np.sqrt( (popt[2]/np.sqrt(len(oneDiv)) )**2 + (popt[2]/np.sqrt(len(oneDiv)) )**2) )

    return popt[4]-popt[1]


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
        popt, pcov = curve_fit(fringeFunctionFit, xValues, yValues, p0=[freqGuess, heightGuess, maxHeightIndex,sigma/3, 120],  maxfev=100000)
        # popt, pcov = curve_fit(fringeFunctionFit, xValues, yValues, p0=[0.1335, 125, 201.3, 67.73, 118.52]) #freqGuess, heightGuess, heightGuess/4, maxHeightIndex,sigma, 120],  maxfev=100000)
        newY = fringeFunctionFit(xValues, *popt)
        ''' End of Curve fit '''

        print(" K is ", popt[0], np.sqrt(np.diag(pcov))[0])

        axis.set_yticks(np.arange(0, 300, 25))
        axis.plot(xValues, yValues,'r.' ,label='fringe pattern amplitude distribution')
        axis.plot(xValues, newY, label='cosine squared with exponential amplitude')
        axis.grid(linewidth=0.1)
        axis.set_title(title)
        axis.set_xlabel('pixels')
        axis.set_ylabel('Amplitude')

        # axis.plot(xValues, gaus(xValues,*popt),'r',  label='fit')
        plt.legend()
        
        figNo += 1

def linearFit(x, m, b):

    return m*x + b



def fourFringePeaks(fringeWaterData):

    ''' Fringe data collected over time'''

    figNo = 8
    fringes =  fringeWaterData['fringes']
    temp = fringeWaterData['temp']
    fourFringe = dict()

    origFig, (origAxis1, origAxis2) = plt.subplots(2, 1, sharey=False)
     
    figure, (axis1, axis2, axis3) = plt.subplots(1, 3, sharey=False)
    # axis = figure.gca()

    title = 'Four Fringes - heated Water'
    yValues = fringes#item 
    xValues = np.linspace(0, len(yValues)-1, len(yValues))

    n = len(xValues)                          #the number of data
    mean = sum(xValues*yValues)/n                   #note this correction
    sigma = sum(yValues*(xValues-mean)**2)/n  

    ''' Fourier Transformed data'''
    fftY = sp.fft(yValues)
    fftX = np.linspace(0, len(fftY)-1, len(fftY))

    bp=fftY[:]
    
    for i in range(len(bp)): # (H-red)
        if i>=200 and i <= len(bp)-200 :bp[i]=0
    fourierYValues=sp.ifft(bp)

    ''' Curve fit '''
    peaks, peakH = sp.signal.find_peaks(fourierYValues, height=0.2, distance = 25)
    peakHeights = peakH['peak_heights']
    heightGuess = np.amax(peakHeights)
    freqGuess = 1/abs(peaks[0]-peaks[1])
    # print('Freq Guess and {0}, Frequencies allowed,'.format(freqGuess, bp.real))


    popt, pcov = curve_fit(tempFringeFunctionFit, xValues, fourierYValues.real, p0=[-np.log(freqGuess), heightGuess-0.08, 0.08, 45],  maxfev=100000)
    newY = tempFringeFunctionFit(xValues, *popt)


    '''Counting Four Fringe Patterns'''
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

    # yValues = np.cos(2*np.pi*xValues)

    '''Plots '''

    
    # origAxis.set_yticks(np.arange(0, 300, 25))
    # origAxis.set_xticks(np.arange(0, 8000, 1))
    origAxis1.set_ylim(0, 0.4)
    origAxis1.plot(xValues[0:1000], yValues[0:1000],'r.' ,label='raw fringe data')
    origAxis1.plot(fftX[0:1000], fourierYValues[0:1000], 'b', label='Fourier transformed data')
    # origAxis.plot(fftX, newY, 'r')
    origAxis1.set_title('Fringes - first 1000 scans')
    origAxis1.grid()
    origAxis1.legend()

   
    origAxis1.set_ylabel('Intensity')

    origAxis2.set_ylim(0, 0.4)
    origAxis2.plot(xValues[-1000:], yValues[-1000:],'r.' , label='raw fringe data')
    origAxis2.plot(fftX[-1000:], fourierYValues[-1000:], 'b', label='Fourier transformed data')
    # origAxis.plot(fftX, newY, 'r')
    origAxis2.set_title('Fringes - last 1000 scans')
    origAxis2.grid()
    origAxis2.legend()

    origAxis2.set_xlabel('scans, 4 scans/sec')
    origAxis2.set_ylabel('Intensity')

    axis1.set_ylim(75,500)
    # axis2.set_ylim(0.01, 0.2)
    axis1.plot(linearX, fringeTimeValues,'r.' ,label='Four_fringes')
    axis2.plot(linearX, tempValues,'bo', label='temperature')
    axis3.plot(tempValues, fringeTimeValues, 'r.',label='Fringe vs Temp')
    axis1.grid(linewidth=0.5)
    axis2.grid(linewidth=0.5)
    axis1.set_title(title)
    axis2.set_title("Temperature")


    axis1.legend()
    axis2.legend()

    # axis.plot(xValues, gaus(xValues,*popt),'r',  label='fit')
    # plt.legend()

    # print(peaks)
    # 

    

    # print("peaks and Freq four-Fringe", peaks, freqGuess)


    # popt, pcov = curve_fit(tempFringeFunctionFit, xValues, yValues, p0=[-np.log(freqGuess), heightGuess, 0.08, 0],  maxfev=100000)
    # newY = tempFringeFunctionFit(xValues, *popt)

    # print(popt)


def convertVolToTemp(data):

    data = data * (10**6)
    t = symbols('t')
    tempEquation = 0
    C_1 = [0, 2.508355e-02, 7.860106e-08, -2.503131e-10, 8.315270e-14, -1.228034e-17, 9.804036e-22, -4.413030e-26, 1.057734e-30, -1.52755e-35]
    # C_2 = [-1.760041e+1, 3.892120497e+1, 1.8558770032e-2, -9.945759287e-5, 3.184094571e-07, -5.6072844889e-10, 5.6075059059e-13, -3.202072003e-16, 9.7151147152e-20, -1.2104721275e-23]
    # A_0 = 1.185976e+2
    # A_1 = -1.183432e-4

    for index in range(len(C_1)):
        tempEquation += sympify(C_1[index]*(t**index))

    # tempEquation += sympify (A_0*exp(A_1 * (t-126.9686)**2 ) )

    errEquation = diff(tempEquation, t)
    tempEquation = lambdify(t, tempEquation, 'mpmath')
    errEquation = lambdify(t, errEquation, 'mpmath')

    actualTemp = [tempEquation(elem) for elem in data]
    errorOnTemp = [errEquation(elem)*1 for elem in data]

    print('tempConv', actualTemp, "\n\n\n\n\n\n\n\n",errorOnTemp)
    # actualTemp = [tempEquation.evalf(subs={t:elem}) for elem in data]
    # errorOnTemp = [errEquation.evalf(subs={t:elem})*0.00001 for elem in data]

    return actualTemp, errorOnTemp

def fitTemp(fringeData):


    figure = plt.figure()
    axis = figure.gca()
    index = 200

    title = 'Temperature Curve Fit'
    yValues = np.array([float(value) for value in fringeData['temp']], dtype=np.float64)
    # print(type(yValues), yValues)
    xValues = np.linspace(0, len(yValues)-1, len(yValues))

    AmplGuess = np.amax(yValues) - np.amin(yValues)
    offsetGuess = np.amin(yValues)
    coeffGuess = (-1*1/index) * np.log(yValues[index]/np.amax(yValues))
    print("coefficient guess", coeffGuess)

    popt, pcov = curve_fit(tempFit, xValues, yValues, p0=[AmplGuess, coeffGuess, offsetGuess], maxfev=100000)

    print("exponential Coeff, " , popt, np.sqrt(np.diag(pcov)))

    # '''Error on temperature - error propagation on temp fit'''



    # # errorOnTempFit = [errorOnTempFitEquation(fringeData['temp'][index]) * fringeData['errorOntemp'][index] for index in range(len(fringeData['temp']))]

    # '''end of error prop'''

    newY = tempFit(xValues, *popt)

    # print('tempFitData', newY, errorOnTempFit)
    axis.plot(xValues, yValues,'r' ,label='raw temperature readings')
    axis.plot(xValues, newY, label='exponential Curve Fit')
    # axis.errorbar(xValues, newY, yerr=errorOnTempFit, label='Curve_Fit')
    axis.grid(linewidth=0.5)
    axis.set_title(title)
    axis.set_xlabel('scans, 4 scans/sec')
    axis.set_ylabel('Temperature in celsius')

    axis.legend()


    return popt, np.sqrt(np.diag(pcov))

def fourTempPlot(fringeData, coeff, conv):

    figure = plt.figure()
    axis = figure.gca()

    fringes =  fringeData['fringes']
    temp = fringeData['temp']
    tempErrors = fringeData['errorOntemp']
    fourFringe = dict()

    title = 'dn/DT vs Temperature for distilled water'
    yValues = temp#item 
    

    ''' Fourier Transformed data'''
    fftY = sp.fft(fringes)
    fftX = np.linspace(0, len(fftY)-1, len(fftY))

    bp=fftY[:]
    
    for i in range(len(bp)): # (H-red)
        if i>=200 and i <= len(bp)-200 :bp[i]=0
       
    fourierYValues=sp.ifft(bp)


    ''' Finding Peaks'''
    peaks, peakH = sp.signal.find_peaks(fourierYValues, height=0.2, distance = 25)
    peakHeights = peakH['peak_heights']

    '''Counting Four Fringe Patterns'''
    index = 0
    fourFringeCount = 1
    inverseTempDiff = list()
    errorOnDt = list()


    '''Error Propagataion '''

    '''error on temp'''
    a, k, t, c = symbols('a, k, t, c')
    tempFitEquation = sympify(a*exp(-k*t) + c)
    errA = (conv[0])
    errK = (conv[1])
    errC = (conv[2])
    errorOnTempFitEquation =  sqrt((diff(tempFitEquation, a)*errA)**2 + (diff(tempFitEquation, k)*errK)**2 + (diff(tempFitEquation, c)*errC)**2 )
    print(latex(errorOnTempFitEquation))
    errorOnTempFitEquation = lambdify([a, k, t, c], errorOnTempFitEquation, 'mpmath')

    '''error on dn/dT'''

    dt, d, errDt, errD = symbols('dt, d, errDt, errD')
    dndt = sympify(4*632.8e-09/ (dt*d))

    errDnDtEq = sqrt((diff(dndt, dt)*errDt)**2 + (diff(dndt, d)*errD)**2)
    errDnDtEq  = lambdify([dt, d, errDt, errD], errDnDtEq)

    t1Array = list()
    errT1Array = list()
    t2Array = list()
    errT2Array = list()
    dndtArray = list()
    
    finDataFrame =pd.DataFrame()
    
    while True:
        if (index + 4) < len(peaks):

            t1 = tempFit(peaks[index+4], *coeff) 
            t2 =  tempFit(peaks[index], *coeff)
            print(index, t1, t2)
            t1Array.append(t1)
            t2Array.append(t2)

            inverseTempDiff.append(1/(t1 - t2))
            errT1 = errorOnTempFitEquation(coeff[0], coeff[1], t1, coeff[2])
            errT2 = errorOnTempFitEquation(coeff[0], coeff[1], t2, coeff[2])

            errT1Array.append(errT1)
            errT2Array.append(errT2)

            errorOnDt.append(np.sqrt(errT1**2  + errT2**2))
            # fourFringeTime = peaks[index + 4] - peaks[index]
            # fourFringeAvgTemp = findAverageTemp(temp, (peaks[index], peaks[index+4]) )
            # fourFringe["fringe{0}".format(fourFringeCount)] = [fourFringeTime, fourFringeAvgTemp]
            # fourFringeCount += 1
            index += 4
        else: break
    # print(inverseTempDiff)
   


    yValues = [value*(4*623.8e-09)/(49.94e-03) for value in inverseTempDiff]
    dndtArray = yValues
    xValues = np.linspace( tempFit(peaks[0], *coeff), tempFit(peaks[-1], *coeff), len(inverseTempDiff))
    yErr = [errDnDtEq( inverseTempDiff[index], 49.94e-03, errorOnDt[index], 0.005e-03) for index in range(len(inverseTempDiff))]

    finDataFrame['T1'] = t1Array
    finDataFrame['T2'] = t2Array
    finDataFrame['InverseDt'] = [1/elem for elem in inverseTempDiff]
    # finDataFrame['dnDt'] = dndtArray
    finDataFrame['errorT1'] = errT1Array
    finDataFrame['errorT2'] = errT2Array
    finDataFrame['errordt'] = errorOnDt
    # finDataFrame['errordNdT'] = yErr

    finDataFrame.to_csv(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Mach Zender Experiment\finalData.csv')
    # yErr = np.full(len(yValues), 1e-05)

    popt, pcov = curve_fit(linearFit, xValues, yValues, maxfev=10000)

    newY = linearFit(xValues, *popt)


    print('linear Coeff',popt, np.sqrt(np.diag(pcov)))


    ''' end of error prop '''
    # axis.plot(xValues, yValues,'r.' ,label='data')
    axis.errorbar(xValues, yValues, yerr=yErr, fmt='o')
    axis.plot(xValues, newY, 'r', label='curve_fit')
    axis.grid(linewidth=0.5)
    axis.set_title(title)

    axis.set_ylabel('dn/dt (change in refractive index with temperature)')
    axis.set_xlabel('Temperature in  celsius')

    axis.legend()




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
    
    actualTemp, errorOnTemp = convertVolToTemp(fringeDataWithWater['temp'])
    fringeDataWithWater['temp'], fringeDataWithWater['errorOntemp'] = actualTemp, errorOnTemp

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

    coeffTempFit, converg = fitTemp(fringeDataWithWater)
    fourTempPlot(fringeDataWithWater, coeffTempFit, converg)


    print("Mesh --- one division = {0} pixels = {1} mm".format(meshOneDivPixels, meshOneDivMM))
    print("grating --- one division on camera = {0} pixels".format(gratingCameraPixels))

    plt.show()
