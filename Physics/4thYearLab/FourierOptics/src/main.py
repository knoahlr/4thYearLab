import matplotlib
matplotlib.use("Qt5Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys, os, re
from scipy import optimize, signal
from pathlib import Path
import cv2, pandas
import scipy as sp
from scipy.optimize import curve_fit, optimize

FREQ = None

def sincFunc(x, a, b, c, d):

	return a*(np.sinc(b*(x - c)) * np.sinc(b*(x - c)) ) + d
	# return a*np.sinc(b*(x - c)) + d

	# return a*( np.sin(np.pi*b*(x-c)) /( np.pi*(x-c) ) ) + d
def fixedSincFunc(x, a, c, d):

	global FREQ
	freq = FREQ
	return a*(np.sinc(freq*(x - c)) * np.sinc(freq*(x - c)) ) + d

def fixedCircFunc(x, a, c, d):

	global FREQ
	freq = FREQ
	return a*(np.sinc(freq*(x - c)) * np.sinc(freq*(x - c)) ) + d

def besseCircFunc(x, a, b, c, d):

	return a*sp.special.jv(1, b(x-c)) +  d
def gauss2(x, a1, x1, sigma1, a2, x2, sigma2, c):

    return a1*np.exp(-(x-x1)**2/(2*sigma1**2)) + a2*np.exp(-(x-x2)**2/(2*sigma2**2)) + c

def curveFit(xValues, yValues, factor, case="slit", fitType='fixedFreq'):

	''' Curve fitting Slit Transform '''
	global FREQ
	peaks, peakH= sp.signal.find_peaks(yValues, height=20, distance = 20)
	peakHeights = peakH['peak_heights']
	heightGuess = np.amax(peakHeights)
	peakCenterIndex = np.where(peakHeights==heightGuess)[0]



	''' Peak Prominence '''

	peakProminence = signal.peak_prominences(yValues, peaks)[0]
	peakProminenceMean = np.mean(peakProminence)
	peakProminenceStd = np.std(peakProminence)
	newPeaks = peaks
	
	index = 0
	while True:
		if index < len(newPeaks):
			if peakProminence[index] < 10 and peakHeights[index] < 40: #abs((peakProminenceMean - peakProminenceStd))) < 0

				newPeaks = np.delete(newPeaks, index)
				peakProminence = np.delete(peakProminence, index)
				peakHeights = np.delete(peakHeights, index)
			else: 
				index  += 1
		else: break


	if case == 'slit': 
		
		errorPeakValue = np.sort(peakHeights)[::-1][3]
		peakOffset = np.sort(peakHeights)[::-1][0]
		peakOffsetIndex = np.where(peakHeights==peakOffset)[0]
		freqGuess = FREQ#1/abs(xValues[peaks[peakOffsetIndex][0]] - xValues[peaks[peakOffsetIndex+1][0]] )

		centerGuess = xValues[peaks[peakCenterIndex]][0]
		offset = np.amin(yValues)
		# heightGuess = float(abs(yValues[peaks[peakCenterIndex]]-yValues[peaks[peakCenterIndex+1]])**2)
		heightGuess = freqGuess*4#float(abs(yValues[peaks[peakCenterIndex]]-yValues[peaks[peakCenterIndex+1]])**2)

	else: 

		if fitType == 'fixedFreq': 
			errorPeakValue = np.sort(peakHeights)[::-1][1]
		else: 
			errorPeakValue = np.sort(peakHeights)[::-1][3]

		peakOffset = np.sort(peakHeights)[::-1][0]
		peakOffsetIndex = np.where(peakHeights==peakOffset)[0]
		freqGuess = 1/abs(xValues[peaks[peakOffsetIndex]] - xValues[peaks[peakOffsetIndex+1]] )
		centerGuess = xValues[peaks[peakCenterIndex]]
		offset = np.amin(yValues)
		heightGuess = float(abs(yValues[peaks[peakCenterIndex]]-yValues[peaks[peakCenterIndex+1]])**2)

	print('Error Peak Value', errorPeakValue, np.sort(peakHeights)[::-1])


	cFitYValues = list()
	cFitXValues = list()
	corrYValues = list()
	corrXValues = list()

	index = 0
	for value in yValues:
		if not value > errorPeakValue: 

			cFitYValues.append(value)
			cFitXValues.append(index)

			corrYValues.append(value)
			corrXValues.append(index)
			index += 1

		else:

			corrYValues.append(np.nan)
			corrXValues.append(np.nan)
			index += 1

	''' Curve Fit - Sinc function'''



	corrXValues = np.array(corrXValues) * factor
	cFitXValues = np.array(cFitXValues) * factor

	print(heightGuess, freqGuess, centerGuess, offset)
	if fitType == 'fixedFreq':
		if case == 'slit':
			popt, pcov = curve_fit(fixedSincFunc, cFitXValues, cFitYValues, check_finite=False, p0=[ heightGuess, centerGuess, offset], method='lm' , maxfev=100000)
			newy = fixedSincFunc(corrXValues, *popt)
		else:
			popt, pcov = curve_fit(fixedCircFunc, cFitXValues, cFitYValues, check_finite=False, p0=[ heightGuess, centerGuess, offset], method='lm' , maxfev=100000)

			newy = fixedCircFunc(corrXValues, *popt)

	else:
		popt, pcov = curve_fit(sincFunc, cFitXValues, cFitYValues, p0=[ heightGuess, freqGuess, centerGuess, offset] , maxfev=1000000)
		
		print('The observed circular frequency is {0}'.format(popt[1]))

		newy = sincFunc(corrXValues, *popt)

	# print("fit",*popt)

	

	return (corrXValues, corrYValues, newy)





def viewplots(Data, slitWidth, case='slit'):

	''' Multiple fringe pattern images'''
	global FREQ

	pixel_mm = 332.46
	lambdaPixels = 632.8e-09 * pixel_mm *1e3
	focal_pixels = 382 * pixel_mm
	factor = 1/(lambdaPixels * focal_pixels)
	print('factor', factor)

	for key, item in Data.items():
		
		print('\n',key,'\n')
		if case == 'slit': 

			figure = plt.figure()
			axis = figure.gca()
			title = key
			
			yValues = item
			xValues = np.linspace(0, len(yValues)-1, len(yValues))
		

			xValues= xValues * factor
			# print(lambdaPixels, focal_pixels, factor)

			FREQ = slitWidth
			
			corrXValues, corrYValues, newy = curveFit(xValues, yValues, factor, case=case, fitType='fixedFreq')
			corrXValues2, corrYValues2, newy2 = curveFit(xValues, yValues, factor, case=case, fitType='')

			axis.plot(xValues, yValues,'r.', label="Slice of {0}".format(key))
			axis.plot(corrXValues, newy,'b', label="Curve fit with Fixed Frequency of Slice of {0}".format(key))
			axis.plot(corrXValues2, newy2,'g', label="Curve fit  of Slice of {0}".format(key))
			# axis.plot(corrXValues, corrYValues,'c', label="Curve fit of Slice of {0}".format(key))

			axis.grid(b=True) 
			axis.set_title(title)
			axis.set_xlabel('pixels')
			axis.set_ylabel('Amplitude')

		else: 

			FREQ = slitWidth
			
			figure = plt.figure()
			axis = figure.gca()
			yValues = item
			xValues = np.linspace(0, len(yValues)-1, len(yValues))
			xValues= xValues * factor

			corrXValues, corrYValues, newy = curveFit(xValues, yValues, factor, case=case, fitType='fixedFreq')
			corrXValues2, corrYValues2, newy2 = curveFit(xValues, yValues, factor, case=case, fitType='')

			

			axis.plot(xValues, yValues,'r.', label="Slice of {0}".format(key))
			axis.plot(corrXValues, newy,'b', label="Curve fit with Fixed Frequency of Slice of {0}".format(key))
			axis.plot(corrXValues2, newy2,'g', label="Curve fit of Slice of {0}".format(key))
			# axis.plot(corrXValues, corrYValues, 'g', label="Curve fit of Slice of {0}".format(key))
			axis.set_xlabel('pixels')
			axis.set_ylabel('Amplitude')
			axis.grid(b=True)

		plt.legend()

	return  axis, factor, len(yValues)

def conversionFactor(etalonStrip, etalon):

	figure, (axis1, axis2) = plt.subplots(1, 2, sharey=False)
	yValues = etalonStrip
	xValues = np.linspace(0, len(yValues)-1, len(yValues))
	
	if etalon == 'etalon':
		peaks, peakH = sp.signal.find_peaks(yValues, height=0.0055, distance=100)

		peakHeights = peakH['peak_heights']
	
		oneDiv = yValues[0:  int(np.round( (peaks[1] + peaks[2])/2, 0) )]
		x_oneDiv = xValues[0:int(np.round( (peaks[1] + peaks[2])/2, 0) )]

		popt, pcov = curve_fit(gauss2, x_oneDiv, oneDiv, p0=[np.mean(peakHeights), peaks[0], 50, np.mean(peakHeights), peaks[1], 50, np.amin(oneDiv)] , maxfev=100000)

	else: 
		peaks, peakH = sp.signal.find_peaks(yValues, height=0.015)

		peakHeights = peakH['peak_heights']
	
		oneDiv = yValues[0:  int(np.round( (peaks[1] + peaks[2])/2, 0) )]
		x_oneDiv = xValues[0:int(np.round( (peaks[1] + peaks[2])/2, 0) )]

		popt, pcov = curve_fit(gauss2, x_oneDiv, oneDiv, p0=[np.amax(oneDiv), peaks[0], 10, np.amax(oneDiv), peaks[1], 10, np.amin(oneDiv)] , maxfev=100000)

		print("popt before",np.amax(oneDiv), peaks[0], 50, np.amax(oneDiv), peaks[1], 50, np.amin(oneDiv) )
	
	print("popt After", popt)
	newY = gauss2(x_oneDiv, *popt)

	print(" In {0}, 1mm = {1} pixels".format(etalon, popt[4]-popt[1]) )
	axis1.set_title(etalon)
	axis1.grid(b=True)
	axis2.set_title("One Division curve fit - {0}".format(etalon))

	axis1.plot(xValues, yValues,'r' ,label='{0}, oneDiv = 1mm'.format(etalon))
	axis2.plot(x_oneDiv, oneDiv, label='oneDiv')
	axis2.plot(x_oneDiv, newY, label='Gaussian Curve Fit')
	axis1.set_xlabel('distance in pixels')
	axis2.set_xlabel('distance in pixels')
	axis1.set_ylabel('Amplitude')
	axis2.set_ylabel('Amplitude')
	axis2.grid(b=True)

	axis1.legend()
	axis2.legend()
	# print('Etalon', popt[1], popt[4], popt[4]-popt[1], popt[2]/np.sqrt(len(oneDiv)), popt[5]/np.sqrt(len(oneDiv)),  np.sqrt( (popt[2]/np.sqrt(len(oneDiv)) )**2 + (popt[2]/np.sqrt(len(oneDiv)) )**2) )
def apertureSize(apertureProfile, aperture):

	if aperture == 'circ':

		yValues = apertureProfile
		xValues = np.linspace(0, len(yValues)-1, len(yValues))

		peakY = yValues[900:1200]

		peaks, peakH= sp.signal.find_peaks(peakY, height=200, distance = 20)

		peakHeights = peakH['peak_heights']
		sortedPeaks = np.sort(peakHeights)[::-1]
		Width = abs(peaks[np.where(peakHeights==sortedPeaks[0])[0]] - peaks[np.where(peakHeights==sortedPeaks[1])[0]])

		figure = plt.figure()
		axis = figure.gca()

		# print(np.where(yValues == np.sort(yValues)[::-1][-1]), np.sort(yValues)[::-1][-2])
		# circWidth = np.where(yValues==np.sort(yValues)[::-1][-1]) - np.where(yValues==np.sort(yValues)[::-1][-2])
		# circWidth = abs(peaks[1] - peaks[0])
		# threshold = 0.0048 #0.0048
		# Width = 0
		# for value in yValues:
			# if value > threshold: Width += 1
		print('Cicular Aperture Width is {0} pixels'.format(Width))
		axis.plot(xValues, yValues,'r.', label="Profile of {0}".format(aperture))


		axis.grid(b=True) 
		axis.set_title(aperture)
		axis.set_xlabel('pixels')
		axis.set_ylabel('Amplitude')
	else:

		yValues = apertureProfile
		xValues = np.linspace(0, len(yValues)-1, len(yValues))
		threshold = 103

		Width = 0
		for value in yValues:
			if value > threshold: Width += 1


		yGrad = abs(np.gradient(yValues))
		peaks, peakH= sp.signal.find_peaks(yGrad, height=40, distance = 40)
		peakHeights = peakH['peak_heights']
		sortedPeaks = np.sort(peakHeights)[::-1]

		center1Guess = peaks[np.where(peakHeights==sortedPeaks[0])[0]]
		center2Guess = peaks[np.where(peakHeights==sortedPeaks[1])[0]]
		# peaks, peakH= sp.signal.find_peaks(yValues, height=216, distance = 20)

		# print('Slit Width is {0} pixels'.format(abs(center1Guess-center2Guess)))
		figure, (axis1, axis2) = plt.subplots(1, 2, sharey=False)

		# print(np.where(yValues == np.sort(yValues)[::-1][-1]), np.sort(yValues)[::-1][-2])
		# circWidth = np.where(yValues==np.sort(yValues)[::-1][-1]) - np.where(yValues==np.sort(yValues)[::-1][-2])
		# slitWidth = abs(peaks[1] - peaks[0])
		# print('Cicular Aperture Width is {0} pixels'.format(slitWidth))
		print('sit Aperture Width is {0} pixels'.format(Width))
		axis1.plot(xValues, yValues,'r.', label="Slice of {0}".format(aperture))
		axis1.grid(b=True) 
		axis1.set_title('{0} Aperture'.format(aperture))
		axis1.set_xlabel('pixels')
		axis1.set_ylabel('Amplitude')


		axis2.plot(xValues, yGrad,'r.', label=" Gradient of Slice of {0}".format(aperture))
		axis2.grid(b=True) 
		axis2.set_title('Gradient of {0} Aperture'.format(aperture))
		axis2.set_xlabel('pixels')
		axis2.set_ylabel('Amplitude')



def view3Dplot(Data, factor, case='slit'):
	
	for key, item in Data.items():
		
		'''Normalization'''
		# yValues = [value/heightGuess for value in yValues]
		# yValues -= yValues - np.amin(yValues)

		# axis.set_yticks(np.arange(0, 300, 25))
		title = key
		print('\n',key,'\n')
		if case == 'slit': 

			ampl = item
			figure = plt.figure()
			axis = figure.add_subplot(2,1,1,projection='3d')

			''' 3D Plot '''
			print('dim', ampl.shape[0], ampl.shape[1])
			xx, yy = np.mgrid[0:ampl.shape[0], 0:ampl.shape[1]]
			
			xx = xx * factor
			yy = yy * factor
			axis.plot_surface(xx, yy, ampl ,rstride=30, cstride=30)
			axis.set_title(title)
			axis.set_xlabel('pixels')
			axis.set_ylabel('Amplitude')


		else: 

			ampl = item
			figure = plt.figure()
			axis = figure.gca(projection='3d')
			# corrXValues, corrYValues, newy = curveFit(xValues, yValues, case=case)

			''' 3D Plot '''
			print('dim', ampl.shape[0], ampl.shape[1])
			xx, yy = np.mgrid[0:ampl.shape[0], 0:ampl.shape[1]]

			xx = xx * factor
			yy = yy * factor
			
			axis.plot_surface(xx, yy, ampl ,rstride=3, cstride=3)
			axis.set_title(title)
			axis.set_xlabel('pixels')
			axis.set_ylabel('Amplitude')

		# axis.plot(xValues, gaus(xValues,*popt),'r',  label='fit')
		plt.legend()

		mng = plt.get_current_fig_manager()
		mng.window.showMaximized()

	
	return figure

def profile(profileData, case):


	for key, item in profileData.items():


		
		if case == 'slit':
			
			figure, (axis1, axis2) = plt.subplots(1, 2, sharey=False)
			title = key
			yValues = item
			xValues = np.linspace(0, len(yValues)-1, len(yValues))


			yGrad = abs(np.gradient(yValues))
			peaks, peakH= sp.signal.find_peaks(yGrad, height=3, distance = 20)
			peakHeights = peakH['peak_heights']
			sortedPeaks = np.sort(peakHeights)[::-1]

			center1 = peaks[np.where(peakHeights==sortedPeaks[0])[0]]
			center2 = peaks[np.where(peakHeights==sortedPeaks[1])[0]]


			print('sorted Peaks', sortedPeaks)
			print('Peaks', peaks, peakHeights)
			print(' Slit width is {0}'.format(abs(center1-center2)))

			axis1.plot(xValues, yValues,'r.', label="Slice of {0}".format(key))
			axis1.grid(linewidth=0.01) 
			axis1.set_title(title)
			axis1.set_xlabel('pixels')
			axis1.set_ylabel('Amplitude')


			axis2.plot(xValues, yGrad,'r.', label=" Gradient of Slice of {0}".format(key))
			axis2.grid(linewidth=0.01) 
			axis2.set_title(title)
			axis2.set_xlabel('pixels')
			axis2.set_ylabel('Amplitude')
		else:
			# yGrad = abs(np.gradient(yValues))

			# figure, (axis1, axis2) = plt.subplots(1, 2, sharey=False)
			figure = plt.figure()
			axis = figure.gca()

			title = key

			yValues = item
			xValues = np.linspace(0, len(yValues)-1, len(yValues))

			peaks, peakH= sp.signal.find_peaks(yValues, height=3, distance = 20)
			peakHeights = peakH['peak_heights']
			sortedPeaks = np.sort(peakHeights)[::-1]

			heightGuess = np.amax(peakHeights)
			center1Guess = peaks[np.where(peakHeights==sortedPeaks[0])[0]]
			center2Guess = peaks[np.where(peakHeights==sortedPeaks[1])[0]]

			print('center' ,center1Guess, center2Guess)

			popt, pcov = curve_fit(gauss2, xValues, yValues, p0=[heightGuess, center1Guess, 10, heightGuess, center2Guess, 10, 30 ])
			newy = gauss2(xValues, *popt)

			print('Pinhole Width is {0}'.format(popt[4]-popt[1]))
			

			axis.plot(xValues, yValues,'r.', label="Slice of {0}".format(key))
			axis.plot(xValues, newy,'b', label="Curve fit of Slice of {0}".format(key))

			axis.grid(linewidth=0.01) 
			axis.set_title(title)
			axis.set_xlabel('pixels')
			axis.set_ylabel('Amplitude')


			# axis2.plot(xValues, yProfile,'r.', label=" Gradient of Slice of {0}".format(key))
			# axis2.grid(linewidth=0.01) 
			# axis2.set_title(title)
			# axis2.set_xlabel('pixels')
			# axis2.set_ylabel('Amplitude')



def simulation(dimensions, figure, factor, case):
	print('factor', factor)

	if case == 'slit':
		# figure = plt.figure()
		# axis = figure.gca(projection='3d')
		# dimensions = [dimensions[0] * factor, dimensions[1] * factor]
		axis = figure.add_subplot(2,1,2,projection='3d')
		x = np.linspace(0, dimensions[0]* factor, dimensions[0])
		y = np.linspace(0, dimensions[1]* factor, dimensions[1])
		print(x, y)
		xx, yy = np.meshgrid(x, y)
		xCenter = dimensions[0]/2 * factor
		yCenter = dimensions[1]/2 * factor

		print(xCenter, yCenter)
	
		Z = np.sinc(658*(xx - xCenter))**2  @ np.sinc(658*(yy - yCenter))**2  
		# Z = np.sinc((xx**2+yy**2) - dimensions[0])
		# print(Z.shape,len(yy), len(xx))
		# print('3D Data', Z)
		axis.plot_surface(xx, yy, Z ,rstride=30, cstride=30) #, cmap=plt.cm.coolwarm
		axis.set_xlabel('pixels')
		axis.set_ylabel('pixels')
	else:
		pass

if __name__ == "__main__":

	# global FREQ
	# FREQ = None
	# logFile = open("../../Logs/FourierLog.log", 'w')
	# sys.stdout = logFile

	slitTransform = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part1\JosephSlit.png')
	smallSlitTransform = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part2\slit.png')
	circlularTransform = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part2\pinholeCropped.png')
	circlularTransform3D = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part2\pinholeCropped2.png')
	circlularTransform3D  = cv2.resize(circlularTransform3D, (300,300), interpolation=cv2.INTER_CUBIC)

	slitPhone = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\Phone Calibration\slitCalibration.jpg')
	aperturePhone = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\Phone Calibration\circularCalibration.jpg')
	smallSlitPhone = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\Phone Calibration\smallSlitCropped.jpg')
	
	slitProfile = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part2\JosephSquareSlitWFilter.png')
	pinholeProfile = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part2\pinholeHighPass.png')
	etalonStrip = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\calibration\rotatedImage.png')

	slitCenterY = int(slitTransform.shape[0]/2)
	slitCenterX = int(slitTransform.shape[1]/2)

	smallSlitCenterY = int(smallSlitTransform.shape[0]/2)
	smallSlitCenterX = int(smallSlitTransform.shape[1]/2)

	smallSlitTransformXArray = smallSlitTransform[smallSlitCenterY-20:smallSlitCenterY+20,:]
	smallSlitTransformYArray = smallSlitTransform[:,smallSlitCenterX-20:smallSlitCenterX+20]


	slitTransformXArray = slitTransform[slitCenterY-20:slitCenterY+20,:]
	slitTransformYArray = slitTransform[:,slitCenterX-20:slitCenterX+20]

	circlularTransformCenter = int(circlularTransform.shape[0]/2)
	circlularTransformArray = circlularTransform[circlularTransformCenter-50:circlularTransformCenter+50,:]

	slitProfileCenter = int(slitProfile.shape[0]/2)
	slitProfileCenterY = int(slitProfile.shape[1]/2)
	slitProfileArray = slitProfile[slitProfileCenter-20:slitProfileCenter+20,:]
	slitProfileYArray = slitProfile[:,slitProfileCenterY-20:slitProfileCenterY+20]
	
	etalonCenter = int(etalonStrip.shape[0]/2)
	print('etalonCenter',etalonCenter)
	etalonArray = etalonStrip[etalonCenter-100:etalonCenter-80,150:]
	print(etalonArray)
	
	slitPhoneCenter = int(slitPhone.shape[0]/2)
	smallSlitPhoneCenterX = int(smallSlitPhone.shape[0]/2)
	smallSlitPhoneCenterY = int(smallSlitPhone.shape[1]/2)
	aperturePhoneCenter = int(aperturePhone.shape[0]/2)

	circApertureArray = aperturePhone[aperturePhoneCenter-200:aperturePhoneCenter+180,:]
	slitPhoneArray = slitPhone[slitPhoneCenter+10:slitPhoneCenter+50,:]
	smallSlitPhoneArray = smallSlitPhone[smallSlitPhoneCenterX-10:smallSlitPhoneCenterX+10,smallSlitPhoneCenterY-500:smallSlitPhoneCenterY+500]
	
	rulerArray = slitPhone[slitPhoneCenter+300-50:slitPhoneCenter+300,:]

	pinholeProfileCenter = int(pinholeProfile.shape[0]/2)
	pinholeProfileArray =pinholeProfile[pinholeProfileCenter-20:pinholeProfileCenter+20,:]


	slitTransformXArray = slitTransformXArray.mean(axis=(2,0))
	slitTransformYArray = slitTransformYArray.mean(axis=(1,2))   

	smallSlitTransformXArray = smallSlitTransformXArray.mean(axis=(2,0))
	smallSlitTransformYArray = smallSlitTransformYArray.mean(axis=(1,2))  

	circlularTransformArray = circlularTransformArray.mean(axis=(2,0))
	etalonArray =  etalonArray.mean(axis=(2,0))
	rulerArray = rulerArray.mean(axis=(2,0))
	circApertureArray = circApertureArray.mean(axis=(2,0))
	slitPhoneArray = slitPhoneArray.mean(axis=(2,0))
	smallSlitPhoneArray = smallSlitPhoneArray.mean(axis=(2,0))
	

	slitProfileArray = slitProfileArray.mean(axis=(2,0))
	slitProfileYArray = slitProfileYArray.mean(axis=(1,2))
	pinholeProfileArray =pinholeProfileArray.mean(axis=(2,0))

	slitData = {"yTransform":slitTransformYArray, "xTransform":slitTransformXArray}
	smallSlitData = {"smallSlit_xTransform":smallSlitTransformXArray}
	circularData = {"sliceOfCircularTransform":circlularTransformArray}

	circularData3D = {"sliceOfCircularTransform": circlularTransform3D.mean(axis=(2))}
	slitData3D = {"slitTransform":slitTransform.mean(axis=(2))}

	slitProfileData = {"slit Profile X": slitProfileArray, "slit Profile Y": slitProfileYArray}
	pinholeProfileData  = {'pinholeProfile':pinholeProfileArray}

	cv2.imwrite("../../Logs/XArrtest.jpg", circlularTransformArray)
	cv2.imwrite("../../Logs/XArrtest.jpg", slitTransformXArray)
	cv2.imwrite("../../Logs/YArrtest.jpg", slitTransformYArray)
	cv2.imwrite("../../Logs/etalon.jpg", etalonArray)
	cv2.imwrite("../../Logs/circAperture.jpg",circApertureArray)
	cv2.imwrite("../../Logs/slitAperture.jpg",slitPhoneArray)
	cv2.imwrite("../../Logs/ruler.jpg",rulerArray)


	slitAxis, factor, dimension = viewplots(slitData, 714.56)
	slitAxis, factor, dimension = viewplots(smallSlitData, 397.8)
	circAxis, factor, dimension = viewplots(circularData, 864,'circ')

	inverseEtalonStrip = 1/etalonArray
	inverseRuler = 1/rulerArray
	inverseCirc = 1/circApertureArray

	print('Inverse', inverseEtalonStrip)
	print(len(inverseEtalonStrip), len(etalonStrip))
	conversionFactor(inverseEtalonStrip, 'etalon')
	conversionFactor(inverseRuler, 'Ruler')


	# slitAxis = view3Dplot(slitData3D, factor)
	# circAxis = view3Dplot(circularData3D, factor, 'circ')

	# profile(slitProfileData, 'slit')
	# profile(pinholeProfileData, 'circ')
	# apertureSize(inverseCirc, 'circ')
	apertureSize(circApertureArray, 'circ')
	apertureSize(slitPhoneArray, 'slit')
	apertureSize(smallSlitPhoneArray, 'slit')
	# simulation((slitTransform.shape[1],slitTransform.shape[1]), slitAxis, factor, 'slit')


	plt.show()