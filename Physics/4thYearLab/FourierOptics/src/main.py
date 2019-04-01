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


def sincFunc(x, a, b, c, d):

	return a*(np.sinc(b*(x - c)) * np.sinc(b*(x - c)) ) + d
	# return a*np.sinc(b*(x - c)) + d

	# return a*( np.sin(np.pi*b*(x-c)) /( np.pi*(x-c) ) ) + d

def gauss2(x, a1, x1, sigma1, a2, x2, sigma2, c):

    return a1*np.exp(-(x-x1)**2/(2*sigma1**2)) + a2*np.exp(-(x-x2)**2/(2*sigma2**2)) + c

def curveFit(xValues, yValues, factor, case="slit"):

	''' Curve fitting Slit Transform '''

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
		
		errorPeakValue = np.sort(peakHeights)[::-1][1]
		peakOffset = np.sort(peakHeights)[::-1][0]
		peakOffsetIndex = np.where(peakHeights==peakOffset)[0]
		freqGuess = 1/abs(xValues[peaks[peakOffsetIndex]] - xValues[peaks[peakOffsetIndex+1]] )
		centerGuess = xValues[peaks[peakCenterIndex]]
		offset = np.amin(yValues)
		heightGuess = float(abs(yValues[peaks[peakCenterIndex]]-yValues[peaks[peakCenterIndex+1]])**2)

	else: 

		errorPeakValue = np.sort(peakHeights)[::-1][4]
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
	popt, pcov = curve_fit(sincFunc, cFitXValues, cFitYValues, check_finite=False, p0=[ heightGuess, freqGuess, centerGuess, offset], method='lm' , maxfev=100000)
	print("fit",*popt)

	newy = sincFunc(corrXValues, *popt)

	return (corrXValues, corrYValues, newy)





def viewplots(Data, case='slit'):

	''' Multiple fringe pattern images'''

	pixel_mm = 307
	lambdaPixels = 632.8e-09 * pixel_mm *1e3
	focal_pixels = 400 * pixel_mm

	factor = 1/(lambdaPixels * focal_pixels)

	for key, item in Data.items():
		
		'''Normalization'''
		# yValues = [value/heightGuess for value in yValues]
		# yValues -= yValues - np.amin(yValues)

		# axis.set_yticks(np.arange(0, 300, 25))
		print('\n',key,'\n')
		if case == 'slit': 

			figure = plt.figure()
			axis = figure.gca()
			title = key
			
			yValues = item
			xValues = np.linspace(0, len(yValues)-1, len(yValues))
		

			xValues= xValues * factor
			# print(lambdaPixels, focal_pixels, factor)
			
			corrXValues, corrYValues, newy = curveFit(xValues, yValues, factor, case=case)

			axis.plot(xValues, yValues,'r.', label="Slice of {0}".format(key))
			axis.plot(corrXValues, newy,'b', label="Curve fit of Slice of {0}".format(key))
			axis.plot(corrXValues, corrYValues, 'g', label="Curve fit of Slice of {0}".format(key))

			axis.grid(linewidth=0.01) 
			axis.set_title(title)
			axis.set_xlabel('pixels')
			axis.set_ylabel('Amplitude')




		else: 
			
			figure = plt.figure()
			axis = figure.gca()
			yValues = item
			xValues = np.linspace(0, len(yValues)-1, len(yValues))
			xValues= xValues * factor

			corrXValues, corrYValues, newy = curveFit(xValues, yValues, factor, case=case)

			# ''' 3D Plot '''
			# print('dim', ampl.shape[0], ampl.shape[1])
			# xx, yy = np.mgrid[0:ampl.shape[0], 0:ampl.shape[1]]
			
			axis.plot(xValues, yValues,'r.', label="Slice of {0}".format(key))
			axis.plot(corrXValues, newy,'b', label="Curve fit of Slice of {0}".format(key))
			axis.plot(corrXValues, corrYValues, 'g', label="Curve fit of Slice of {0}".format(key))
			axis.set_xlabel('pixels')
			axis.set_ylabel('Amplitude')

		# axis.plot(xValues, gaus(xValues,*popt),'r',  label='fit')
		plt.legend()

		

def view3Dplot(Data, case='slit'):
	
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


def simulation(dimensions, figure, case):

	if case == 'slit':
		# figure = plt.figure()
		# axis = figure.gca(projection='3d')
		axis = figure.add_subplot(2,1,2,projection='3d')
		xx, yy = np.mgrid[0:dimensions[0], 0:dimensions[1]] 
		Z = np.sinc(0.026*(xx - dimensions[0]/2) )**2  @ np.sinc(( 0.026*(yy - dimensions[1]/2) ))**2  
		# Z = np.sinc((xx**2+yy**2) - dimensions[0])
		print(Z.shape,len(yy), len(xx))
		axis.plot_surface(xx, yy, Z ,rstride=30, cstride=30) #, cmap=plt.cm.coolwarm
		axis.set_xlabel('pixels')
		axis.set_ylabel('pixels')
	else:
		pass

if __name__ == "__main__":



	slitTransform = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part1\JosephSlit.png')
	circlularTransform = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part2\pinholeCropped.png')
	circlularTransform3D = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part2\pinholeCropped2.png')
	circlularTransform3D  = cv2.resize(circlularTransform3D, (300,300), interpolation=cv2.INTER_CUBIC)
	
	slitProfile = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part2\JosephSquareSlitWFilter.png')
	pinholeProfile = cv2.imread(r'D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Fourier Optics\Data\FP_March_15_Part2\pinholeHighPass.png')

	slitCenterY = int(slitTransform.shape[0]/2)
	slitCenterX = int(slitTransform.shape[1]/2)
	slitTransformXArray = slitTransform[slitCenterY-20:slitCenterY+20,:]
	slitTransformYArray = slitTransform[:,slitCenterX-20:slitCenterX+20]

	circlularTransformCenter = int(circlularTransform.shape[0]/2)
	circlularTransformArray = circlularTransform[circlularTransformCenter-20:circlularTransformCenter+20,:]

	slitProfileCenter = int(slitProfile.shape[0]/2)
	slitProfileCenterY = int(slitProfile.shape[1]/2)
	slitProfileArray = slitProfile[slitProfileCenter-20:slitProfileCenter+20,:]
	slitProfileYArray = slitProfile[:,slitProfileCenterY-20:slitProfileCenterY+20]

	pinholeProfileCenter = int(pinholeProfile.shape[0]/2)
	pinholeProfileArray =pinholeProfile[pinholeProfileCenter-20:pinholeProfileCenter+20,:]


	slitTransformXArray = slitTransformXArray.mean(axis=(2,0))
	slitTransformYArray = slitTransformYArray.mean(axis=(1,2))   
	circlularTransformArray = circlularTransformArray.mean(axis=(2,0))

	slitProfileArray = slitProfileArray.mean(axis=(2,0))
	slitProfileYArray = slitProfileYArray.mean(axis=(1,2))
	pinholeProfileArray =pinholeProfileArray.mean(axis=(2,0))

	slitData = {"xTransform":slitTransformXArray, "yTransform":slitTransformYArray}
	circularData = {"sliceOfCircularTransform":circlularTransformArray}

	circularData3D = {"sliceOfCircularTransform": circlularTransform3D.mean(axis=(2))}
	slitData3D = {"slitTransform":slitTransform.mean(axis=(2))}

	slitProfileData = {"slit Profile X": slitProfileArray, "slit Profile Y": slitProfileYArray}
	pinholeProfileData  = {'pinholeProfile':pinholeProfileArray}

	cv2.imwrite("../../Logs/XArrtest.jpg", circlularTransformArray)
	cv2.imwrite("../../Logs/XArrtest.jpg", slitTransformXArray)
	cv2.imwrite("../../Logs/YArrtest.jpg", slitTransformYArray)

	viewplots(slitData)
	viewplots(circularData, 'circ')

	# slitAxis = view3Dplot(slitData3D)
	# circAxis = view3Dplot(circularData3D, 'circ')

	# profile(slitProfileData, 'slit')
	# profile(pinholeProfileData, 'circ')

	# simulation((slitTransform.shape[1],slitTransform.shape[1]), slitAxis, 'slit')


	plt.show()