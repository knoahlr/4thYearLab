import matplotlib.pyplot as plt
import numpy as np
import sys, os, re
from scipy import optimize

import cv2

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


if __name__ == "__main__":

    logFile = open("../Logs/ZeemanLog.log", 'w')
    # sys.stdout = logFile

    volume = np.linspace(5,1000, 100000)

    yValues1 = [float(1/pow(elem, 1)) for elem in volume]
    yValues2 = [float(1/pow(elem, 5/3) )for elem in volume]


    plt.plot(volume, yValues1, 'r')
    plt.plot(volume, yValues2, 'b')

    plt.show()


    







