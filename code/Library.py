from pylab import *
import numpy as np
import scipy.signal
import cjlib
from sklearn import decomposition
from sklearn.decomposition import PCA
from scipy import ndimage
import cv2
from PIL import Image, ImageDraw

def roipolyny(data, XY):
    img = Image.new('L', (data.shape[-1], data.shape[-2]), 0)
    ImageDraw.Draw(img).polygon(XY, outline=1, fill=1)
    mask = np.array(img)
    return mask


def image_filter(data, number):
	figure_size = 5
	if data.ndim == 3:
		data_filter = zeros((data.shape[0],data.shape[1],data.shape[2]))
		for ii in range(data.shape[0]):
			data_filter[ii] = ndimage.median_filter(data[ii], number)
	elif data.ndim == 4:
		data_filter = zeros((data.shape[0],data.shape[1],data.shape[2],data.shape[3]))
		for ii in range(data.shape[0]):
			for ij in range(data.shape[1]):
				data_filter[ii,ij] = scipy.signal.medfilt2d(data[ii,ij], kernel_size=number)
				# data_filter[ii,ij] = ndimage.median_filter(data[ii,ij], figure_size)
	return data_filter

def applyMask(data, dataMask):
	dataLength = data.shape[0]
	Ints = np.zeros((dataLength))
	for ii in range(dataLength):
	    intsTemp = np.ma.array(data[ii], mask=(dataMask==False))
	    Ints[ii] = intsTemp.mean()
	return Ints


from scipy import optimize
def freqDomainShift(freqdata, Zspectrum):
    
    def lorfitfunc(p, x): 
        tfid = np.ones( np.size(x) ) + p[0]
        for pooli in range( int((len(p)-1)/3) ):
            firstVal = pooli * 3 + 1
            tfid = tfid - p[firstVal] * ( p[firstVal+2]**2 / ( p[firstVal+2]**2 + (p[firstVal+1]-x)**2 ) )
        return tfid
    
    def cestFit( freq, mm, p0):
        errfunc = lambda p, x, y: lorfitfunc(p, x) - y
        p1, success = optimize.leastsq(errfunc, p0[:], args=(freq, mm))
        return p1

    fitInds = np.hstack((np.where( (freqdata >= -0.5) * (freqdata <= 0.5) ),
                        np.where( (freqdata >= freqdata.max()-0.5) )))[0]

    p0 = [0.1, 1., 0.05, 0.5]
    fitParameters = cestFit( freqdata[fitInds], Zspectrum[fitInds], p0)
    fit = lorfitfunc(fitParameters, freqdata)
    return fit, fitParameters

def asymAnalysis(freq, mm):
    negInds, posInds = where(freq<0)[0], where(freq>0)[0]
    Asym =  flipud(mm[negInds]) - mm[posInds]
    return freq[posInds], Asym

def asymMap(freq, cestData, lowFreq, highFreq):
	cestData = squeeze(cestData)
	negInds, posInds = where(freq<0)[0], where(freq>0)[0]
	Asym =  flipud(cestData[negInds,:]) - cestData[posInds,:]
	posFreq = freq[posInds]
	Inds = where((posFreq>=lowFreq) & (posFreq<=highFreq))
	asymCEST = Asym[Inds]
	freqInds = posFreq[Inds]

	return asymCEST, freqInds

def LDMap(freq, cestData, lowFreq, highFreq):
    # select frequency
    inds_sat = np.nonzero( freq < (freq.max()+1) )[0] 

    fitInds = np.hstack((np.where( (freq >= -0.5) * (freq <= 0.5) ),
                        np.where( (freq >= freq.max()-0.5) )))[0]

    lorentz_fitting = zeros((cestData.shape[0], cestData.shape[1], cestData.shape[2]))
    lorDiff = zeros((cestData.shape[0], cestData.shape[1], cestData.shape[2]))

    for ii in range(cestData.shape[1]):
        for ij in range(cestData.shape[2]):
            data = cestData[:,ii,ij]
            lorentz_fitting[:,ii,ij] = freqDomainShift(freq, data)[0]
            lorDiff[:,ii,ij] = lorentz_fitting[:,ii,ij] - data

    Inds = where((freq>=lowFreq) & (freq<=highFreq))
    LD = lorDiff[Inds]
    freqInds = freq[Inds]

    return LD, freqInds