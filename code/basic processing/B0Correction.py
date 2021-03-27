import numpy as np
from matplotlib import mlab
from scipy import optimize
import nylib, nylib2, cjlib, cxlib

# def wassrProcessing(dataPath, folders, mask):
# 	wassrFreq = nylib.BrukerPar('%s/%s'%(dataPath, folders), 'method', 'MT_Offsets_NoM0=')
# 	wassrData = nylib.Paravision2dseqNew('%s/%s'%(dataPath, folders))                

#     # sort data
# 	inds = np.argsort(wassrFreq)
# 	wassrFreq = wassrFreq[inds]
# 	wassrData = wassrData[:-1][inds]
    
# 	[cestB0Map, delta_freq_map] = getB0Map(wassrFreq, wassrData, mask)
# 	return cestB0Map, delta_freq_map

def wassrProcessing(dataPath, folders, mask):
	# wassrFreq = nylib.BrukerPar('%s/%s'%(dataPath, folders), 'method', 'MT_Offsets_NoM0=')
	wassrFreq = nylib.BrukerPar('%s/%s'%(dataPath, folders), 'method', 'PVM_ppgFreqList1=')
	wassrData = nylib.Paravision2dseqNew('%s/%s'%(dataPath, folders))                
	# print wassrFreq
    # sort data
	allS0 = np.where(wassrFreq >= 20000)[0]
	if len(allS0) == 0:
		tagS0 = -1
		inds = np.argsort(wassrFreq)
		wassrFreq = wassrFreq[inds]
		wassrData = wassrData[:tagS0][inds]
	elif len(allS0) > 0:
		tagS0 = np.min(allS0)
		# print tagS0
		inds = np.argsort(wassrFreq[:tagS0])
		wassrFreq[:tagS0] = wassrFreq[:tagS0][inds]
		wassrData[:tagS0] = wassrData[:tagS0][inds]
    
	[cestB0Map, delta_freq_map] = getB0Map(wassrFreq, wassrData, mask)
	return cestB0Map, delta_freq_map


def getB0Map(freq, data, mask):
	gyr = 42.58 *(10**6)
	B0map = np.zeros((data.shape[1],data.shape[2]))
	delta_freq_map = np.zeros((data.shape[1],data.shape[2]))
	# fit_freq_range = np.where( (freq > -2500) * (freq < 2500) )
	fit_freq_range = np.where( (freq >= -1000) * (freq <= 1000) )
	# print fit_freq_range
	for ii in range(data.shape[1]):
		for ij in range(data.shape[2]):
			if mask[ii, ij]==0:
				continue
			data2fit = 1-data[:, ii, ij] / np.max(data[:, ii, ij])
			Zspectrum = data[:,ii, ij]
			para = lorentzfit(freq[fit_freq_range], data2fit[fit_freq_range], (1e5, 0, 1e5, 0.1))
			delta_freq = para[1]
			delta_freq_map[ii, ij] = delta_freq

			B0map[ii,ij] = -para[1]/gyr
	return B0map, delta_freq_map

def B0correct(dataPath, folders, data, mask, delta_freq_map, freq):
	if freq is None:
		freq = nylib.BrukerPar('%s/%s'%(dataPath, folders[0]), 'method', 'MT_Offsets_NoM0=')
	else:
		freq = freq * 500


	gyr = 42.58 *(10**6)
	if data.ndim == 4:
		correctData = np.zeros((data.shape[0], data.shape[1], data.shape[2], data.shape[3]))
		for di in range(correctData.shape[0]):
			for ii in range(correctData.shape[2]):
				for ij in range(correctData.shape[3]):
					if mask[ii, ij]==0:
						continue
					correctData[di, :, ii, ij] = b0SignalCorr(freq, data[di, :,ii,ij], -1*delta_freq_map[ii, ij])
					# correctData[di, -1, ii, ij] = data[di,-1,ii,ij]
	else:
		correctData = np.zeros((data.shape[0], data.shape[1], data.shape[2]))
		for ii in range(correctData.shape[1]):
			for ij in range(correctData.shape[2]):
				if mask[ii, ij]==0:
					continue
				correctData[:, ii, ij] = b0SignalCorr(freq, data[:,ii,ij], -1*delta_freq_map[ii, ij])
				# correctData[-1, ii, ij] = data[-1,ii,ij]
	return correctData

# def B0correct(dataPath, folders, data, mask, delta_freq_map, freq):
# 	if freq is None:
# 		freq = nylib.BrukerPar('%s/%s'%(dataPath, folders[0]), 'method', 'PVM_ppgFreqList1=')
# 		sort_inds = np.argsort(freq)
# 		freq = freq[sort_inds]
# 	else:
# 		freq = freq * 500

# 	allS0 = np.where(freq >= 20000)[0]
# 	if len(allS0) == 1:
# 		tagS0 = allS0[0]
# 	elif len(allS0) > 1:
# 		tagS0 = np.min(allS0)
	
# 	gyr = 42.58 *(10**6)
# 	if data.ndim == 4:
# 		correctData = np.zeros((data.shape[0], data.shape[1], data.shape[2], data.shape[3]))
# 		for di in range(correctData.shape[0]):
# 			for ii in range(correctData.shape[2]):
# 				for ij in range(correctData.shape[3]):
# 					if mask[ii, ij]==0:
# 						continue
# 					correctData[di, :tagS0, ii, ij] = b0SignalCorr(freq[:tagS0], data[di, :tagS0,ii,ij], -1*delta_freq_map[ii, ij])
# 					correctData[di, tagS0:, ii, ij] = data[di,tagS0:,ii,ij]
# 	else:
# 		correctData = np.zeros((data.shape[0], data.shape[1], data.shape[2]))
# 		for ii in range(correctData.shape[1]):
# 			for ij in range(correctData.shape[2]):
# 				if mask[ii, ij]==0:
# 					continue
# 				correctData[:tagS0, ii, ij] = b0SignalCorr(freq[:tagS0], data[:tagS0,ii,ij], -1*delta_freq_map[ii, ij])
# 				correctData[tagS0:, ii, ij] = data[tagS0,ii,ij]
# 	return correctData

def b0SignalCorr(freq, zspectrum, b0):
	# print freq
	# print b0
	correctedSpectrum = np.zeros((len(freq)))
	for ii in range(correctedSpectrum.shape[0]):
		freqIndex = np.argmin((freq + b0 - freq[ii])**2)
		# freqIndex = find_nearest(freq, (freq[ii]-b0))[1]

		correctedSpectrum[ii] = zspectrum[freqIndex]
	return correctedSpectrum

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx

def lorentzfit(freqdata,Zspectrum,varargin):
	def lfun3c(p0,x):
		F = p0[0] / ((x-p0[1])**2 + p0[2]) + p0[3]
		return F

	def cestFit( freq, mm, p0):
		errfunc = lambda p, x, y: lfun3c(p, x) - y
		p1, success = optimize.leastsq(errfunc, p0[:], args=(freq, mm))
		return p1

	p0 = varargin
	fitParameters = cestFit( freqdata, Zspectrum, p0)
	return fitParameters























