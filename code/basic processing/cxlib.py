from pylab import *
import numpy as np
import scipy.signal
import cjlib
from sklearn import decomposition
from sklearn.decomposition import PCA
from scipy import ndimage
import cv2

def predict_map(n_classes, predict, mask):
	pred_map = zeros((n_classes, mask.shape[0], mask.shape[1]))
	count = 0
	for ii in range(mask.shape[0]):
		for ij in range(mask.shape[1]):
			if mask[ii, ij] == 1:				
				for ni in range(n_classes):
					if predict[count] == ni:
						pred_map[ni, ii, ij] = 1
				count = count + 1
	return pred_map

def save_cest_mask(cestData, mask):
	index = []
	for n in range(cestData.shape[0]):
		for ii in range(mask.shape[0]):
			for ij in range(mask.shape[1]):
				if mask[ii, ij] == True:
					index.append(cestData[n,:,ii,ij])
	index = array(index)
	return index

def asymAnalysis(freq, mm):
    negInds, posInds = where(freq<0)[0], where(freq>0)[0]
    Asym =  flipud(mm[negInds]) - mm[posInds]
    return freq[posInds], Asym

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
#     fitInds = np.hstack((np.where(freqdata<-0.8),np.where(np.abs(freqdata)<1.1),np.where(freqdata>0.8)))[0]
    # fitInds = np.hstack((np.where( (freqdata > 0.01) * (freqdata < 1.5) ),
    #                      np.where( (freqdata > -1.5) * (freqdata < -0.01) ),
    #                     np.where( (freqdata > freqdata.max()-0.6) )))[0]

    fitInds = np.hstack((np.where( (freqdata >= -0.5) * (freqdata <= 0.5) ),
                        np.where( (freqdata >= freqdata.max()-0.5) )))[0]
    # fitInds = np.hstack((np.where(np.abs(freqdata)<lowFreq),np.where(freqdata<(-1*highFreq))))[0]
    # fitInds = np.hstack((np.where(np.abs(freqdata)<lowFreq),np.where(np.abs(freqdata)>highFreq)))[0]
#     fitInds = np.hstack((np.abs(freqdata)>4),(np.abs(freqdata)<1))[0]
    # print freqdata[fitInds]
    p0 = [0.1, 1., 0.05, 0.5]
    fitParameters = cestFit( freqdata[fitInds], Zspectrum[fitInds], p0)
    fit = lorfitfunc(fitParameters, freqdata)
    return fit, fitParameters

def pre_post(cestPre, cestPost, freq, lowFreq, highFreq):
    Inds = where((freq>=lowFreq) & (freq<=highFreq))
    cestPre_Post = []
    indsFreq = []
    for ii in Inds[0]:
        temp = cestPre[ii] - cestPost[ii]
        cestPre_Post.append(temp)
        
        indsFreq.append(freq[ii])
    
    cestPre_Post = array(cestPre_Post)
    cestPre_Post = squeeze(cestPre_Post)
    indsFreq = array(indsFreq)
    indsFreq = squeeze(indsFreq)
    
    return cestPre_Post, indsFreq

def pre_post_Avg(cestPre_Post, indsFreq, beginFreq, endFreq):
    Inds = where((indsFreq>=beginFreq) & (indsFreq<=endFreq))
    cestPre_PostAvg = zeros((cestPre_Post.shape[1],cestPre_Post.shape[2]))
   
    for ii in Inds[0]:
        cestPre_PostAvg = add(cestPre_PostAvg, cestPre_Post[ii])
    cestPre_PostAvg  = cestPre_PostAvg / len(Inds[0])
    return cestPre_PostAvg

def asymMap(freq, cestData, lowFreq, highFreq):
	cestData = squeeze(cestData)
	negInds, posInds = where(freq<0)[0], where(freq>0)[0]
	Asym =  flipud(cestData[negInds,:]) - cestData[posInds,:]
	posFreq = freq[posInds]
	Inds = where((posFreq>=lowFreq) & (posFreq<=highFreq))
	asymCEST = Asym[Inds]
	freqInds = posFreq[Inds]

	return asymCEST, freqInds

def asymMapNOE(freq, cestData, lowFreq, highFreq):
	cestData = squeeze(cestData)
	negInds, posInds = where(freq<0)[0], where(freq>0)[0]
	Asym =  flipud(cestData[posInds,:]) - cestData[negInds,:]
	negFreq = freq[negInds]
	Inds = where((negFreq>=lowFreq) & (negFreq<=highFreq))
	asymCEST = Asym[Inds]
	freqInds = negFreq[Inds]

	return asymCEST, freqInds

def LDMap(freq, cestData, lowFreq, highFreq):
    # select frequency
    inds_sat = np.nonzero( freq < (freq.max()+1) )[0] 
    # fitinds = np.hstack((np.where( (freq > 0.01) * (freq < 1.) ),
    #                         np.where( (freq > -1.) * (freq < -0.01) ),
    #                         np.where( (freq > freq.max()-1.) )))[0]

    # fitinds = np.hstack((np.where( (freq > 0.01) * (freq < 1.5) ),
    #                         np.where( (freq > -1.5) * (freq < -0.01) ),
    #                         np.where( (freq > freq.max()-.4) )))[0]

    fitInds = np.hstack((np.where( (freq >= -0.5) * (freq <= 0.5) ),
                        np.where( (freq >= freq.max()-0.5) )))[0]

    lorentz_fitting = zeros((cestData.shape[0], cestData.shape[1], cestData.shape[2]))
    lorDiff = zeros((cestData.shape[0], cestData.shape[1], cestData.shape[2]))

    for ii in range(cestData.shape[1]):
        for ij in range(cestData.shape[2]):
            data = cestData[:,ii,ij]
            lorentz_fitting[:,ii,ij] = freqDomainShift(freq, data)[0]
            # newfreq, mm_fixed, lorentz_fitting[:,ii,ij], At, x0t, wt, bt, kt = cjlib.cestFit( freq[inds_sat], data[inds_sat], fitinds, freq[inds_sat])
            lorDiff[:,ii,ij] = lorentz_fitting[:,ii,ij] - data

    Inds = where((freq>=lowFreq) & (freq<=highFreq))
    LD = lorDiff[Inds]
    freqInds = freq[Inds]
    # print cestData.shape

    return LD, freqInds

def PCA_denoising(data, pca_method, pect):
	# method1: Malinowskis empirical indicator
	# method2: choose 95%
	# method3: Nelson criterion
	# method4: Median criterion
	def select_k(eigen_value, eigen_vectors, var, pca_method, pect):
	    if pca_method == 1:
	    # method1: select component k
	        ref = 2*np.median(np.sqrt(eigen_value))
	        lam_t = []

	        for ii in range(len(eigen_value)):
	            if np.sqrt(eigen_value[ii]) < ref:
	                lam_t.append(eigen_value[ii])
	        lam_t = array(lam_t)

	        belt = 2.
	        # 1.29
	        lam_selt = []
	        index_selt = []
	        for ii in range(len(eigen_value)):
	            if eigen_value[ii] >= belt*belt * np.median(lam_t):
	                lam_selt.append(eigen_value[ii])
	                index_selt.append(ii)
	        lam_selt = array(lam_selt)
	        index_selt = array(index_selt)

	        k = np.argmax(index_selt)
	    elif pca_method == 2:
	    # method2: select component k
	        # fitinds = (np.where( ( var <= pect )))[0]
	        # k = max(fitinds)
	        k = pect

	    elif pca_method == 3:
	    	n = len(eigen_value)
	    	m = eigen_vectors.shape[0]
	    	RE = np.zeros((n-1,1))
	    	for ki in range(n-1):
	    		data = eigen_value[ki:]
	    		RE[ki] = np.sqrt(np.sum(data) / (m * (n - ki)))

	    	min_RE = np.zeros((n-1,1))
	    	for ki in range(n-1):
	    		min_RE[ki] = RE[ki] / (n-ki)**2

	    	k = np.argmin(min_RE)

	    elif pca_method == 4:
	    	n = len(eigen_value)
	    	r2 = np.zeros((n-1,1))
	    	for ki in range(n-1):
	    		sum1 = 0
	    		for l in range(ki+1,n):
	    			sum1 = sum1 + l*eigen_value[l]
	    		
	    		I = np.array(range(ki+1,n))
	    		sum1 = sum1 * (n-ki)
	    		lam = eigen_value[ki:]
	    		sum_lam = np.sum(lam)
	    		sum_I = np.sum(I)
	    		sum_I2 = np.sum(I**2)
	    		sum2_I = (np.sum(I))**2
	    		sum_lam2 = np.sum(lam**2)
	    		sum2_lam = (np.sum(lam))**2

	    		numerator = sum1 - sum_lam * sum_I
	    		denominator = np.sqrt((n-ki) * sum_I2 - sum2_I) * np.sqrt((n-ki) * sum_lam2 - sum2_lam)
	    		r2[ki] = (numerator / denominator) ** 2

	    	# print r2
	    	fitinds = (np.where( ( r2 >= 0.8 )))[0]
	    	# print fitinds
	    	k = np.min(fitinds)
	    	# print k


	    
	    return k

	covar_matrix = PCA()
	# Calculate Eigenvalues
	covar_matrix.fit(data)
	eigen_vectors = covar_matrix.components_
	eigen_value = covar_matrix.explained_variance_#calculate variance ratios


	var=np.cumsum(np.round(covar_matrix.explained_variance_ratio_, decimals=4)*100)

	k = select_k(eigen_value, eigen_vectors, var, pca_method, pect)
	print (k)

	covar_matrix_new = PCA(n_components = k)
	covar_matrix_new.fit(data)
	data_pca = covar_matrix_new.transform(data)
	data_new = covar_matrix_new.inverse_transform(data_pca)

	return data_new, eigen_value, var

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

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def smooth_data_np_average(arr, span):  # my original, naive approach
    return [np.average(arr[val - span:val + span + 1]) for val in range(len(arr))]

def smooth_data_np_convolve(arr, span):
    return np.convolve(arr, np.ones(span * 2 + 1) / (span * 2 + 1), mode="same")

def smooth_data_np_cumsum_my_average(arr, span):
    cumsum_vec = np.cumsum(arr)
    moving_average = (cumsum_vec[2 * span:] - cumsum_vec[:-2 * span]) / (2 * span)

    # The "my_average" part again. Slightly different to before, because the
    # moving average from cumsum is shorter than the input and needs to be padded
    front, back = [np.average(arr[:span])], []
    for i in range(1, span):
        front.append(np.average(arr[:i + span]))
        back.insert(0, np.average(arr[-i - span:]))
    back.insert(0, np.average(arr[-2 * span:]))
    return np.concatenate((front, moving_average, back))

from scipy.signal import savgol_filter
import pandas as pd

def smooth_Zspec(freq, data, window_size, poly_order):
	# inds = np.hstack((np.where( (freq >= -5.1) * (freq <= -0.5) ),
	# 	np.where( (freq >= 0.5) * (freq <= 5.1) )))[0]

	inds = np.hstack((np.where( (freq >= -5.1) * (freq <= -0.51) ),
		np.where( (freq >= 0.51) * (freq <= 5.1) )))[0]
	inds1 = np.where( (freq >= -5.1) * (freq <= 5.1) )[0]

	print (freq[inds])
	if data.ndim == 3:
		data_filter = zeros(data.shape)
		for mm in range(data.shape[1]):
			for nn in range(data.shape[2]):
				y = data[inds1, mm, nn]
				# data_filter[inds1, mm, nn] = pd.Series(y).rolling(window=window_size).mean()
				data_filter[inds1, mm, nn] = savgol_filter(y, window_size, poly_order, mode='nearest')

		for ii in range(len(freq)):
			if ii in inds:
				continue
			else:
				data_filter[ii] = data[ii]
	elif data.ndim == 4:
		data_filter = zeros(data.shape)
		for ii in range(data.shape[0]):
			for mm in range(data.shape[2]):
				for nn in range(data.shape[3]):
					y = data[ii, inds1, mm, nn]
					# print y.shape
					# data_filter[ii, inds1, mm, nn]  = pd.Series(y).rolling(window=window_size).mean()
					data_filter[ii, inds1, mm, nn] = savgol_filter(y, window_size, poly_order, mode='nearest')

		for jj in range(data.shape[0]):
			for ii in range(len(freq)):
				if ii in inds:
					# print freq[ii]
					continue
				else:
					data_filter[jj,ii] = data[jj,ii]
	return data_filter





