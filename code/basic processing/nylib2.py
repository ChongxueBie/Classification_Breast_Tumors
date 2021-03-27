import cjlib
from PIL import Image, ImageDraw
import numpy as np
from matplotlib.pyplot import *
from scipy import optimize


def roipolyny(data, XY):
    img = Image.new('L', (data.shape[-1], data.shape[-2]), 0)
    ImageDraw.Draw(img).polygon(XY, outline=1, fill=1)
    mask = np.array(img)
    return mask

def applyMask(data, dataMask):
	dataLength = data.shape[0]
	Ints = np.zeros((dataLength))
	for ii in range(dataLength):
	    intsTemp = np.ma.array(data[ii], mask=(dataMask==False))
	    Ints[ii] = intsTemp.mean()
	return Ints

def getRoi(data, slice, inputTime = -1):
	figure(1)
	clf()
	figure(1).subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0)
	# imshow( np.rot90(data[slice]))
	# gray()
	cjlib.mimage( data[slice] )
	
	XY = ginput(inputTime)
	mask = roipolyny(data, XY)
	print (XY)
	
	# plot the roi
	XY=np.array(XY)
	XY = np.concatenate ((XY, [XY[0,:]]), axis=0 )
	plot(XY.transpose()[0], XY.transpose()[1], 'r-*', linewidth = 4)
	figure(1).canvas.draw()
	
	return XY, mask

def sector_mask(shape,centre,radius,angle_range):
    """
    taken from http://stackoverflow.com/questions/18352973/mask-a-circular-sector-in-a-numpy-array
    Return a boolean mask for a circular sector. The start/stop angles in  
    `angle_range` should be given in clockwise order.
    """

    x,y = np.ogrid[:shape[0],:shape[1]]
    cx,cy = centre
    tmin,tmax = np.deg2rad(angle_range)

    # ensure stop angle > start angle
    if tmax < tmin:
            tmax += 2*np.pi

    # convert cartesian --> polar coordinates
    r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy)
    theta = np.arctan2(x-cx,y-cy) - tmin

    # wrap angles between 0 and 2*pi
    theta %= (2*np.pi)

    # circular mask
    circmask = r2 <= radius*radius

    # angular mask
    anglemask = theta <= (tmax-tmin)

    return circmask*anglemask

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
    fitInds = np.hstack((np.where(freqdata<-0.8),np.where(np.abs(freqdata)<1.1),np.where(freqdata>0.8)))[0]
    p0 = [0.1, 1., 0.05, 0.5]
    fitParameters = cestFit( freqdata[fitInds], Zspectrum[fitInds], p0)
    return fitParameters[2]

def printFileContents(filePath):
    f = open(filePath, 'r')
    file_contents = f.read()
    print (file_contents)
    f.close()
