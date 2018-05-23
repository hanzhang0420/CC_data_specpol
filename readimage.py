'''Readimage.py by Han Zhang (2018)

Readimage.py is a function to read CanariCam polarimetry fits.gz files. 
Input: path of the fits file
Output: images list (correspond to images at four different waveplate angles 0, 22.5, 67.5, 45) and instrument position angle

'''
from astropy.io import fits 
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#file_path = '/scratch/hanzhang/YSO_CC/GTC7-14AGCAN/OB0001/object/'
#fits_file = file_path + '0000681098-20140617-CANARICAM-Polarimetry.fits.gz'
def readimage(fits_file):
    image=fits.open(fits_file)
    #image.info()   # gives extensions 
    nsave=image[0].header['NSAVSETS']
    nnods=image[0].header['NNODS']
    nnodsets=image[0].header['NNODSETS']#nod A or B different subtraction order  four different waveplate angles
    ZD1=image[0].header['ZD1'] ; ZD2=image[0].header['ZD2']
    RMA1=image[0].header['ROTATOR1'] ; RMA2=image[0].header['ROTATOR2']
    theta_instr=np.mean([ZD2,ZD2])-np.mean([RMA1,RMA2])
    obj=image[0].header['OBJECT']; f_name=image[0].header['FILTER1'] ; PA=image[0].header['INSTRPA']
    date=image[0].header['DATE']
    print '============================================='
    print 'This is the data of', obj, 'at', date, 'with filter',f_name, 'Position Angle', PA
    print '============================================='
    image_p1=[] ; image_p2=[] ; image_p3=[] ; image_p4=[]  #initial defite four list; list can contain various data types

    for i in range(1,nsave*nnods*nnodsets):
        if image[i].header['NOD'] == 'A':
            fir=0 ; sec=1
        elif image[i].header['NOD'] == 'B':
            fir=1 ; sec=0
        wplate=image[i].header['WPLATE']
        if wplate==float(0.0):
            image_p1.append(image[i].data[fir,:,:]-image[i].data[sec,:,:])
        elif wplate==float(45.0):    
            image_p2.append(image[i].data[fir,:,:]-image[i].data[sec,:,:])
        elif wplate==float(22.5):
            image_p3.append(image[i].data[fir,:,:]-image[i].data[sec,:,:])
        elif wplate==float(67.5):
            image_p4.append(image[i].data[fir,:,:]-image[i].data[sec,:,:])

    print 'for each waveplate, the number of frames', len(image_p1)#image_p1 is a list. len() returns the items in the list 

    image_r=[]
    image_r.append(np.mean(image_p1,axis=0))
    image_r.append(np.mean(image_p2,axis=0))
    image_r.append(np.mean(image_p3,axis=0))
    image_r.append(np.mean(image_p4,axis=0))

    return image_r, theta_instr
