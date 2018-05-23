from scipy import ndimage as snd
import numpy as np
from photutils import centroid_com, centroid_1dg, centroid_2dg 
from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats
import math
import matplotlib.pyplot as plt
from skimage import data 
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from scipy.ndimage import fourier_shift
import numpy as np
import matplotlib.pyplot as plt

def p_theta(Q, U, err_Q, err_U):
    Ip=np.sqrt(U**2+Q**2)
    err_p=np.sqrt((U*err_U)**2+(Q*err_Q)**2)/Ip
    theta=0.5*np.arctan2(U,Q)*180/np.pi
    err_t=0.5*(err_p/Ip)*180/np.pi
    return Ip,err_p,theta,err_t


def rebin(a, fshape):
    if len(a.shape)>1:
        sh = fshape[0],a.shape[0]//fshape[0],fshape[1],a.shape[1]//fshape[1]
        return a.reshape(sh).mean(-1).mean(1)
    else:
        sh = fshape,a.shape[0]//fshape
        return a.reshape(sh).mean(1)
    
def cc_spec(image,width,p_instr,theta_instr): 
   # global make_source_mask
     #crop the images with separate o and e ray sizex is the cropped half xcoordinate sizey1,sizey2 the boundary in y direction
    o_ray=[] ; e_ray=[]; noise=[]; ny=240
    for img in image:
        img=rebin(img[:,20:280],[240,130])
        #print 'size of the image', img.shape
        cen_o=int(round(101.0/(240/ny)))
        cen_e=int(round(135/(240/ny)))
        o_ray.append(img[(cen_o-width):(cen_o+width),:])
        e_ray.append(img[(cen_e-width):(cen_e+width),:])
        mask = make_source_mask(img, snr=1.5, npixels=10, dilate_size=10)
        maskr= np.invert(mask)
        mean, median, std = sigma_clipped_stats(img, sigma=2.3, iters=8, mask=mask)
        noise.append(std)
        
    xc=[];yc=[];crop_image_o=[];crop_image_e=[]
    I_1=np.mean(o_ray[0]+e_ray[0],axis=0); Q_1=np.mean(o_ray[0]-e_ray[0],axis=0)
    I_2=np.mean(o_ray[1]+e_ray[1],axis=0); Q_2=np.mean(o_ray[1]-e_ray[1],axis=0)
    I_3=np.mean(o_ray[2]+e_ray[2],axis=0); U_1=np.mean(o_ray[2]-e_ray[2],axis=0)
    I_4=np.mean(o_ray[3]+e_ray[3],axis=0); U_2=np.mean(o_ray[3]-e_ray[3],axis=0)
    
    Q=(Q_1-Q_2)/2.0
    U=(U_1-U_2)/2.0
    I=(I_1+I_2+I_3+I_4)/4.0
    n=I_1.shape
    ###Instrument correction
    U_c=U-p_instr*np.cos(2.0*theta_instr*np.pi/180.0)*I
    Q_c=Q-p_instr*np.sin(2.0*theta_instr*np.pi/180.0)*I
    
    factor=((2.0*width+1.0)**0.5)
    err_Q=(noise[0]+noise[1])/(2.0*factor)+np.zeros((n[0]))
    #Q_1 noise =sqrt(2)*noise[0]  Q noise=Q_1 noise/sqrt(2)
    err_U=(noise[2]+noise[3])/(2.0*factor)+np.zeros((n[0]))                          
    #same 
    err_I=np.mean(noise)/((2.0)**0.5*factor)+np.zeros((n[0]))
    return Q, U, I, err_Q, err_U, err_I
