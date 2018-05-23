import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib.ticker import MultipleLocator
from matplotlib import rc, font_manager
from  math import *

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.labelsize'] = 14.
mpl.rcParams['xtick.labelsize'] = 11.
mpl.rcParams['ytick.labelsize'] = 11.

def plot_pol_data(x,y,err_y,xt,yt):
    plt.plot(x,y,'k',label='data')
    plt.errorbar(x,y,yerr=err_y,fmt='k.')
    #plt.xlim([8.0,13.])
    #plt.ylim([0,15])            
    plt.xlabel(xt)
    plt.ylabel(yt)



from cc_spec import p_theta



def plot_pol_decom(q,u,err_q,err_u,q_a,q_e,u_a,u_e,wavet,um):
    fig_size=plt.rcParams["figure.figsize"]

    fig_size[0]=6.0
    fig_size[1]=6.0

    plt.rcParams["figure.figsize"]=fig_size
    p,dp,theta,dt=p_theta(q,u,err_q,err_u)
    p_model=((u_a+u_e)**2.0+(q_a+q_e)**2.0)**(0.5)
    p_a = np.sqrt(q_a**2 + u_a**2) 
    p_e = np.sqrt(q_e**2 + u_e**2) 
    
    theta_a = 0.5*np.arctan(u_a/q_a)*180./np.pi
    theta_e = 0.5*np.arctan(u_e/q_e)*180./np.pi
    
    plt.figure()

    ax=subplot(221)
    plt.plot(wavet,q_a,'r--',label='absorption')
    plt.plot(wavet,q_e,'b:',label='emission')
    plt.plot(wavet,q_a+q_e,'k-',label='P')
    plt.errorbar(wavet,q,yerr=err_q,fmt='k.')
    
    plt.xlim([8.0,13])
    plt.ylim([-um,um])
    plt.setp(ax.get_xticklabels(), visible=False)
    text(11.0,0.8*um,'$q (Q/I) (\%)$',fontsize=13)
    plt.setp(ax.spines.values(), linewidth=2.0)

    
    ax=subplot(222)
    plt.plot(wavet,u_a,'r--',label='absorption')
    plt.plot(wavet,u_e,'b:',label='emission')
    plt.plot(wavet,u_a+u_e,'k-',label='P')
    plt.errorbar(wavet,u,yerr=err_u,fmt='k.')

    plt.xlim([8.0,13])
    plt.ylim([-um,um])
    text(11.0,0.8*um,'$u (U/I) (\%)$',fontsize=13)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_label_position('right')
    plt.setp(ax.spines.values(), linewidth=2.0)

    ax=subplot(223)
    plt.plot(wavet,p_a,'r--',label='absorption')
    plt.plot(wavet,p_e,'b:',label='emission')
    plt.plot(wavet,p_model,'k-',label='P')
    plt.errorbar(wavet,p,yerr=dp,fmt='k.')

    plt.xlabel("Wavelength ($\mu$m)")
    plt.xlim([8.0,13.])
    plt.ylim([0,um])
    text(12.,0.8*um,'$p (\%)$',fontsize=13)
    plt.setp(ax.spines.values(), linewidth=2.0)

    ax=subplot(224)
    plt.plot(wavet,theta_a,'r--',label='Absorption')
    plt.plot(wavet,theta_e,'b:',label='Emission')
      
    plt.errorbar(wavet,theta,yerr=dt,fmt='k.')
    plt.plot(wavet,theta,'k-',label='Total')
    plt.xlabel("Wavelength ($\mu$m)")
    plt.xlim([8.0,13])
    plt.ylim([-100,100])
    text(12.,80,'$\Theta (^{\circ})$',fontsize=13)
    ax.yaxis.set_label_position('right')
    ax.yaxis.set_ticks_position('right')
    legend = plt.legend(loc='upper left',prop={'size': 9})

    subplots_adjust(wspace=0.05,hspace=0.05)
    plt.setp(ax.spines.values(), linewidth=2.0)
    plt.savefig('output/Plot_decom.png',bbox_inches='tight')
    plt.show()
    return p_a,p_e,theta_a,theta_e

    
