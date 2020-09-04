# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 16:21:46 2018

@author: Miguel Castillo
Developed at INL (International Iberian Nanotechnology Laboratory)

>This code is able to implement the transfer matrix method of a 1D photonic 
crystal.
>It finds the Reflectivity and Trasmissivity for any refractive indices for 
wavelength, angle of incidence, or both.
>It can also find the electric field along that 1D photonic crystal but only 
for a single angle and wavelegth.

"""

"""Importing stuff"""
from pylab import *
from numpy import * #numpy
from matplotlib import *
from scipy import interp, arange, exp #to use interpolation
from scipy.interpolate import interp1d #to use interpolation between points
import os #to see and change directory
from collections import OrderedDict
import sys 
# import time #to measure how long it takes to run code
import cmath #to find arcsin of complex numbers
from operator import mul #To get product of all elements of list
from functools import * #so we can find product of list
import warnings #To ignore "ComplexWarning" message
warnings.filterwarnings('ignore') #Continuation of previous line
from random import * #To create random numbers
import matplotlib.pyplot as plt #To be able to contour plot
import matplotlib.mlab as ml #To be able to contour plot
import scipy as sp #To get arcsin of numbers greater than 1
import scipy.constants #To get physical constants

"""
>This file path in the refractive index function might need adjusting
>Always be careful with the complex refractive indices. Here we have used the 
convention of n-i.k
"""
ndirectory = 'C:\\Users\\mcastillo50319\\Desktop\\work stuff\\Code\\refractive_indices' #FIX ME, this might need to be fixed in order to get to the folder with the refractive index data

"""Functions"""

"Refractive Index"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#can make this code shorter by putting all the real part stuff  outside if statments
def refractive_index(onoff, size, units, realpart, imgpart = None):
    """
    >This functions reads off the text files of the real and imaginary part (if
    applicable) of the refractive indices and creates an interpolation between 
    the points. These text files should have the wavelength. 
    >size gives the number of points the user wishes to have
    >units: specify units used for the wavelength since most times it's either 
    in micrometers or nanometers.
    >realpart: insert the name of the file (WITH EXTENSION) with the real part
    of the refractive index.
    >imgpart: insert the name of the file (WITH EXTENSION) with the imaginary
    part of the refractive index, if it exists
    >the last two MUST BE A STRING, i.e. they need to be in "quotes"
    >Write "on" if we want to see the plot of the refractive indices or "off" 
    or "" otherwise. if onoff == float then it's going output a single value 
    of the refractive index at the specified wavelength.
    >Returns an array with the real part of the refractive index, another with
    the imaginary component (if applicable). It can also return an array of 
    wavelengths between min and max wavelenght used in file with as many points
    as size. See example below
    >NOTE THAT WE ARE USING THE CONVENTION n --> n-ik, SO BE CAREFUL WITH HOW
    YOU HAVE YOUR FILES
    
    Examples of how to use it:
        n = refractive_index(820e-9, 500, 1e-6,"GaAs_real.dat", "GaAs_img.dat")
        n = refractive_index(820e-9, 500, 1e-6, "GaAs_real.dat")
        n , wl = refractive_index("", 500, 1e-6, "GaAs_real.dat", "GaAs_img.dat")
        n , wl = refractive_index("on", 500, 1e-6, "GaAs_real.dat", "GaAs_img.dat")
        
    Note that #size# to get a single value of refractive index will only affect
    on how close we are to the specified wavelength so it makes not much 
    difference
        """
    #changing directory to read data set
    current = os.getcwd() # to get current directory
    try:
        os.chdir("/refractive_indices") #if this works, it makes the code slightly faster...
    except:
        os.chdir(ndirectory) 
    if imgpart is None:
        #Real part
        realnfile = open(realpart, "r") #opening file
        column = numpy.loadtxt(realnfile) #extracting arrays
        wlrealn= column[:,0]*units #defining first column as wavelength and making it SI
        nreal = column[:,1] #defining second column as real n
        #interpolating to have a continuos function
        wlnew=linspace(min(wlrealn),(max(wlrealn)),size) #Creates a "continuos" array between the smallest and largest wavelength
        realn=interp(wlnew,wlrealn,nreal) 
        n=realn
        if onoff == "on":
            #plotting real part
            plt.figure()
            plt.plot(wlnew, n.real,label="Real part of refractive index")
            plt.xlabel('${\lambda}/m$')
            plt.title('Real n')
            plt.grid()
            plt.legend(loc = "best")
            plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0,0)) #force exponential in x axis
            os.chdir(current) #to get one folder back in the directory
            return n, wlnew
        elif onoff == "off" or onoff == "":
            os.chdir(current) #to get one folder back in the directory
            return n, wlnew
        elif isinstance(onoff,int) or isinstance(onoff,float):
            if onoff > max(wlnew) or onoff < min(wlnew):
                os.chdir(current) #to get one folder back in the directory
                sys.exit("The wavelength given is out of the range this refractive index works with")
            else:
                pass
            wlwanted = nearest(wlnew, onoff)
            p = [i for i,x in enumerate(wlnew) if x == wlwanted]
            nwanted = n[p][0]
            os.chdir(current) #to get one folder back in the directory
            return nwanted
        else:
            os.chdir(current) #to get one folder back in the directory
            sys.exit("Your first argument of this function should either be 'on', 'off', '', or a number to give a plot, no plot, no plot or a specific refractive index at the wavelength given")
    else:
        #Real part
        realnfile = open(realpart, "r") #opening file
        column = numpy.loadtxt(realnfile) #extracting arrays
        wlrealn= column[:,0]*units #defining first column as wavelength and making it SI
        nreal = column[:,1] #defining second column as real n
        #Imaginary part
        imgnfile = open(imgpart, "r") 
        column = numpy.loadtxt(imgnfile)
        wlimgn = column[:,0]*units 
        nimg = column[:,1]
        #interpolating to have a continuos function
        wlnew = concatenate((wlimgn, wlrealn), axis=None) #joining both arrays into a 1D one as to find the largest and smallest wavelength
        wlnew = linspace(min(wlnew),(max(wlnew)),size) #Creates a "continuous" array between the smallest and largest wavelength
        realn = interp(wlnew,wlrealn,nreal) #interpolation for real part
        imgn = interp(wlnew,wlimgn,nimg) #interpolation for imaginary part
        n = realn - imgn*1j #defining n using n --> n-ik convention
        if onoff == "on":
            #plotting real and imaginary part seperatly
            fig, ax1 = plt.subplots()
            color = 'tab:red'
            ax1.set_xlabel('${\lambda}/m$')
            ax1.set_ylabel('Real n', color=color)
            ax1.plot(wlnew, real(n), color=color)
            ax1.tick_params(axis='y', labelcolor=color)
            
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
            
            color = 'tab:blue'
            ax2.set_ylabel('Imaginary n', color=color)  # we already handled the x-label with ax1
            ax2.plot(wlnew, -imag(n), color=color)
            ax2.tick_params(axis='y', labelcolor=color)
            
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            plt.show()
            plt.title('Real and Imaginary parts of n')
            plt.grid()
            plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0,0)) #force exponential in x axis
            os.chdir(current) #to get one folder back in the directory
            return n, wlnew
        elif onoff == "off" or onoff == "":
            os.chdir(current) #to get one folder back in the directory
            return n, wlnew
        elif isinstance(onoff,int) or isinstance(onoff,float):
            if onoff > max(wlnew) or onoff < min(wlnew):
                os.chdir(current) #to get one folder back in the directory
                sys.exit("The wavelength given is out of the range this refractive index works with")
            else:
                pass
            wlwanted = nearest(wlnew, onoff)
            p = [i for i,x in enumerate(wlnew) if x == wlwanted]
            nwanted = n[p][0]
            os.chdir(current) #to get one folder back in the directory
            return nwanted
        else:
            os.chdir(current) #to get one folder back in the directory
            sys.exit("Your first argument of this function should either be 'on', 'off', '', or a number to give a plot, no plot, no plot or a specific refractive index at the wavelength given")
    
"Overlap between two arrays"""""""""""""""""""""""""""""""""""""""""""""""""""

def overlap (wl1,wl2,size=None): 
    "This function gives the overlap of the range of two arrays, not the overlap of certain points"
    "size: if empty, the function will give only the values in the arrays"
    "If size is not empty, it will be the size of the array between the overlap of the two arrays,"
    "making it semi-continuous. This way is much faster for long arrays" 
    "wl1 and wl2 are simply the two arrays you want to find the overlap"
    if size is None:
        wltest = concatenate((wl1, wl2), axis=None)
        if  min(wl1) < min(wl2):
            wltest = list(filter(lambda x: x >= min(wl2), wltest)) #removes anything smaller than min(wl2)
        else:
            wltest = list(filter(lambda x: x >= min(wl1), wltest)) #similar
        if  max(wl1) < max(wl2):
            wltest = list(filter(lambda x: x <= max(wl1), wltest)) #similar
        else:
            wltest = list(filter(lambda x: x <= max(wl2), wltest)) #similar
        if not wltest: #if list is empty gives error message and exits program
            sys.exit("There is no overlap between the two given arrays")
        else:
            return sort(list(OrderedDict.fromkeys(wltest))) #this removes any duplicate value and sorts the numbers
    else:
        #this makes the code run much faster, but requires that we care only about the initial and final values of the overlap
        wltest = sort([min(wl1),max(wl1),min(wl2),max(wl2)])  
        if  min(wl1) < min(wl2):
            wltest = list(filter(lambda x: x >= min(wl2), wltest)) #removes anything smaller than min(wl2)
        else:
            wltest = list(filter(lambda x: x >= min(wl1), wltest)) #similar
        if  max(wl1) < max(wl2):
            wltest = list(filter(lambda x: x <= max(wl1), wltest)) #similar
        else:
            wltest = list(filter(lambda x: x <= max(wl2), wltest)) #similar
        if not wltest: #if list is empty gives error message and exits program
            sys.exit("There is no overlap between the two given arrays")
        else:
            wltest = list(OrderedDict.fromkeys(wltest)) #this removes any duplicate value
            return linspace(min(wltest),max(wltest),size)
        
"Nearest point of array"""""""""""""""""""""""""""""""""""""""""""""""""""
def nearest(array,value):
    "This function gives the number in the array that's closest to a given value"
    "If function is sorted, it will do a much faster method"
    if all(array==sort(array)) == True:
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
            return array[idx-1]
        else:
            return array[idx]
    else:
        array = np.asarray(array) # in case you do not input an array
        idx = (np.abs(array - value)).argmin() #finds minimum
        return array[idx]

"Check if natural number"""""""""""""""""""""""""""""""""""""""""""""""""""       
def natural(x):
    "this function returns True if x is a natural number, False otherwise"
    return x>0 and isinstance(x,int)

"Put a list to a certain power"""""""""""""""""""""""""""""""""""""""""""""""""""       
def power(list, exponent):
    "this provides and easy wayof having each element of a list into a power"
    return [ x**exponent for x in list ]
  
"Transmittance and Reflectance"""""""""""""""""""""""""""""""""""""""""""""""""""       
def behaviour(size, wl, thetai, md, **options):
    """
    >size gives the number of points (i.e. resolution) you wish to work with
    >wl is the wavelength range you wish to work with, it can either be an 
    array (can have just 2 values: max and min) or a single wavelength
    >thetai is simply either the range or the constant value of the incident 
    angle in degrees. Insert either: const, [initial, final], or [initial, 
    final, number_of_angles]. The number of angles only matters if both wl and
    theta are arrays.
    >if both are arrays, we'll get a contour plot with wavelenght and angle
    versus reflectance and transmittance for each polarisation
    >if wl is an array and thetai is constant, we get a 2D plot of reflectance/
    transmittance vs wavelength and frequency
    >if theta is an array and wl is constant, we get a 2D plot of reflectance/
    transmittance vs thetai
    >if both are constants, we get only one value back (one for each reflection
    and polarisatoin)
    >in md, insert an array with 2 arrays inside, one for refractive indices,
    the other for thicknesses. All of these should be inserted in order. Each 
    one of these can be simplified by multiplying repeated layers. A detailed
    example for the iridoplast is:
        
    size = 300 
    wl = [300e-9,900e-9] #This range can exceed the wl of your refractive indices
    thetai = 30 #incident angle in degrees
    nl = 1.35 #wavelength independent refractive index
    nc = [1e-09, 'Chl_real.dat', 'Chl_img.dat'] #wavelenght dependent and complex refractive index
    M2 = 6.7/2*1e-9 #thicknesses of layers
    M = 6.7*1e-9 #thicknesses of layers
    L = 7.5*1e-9 #thicknesses of layers
    ds = 125*1e-9 #thicknesses of layers
    m = [nl,[nc,[nl,nc]*2,nl,nc,nl]*7,nc,[nl,nc]*2,nl,nc,nl] #structure of the refractive index
    d = [[M2,[L,M]*2,L,M2,ds]*7,M2,[L,M]*2,L,M2] #thicknesses of the layers of iridoplasts
    md = [m,d]
    RsoB, Tso, Rpo, Tpo, wloB, thetaio = behaviour(size, wl, thetai, md, plot = "on", file = "off", frequency = "on", text = "on", cmap = "Blues")
    
    (just run the lines above to see example, you just need to have the refractive index files and have runned all these functions)

    #d, m and ds must always be a list!
    

    >Here, the refractive index is either a constant or a list. In this list,
    we have the units of the wavelength used in the text files, which must be 
    in quotes. These text files give the real and imaginary refractive indices 
    at a certain wavelength. The other values represent thicknesses of layers.
    If the refractive index is real over all wavelenghts, write only the real
    file.
    >Returns the reflectance and transmittance of the s and p waves, the
    wavelength used and the angle used 
    >plot ="on" to show plots
    >file = "on" to enable file saving with data
    >cmap = "name of contouf colour" simply lets the user decide which colour
    to use in the contour. This isnt working!
    >frequency = "on" if we want plots in frequency as well
    >text = "on" to turn on prints
    """
    ###########################################################################
    #First lets deal with the angle
    thetaidummy = thetai #This is just to say which angle is in being used in the plot
    if isinstance(thetai, (list, tuple, np.ndarray)) is False or len(thetai) == 1:
        thetai = linspace(thetai,thetai,size)*pi/180
    elif len(thetai) == 3: #extract number of angles 
        sizetheta = thetai[2]
        thetai = linspace(thetai[0],thetai[1],size)*pi/180
        if natural(sizetheta) is False:
            sys.exit("The number of angles you want to work with (Third element of incident angle) must be a natural number.")
        else:
            pass
    elif len(thetai) > 3:
        sys.exit("The angle must contain only up to three arguments: The first and second tell you the range of angles you want to work in. The third is the number of points in that range you want to work with in the case that both wavelength and incident angle are not constants. If no number is provided, this will be assumed to be 10.")
    elif len(thetai) == 2:
        thetai = linspace(min(thetai),max(thetai),size)*pi/180
        sizetheta = 10 #default size if non is given
    else:
        sys.exit("Something is wrong with the given incident angle.")

    #Check that incident angle is in the right place    
    if max(thetai) > pi/2 or min(thetai) < 0:
        sys.exit("The incident angle must be between 0 and 90 degrees")
    else:
        pass
    
    #extract information from the m arguments
    m = md[0]
    d = md[1]
    #Getting the n separated per layer
    dummy = [1,2,3] #dummy array just to use in while loop
    #making constants in the list into lists to use a the line below it
    while dummy[-1] > dummy[-2]: #stops when length stops increasing i.e. when it's all in a list
        for i in range(len(m)):
            if isinstance(m[i], (list, tuple, np.ndarray)) is False:
                m[i] = [m[i]]               
            else:
                if isinstance(m[i][1], str) is True:
                    m[i] = [m[i]]
                else:
                    pass
        #making list flatter
        m = [item for sublist in m for item in sublist]
        dummy.append(len(m))
    ndummy = m
    #Similarly for thickness, but this time there are no strings so it's simpler
    dummy = [1,2,3] #dummy array just to use in while loop
    #making constants in the list into lists to use a the line below it
    while dummy[-1] > dummy[-2]: #stops when length stops increasing i.e. when it's all in a list
        for i in range(len(d)):
            if isinstance(d[i], (list, tuple, np.ndarray)) is False:
                d[i] = [d[i]]               
            else:
                pass                     
        #making list flatter
        d = [item for sublist in d for item in sublist]
        dummy.append(len(d))
    d = [0] + d + [0] #adding initial and final zero thickness
    
    #Now, we extract the refractive indices and wavelengths we are using
    wldummy = zeros(len(ndummy)).tolist()
    n = zeros(len(ndummy)).tolist()
    for j in range(len(ndummy)):
        if isinstance(ndummy[j], (list, tuple, np.ndarray)) is True:
            if len(ndummy[j]) == 3:
                n[j], wldummy[j] = refractive_index("", size, ndummy[j][0], ndummy[j][1], ndummy[j][2])
            else:
                n[j], wldummy[j] = refractive_index("", size, ndummy[j][0], ndummy[j][1])
        #This will give the refractive index and wavelength. The ndummy[j][n] part gives units, real, and imaginary part respectively
        else:
            n[j], wldummy[j] = ndummy[j], 1 #this gives the constant value for refractive index and an innocuous value of 1 to wl to facilitate the use later 
    ###########################################################################
    #Finding the wavelength range at which we can work
    """
    >In this part, we first check if the dummy array (wldummy) containing all
    the wavelengths has any array in it. Each one of this arrays will be a 
    constraint to the wavelenght we might want to use.
    >Each ith value of wldummy with an array will give the ith value of 
    dummydummy a value of 1, or 0 if it has no array.
    >If dummydummy has no 1's, then there are no constraints and we can just use
    the desired wavelength.
    >If there is only one constraint, then we just find the overlap between the 
    desired wavelength and the constraint
    >If there is more than one, we loop and find the overlap between the desired 
    wavelength and all constraints.
    """
    dummydummy = zeros(len(wldummy)) 
    if isinstance(wl, (list, tuple, np.ndarray)) is False: #if we care only about a constant value, need to make this to avoid error in the overlap process
        wl = [wl,wl]
    else:
        pass
    for i in range(len(wldummy)): #To get an array of 0's and 1's 
        if isinstance(wldummy[i], (list, tuple, np.ndarray)) is True:
            dummydummy[i] = 1
        else:
            dummydummy[i] = 0
    if count_nonzero(dummydummy == 1) == 0:
        wl=linspace(min(wl),max(wl),size)
    elif count_nonzero(dummydummy == 1) == 1:
        dummydummy = dummydummy.tolist() #transform array into list to be able to find position of the value=1
        wl=overlap(wl,wldummy[dummydummy.index(1)],len(wldummy[dummydummy.index(1)])) 
    else:
        indices = [i for i, x in enumerate(dummydummy) if x == 1] #gives an array with the positions of all the elements equal to one in ndummy
        for j in range(len(indices)):
            wl = overlap(wl,wldummy[indices[j]],len(wldummy[indices[j]]))
            
    #Need to change the n's to match in wavelength we can work in
    for i in range(len(n)):   
        if isinstance(n[i], (list, tuple, np.ndarray)) is True:
            n[i] = interp(wl,wldummy[i],n[i]) #interpolation that cuts down the n's to the wavelength we're using
        else:
            pass

    #Error message if user has put the wrong m, n combination
    if len(m) != len(d):
        sys.exit("Your refractive index vector must have the same length as your as your thicknesses vector plus two")
    else:
        pass

    """at this point we have all the refractive indices and thicknesses of the 
    layers in form of several arrays, where each element of that array 
    corresponds to one layer"""
    
    #now we check which situation we're in i.e wl or thetat, are they constant?
    if (thetai[0] == thetai[1] and wl[0] != wl[1]) or (thetai[0] != thetai[1] and wl[0] == wl[1]): #if one is a constant and the other is an array
        #Defining variebles for each layer. Each layer needs to have a number of points with the length of "size"
        theta = zeros(len(n)).tolist() #ith element gives the angle at the ith layer
        theta[0] = thetai 
        phi = zeros(len(n)).tolist() #This will have the phase of each intermediate layer
        #Making constant refractive indices into same length arrays and finding angles and phases at different layers
        for i in range(len(n)):
            n[i] = n[i]*ones(size) #To make sure constant n values don't bug 
            theta[i] = array(list(map(cmath.asin,(n[0]*sin(thetai)/n[i]))))
            phi[i] = 2*pi*n[i]*d[i]*cos(theta[i])/wl #since thickness is for initial and final is zero, their phi is also zero
            #thetai and wl make sure each layer has "size" points
        #Create dummy arrays to have dynamical and propagation matrices
        Ds = zeros(len(n)).tolist()
        Dp = zeros(len(n)).tolist()
        P = zeros(len(n)).tolist()
        Kls = zeros(len(n)-2).tolist() #This serves as a dummy array to get D*P*D^-1 for each layer
        Klp = zeros(len(n)-2).tolist() 
       #Create dummy arrays to have have the reflectance/transmittance 
        Rs = zeros(size)
        Ts = zeros(size)
        Rp = zeros(size)
        Tp = zeros(size)
        #Actual maths
        for i in range(size): #loops over all points 
            for j in range(len(n)): #loops over all media
                #Dynamical Matrices for s/TE waves
                Ds[j] = matrix([[1,1],[n[j][i]*cos(theta[j][i]),-n[j][i]*cos(theta[j][i])]])
                #Dynamical Matrices for p/TM waves  
                Dp[j] = matrix([[cos(theta[j][i]),cos(theta[j][i])],[n[j][i],-n[j][i]]])
                #Propagation Matrix
                P[j] = matrix([[exp(1j*phi[j][i]),0],[0,exp(-1j*phi[j][i])]]) 
            for j in range(len(n)-2): #need to exclude initial and final 
                #dummy intermediate arrays
                Kls[j] = Ds[j+1]*P[j+1]*Ds[j+1].I
                Klp[j] = Dp[j+1]*P[j+1]*Dp[j+1].I
            Ks = reduce(mul,Kls)
            Kp = reduce(mul,Klp)
            #transfer matrix
            Ms = Ds[0].I*Ks*Ds[-1]
            Mp = Dp[0].I*Kp*Dp[-1]
            #Reflectance and transmittance for s/TE waves
            Rs[i] = (abs(Ms[1,0]/Ms[0,0]))**(2)
            Ts[i] = real(n[-1][i]*cos(theta[-1][i]))/real(n[0][i]*cos(theta[0][i]))*(abs(1/Ms[0,0]))**2
            #Reflectance and transmittance for p/TM waves
            Rp[i] = (abs(Mp[1,0]/Mp[0,0]))**(2)
            Tp[i] = real(n[-1][i]*cos(np.conjugate(theta[-1][i])))/real(n[0][i]*cos(np.conjugate(theta[0][i])))*(abs(1/Mp[0,0]))**2

        if thetai[0] == thetai[1] and wl[0] != wl[1]:
            if options.get("plot") == "on":
                "Wavelength Plots"
                #s/TE waves plots 
                plt.figure()
                if thetaidummy == 0: #simplify stuff if s and p waves are the same
                    plt.title('s/TE and p/TM waves')
                else:
                    plt.subplot(2, 1, 1)
                    plt.title('s/TE waves')
                plt.plot(wl,Rs,label = 'R')
                plt.plot(wl,Ts,label="T")
                plt.plot(wl,Ts+Rs,label="R+T")
                plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #force exponential in x axis
                plt.grid()
                plt.legend(loc="best")
                plt.xlabel('${\lambda}/m$') 
                #p/TM waves plots
                if thetaidummy == 0:
                    pass
                else:
                    plt.subplot(2, 1, 2)
                    plt.plot(wl,Rp,label = 'R')
                    plt.plot(wl,Tp,label="T")
                    plt.plot(wl,Tp+Rp,label="R+T")
                    plt.xlabel('${\lambda}/m$') 
                    plt.title('p/TM waves')
                    plt.grid()  
                    plt.legend(loc="best")
                plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #force exponential in x axis
                plt.figtext(0.5, 0, "${\Theta}_i$ = %.2f$^\circ$" %(thetaidummy), position=(0.137,0.925), fontsize = 12, bbox={"facecolor":"orange", "alpha":0.5, "pad":5}) 
                                
                "Go from wavelength to frequency"
                if options.get("frequency") == "on": 
                    from scipy.constants import c #get proper speed of light
                    omega = 2*pi*c/wl
                    #s/TE waves plots 
                    plt.figure()
                    if thetaidummy == 0: #simplify stuff if s and p waves are the same
                        plt.title('s/TE and p/TM waves')
                    else:                        
                        plt.subplot(2, 1, 1)
                        plt.title('s/TE waves')
                    plt.plot(omega, Rs, label = 'R')
                    plt.plot(omega ,Ts, label = "T")
                    plt.plot(omega, Ts + Rs, label = "R+T")
                    plt.grid()
                    plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0,0)) #force exponential in x axis
                    plt.legend(loc = "best")
                    
                    #p/TM waves plots 
                    if thetaidummy == 0:
                        pass
                    else:
                        plt.subplot(2, 1, 2)
                        plt.plot(omega ,Rp ,label = 'R')
                        plt.plot(omega, Tp, label = "T")
                        plt.plot(omega, Tp + Rp, label = "R+T")
                        plt.title('p/TM waves')
                        plt.grid()
                        plt.legend(loc = "best")
                    plt.xlabel('${\omega}/s$')
                    plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0,0)) #force exponential in x axis
                    plt.figtext(0.5, 0, "${\Theta}_i$ = %.2f$^\circ$" %(thetaidummy), position = (0.137,0.925), fontsize = 12, bbox = {"facecolor":"orange", "alpha":0.5, "pad":5})
                    
                else:
                    pass
            else:
                pass
            "Saving data in file"
            if options.get("file") == "on":
                #Rs
                Rsdummy = list(zeros(len(wl)))
                for i in range(len(wl)):
                    Rsdummy[i] = [0,0]
                    Rsdummy[i][0] = Rs[i]
                    Rsdummy[i][1] = wl[i]
                np.savetxt('Rswl.txt', Rsdummy, header= "the incident angle is %f degrees. First column is Reflectance of s waves, second is wavelength" %thetaidummy )
                #Rp
                Rpdummy = list(zeros(len(wl)))
                for i in range(len(wl)):
                    Rpdummy[i] = [0,0]
                    Rpdummy[i][0] = Rp[i]
                    Rpdummy[i][1] = wl[i]
                np.savetxt('Rpwl.txt', Rpdummy, header= "the incident angle is %f degrees. First column is Reflectance of p waves, second is wavelength" %thetaidummy)
                #Ts
                Tsdummy = list(zeros(len(wl)))
                for i in range(len(wl)):
                    Tsdummy[i] = [0,0]
                    Tsdummy[i][0] = Ts[i]
                    Tsdummy[i][1] = wl[i]
                np.savetxt('Tswl.txt', Tsdummy, header= "the incident angle is %f degrees. First column is Transmittance of s waves, second is wavelength" %thetaidummy)
                #Tp
                Tpdummy = list(zeros(len(wl)))
                for i in range(len(wl)):
                    Tpdummy[i] = [0,0]
                    Tpdummy[i][0] = Tp[i]
                    Tpdummy[i][1] = wl[i]
                np.savetxt('Tpwl.txt', Tpdummy, header= "the incident angle is %f degrees. First column is Transmittance of p waves, second is wavelength" %thetaidummy)
            else:
                pass
        elif thetai[0] != thetai[1] and wl[0] == wl[1]:
            if options.get("plot") == "on":
                "Angle Plots"
                #s/TE waves plots 
                plt.figure()
                plt.subplot(2, 1, 1)
                plt.plot(thetai, Rs, label = 'R')
                plt.plot(thetai, Ts, label="T")
                plt.plot(thetai, Ts + Rs, label="R+T")
                plt.title('s/TE waves')
                plt.grid()
                plt.legend(loc="best")
                
                #p/TM waves plots
                plt.subplot(2, 1, 2)
                plt.plot(thetai, Rp, label = 'R')
                plt.plot(thetai, Tp, label="T")
                plt.plot(thetai, Tp + Rp, label="R+T")
                plt.xlabel('${\Theta}/rad$') 
                plt.title('p/TM waves')
                plt.grid()
                plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #force exponential in x axis
                plt.figtext(0.5, 0, "${\lambda}$ = %.2e m" %wl[0], position = (0.137,0.925), fontsize = 12, bbox={"facecolor":"orange", "alpha":0.5, "pad":5}) 
                plt.legend(loc = "best" )
            else:
                pass
            "Saving data in file"
            if options.get("file") == "on":
                #Rs
                Rsdummy = list(zeros(len(wl)))
                for i in range(len(wl)):
                    Rsdummy[i] = [0,0]
                    Rsdummy[i][0] = Rs[i]
                    Rsdummy[i][1] = thetai[i]
                np.savetxt('Rst.txt', Rsdummy, header= "the wavelength used is %.15f meters. First column is Reflectance of s waves, second is angle" %wl[0])
                #Rp
                Rpdummy = list(zeros(len(wl)))
                for i in range(len(wl)):
                    Rpdummy[i] = [0,0]
                    Rpdummy[i][0] = Rp[i]
                    Rpdummy[i][1] = thetai[i]
                np.savetxt('Rpt.txt', Rpdummy, header= "the wavelength used is %.5e meters. First column is Reflectance of p waves, second is angle" %wl[0])
                #Ts
                Tsdummy = list(zeros(len(wl)))
                for i in range(len(wl)):
                    Tsdummy[i] = [0,0]
                    Tsdummy[i][0] = Ts[i]
                    Tsdummy[i][1] = thetai[i]
                np.savetxt('Tst.txt', Tsdummy, header= "the wavelength used is %.5e meters. First column is Transmittance of s waves, second is angle" %wl[0])
                #Tp
                Tpdummy = list(zeros(len(wl)))
                for i in range(len(wl)):
                    Tpdummy[i] = [0,0]
                    Tpdummy[i][0] = Tp[i]
                    Tpdummy[i][1] = thetai[i]
                np.savetxt('Tpt.txt', Tpdummy, header= "the wavelength used is %.5e meters. First column is Transmittance of p waves, second is angle" %wl[0])
            else:
                pass

    elif thetai[0] == thetai[1] and wl[0] == wl[1]: # If both are constants. Could have done all of this before, but this way it's faster
        #Defining variebles for each layer. Each layer needs to have a number of points with the length of "size"
        theta = zeros(len(n)).tolist() #ith element gives the angle at the ith layer
        theta[0] = thetai 
        phi = zeros(len(n)).tolist() #This will have the phase of each intermediate layer
        #Making constant refractive indices into same length arrays and finding angles and phases at different layers
        for i in range(len(n)):
            n[i] = n[i]*ones(size) #To make sure constant n values don't bug 
            theta[i] = array(list(map(cmath.asin,(n[0]*sin(thetai)/n[i]))))
            phi[i] = 2*pi*n[i]*d[i]*cos(theta[i])/wl #since thickness is for initial and final is zero, their phi is also zero
        #Create dummy arrays to have dynamical and propagation matrices
        Ds = zeros(len(n)).tolist()
        Dp = zeros(len(n)).tolist()
        P = zeros(len(n)).tolist()
        Kls = zeros(len(n)-2).tolist() #This serves as a dummy array to get D*P*D^-1 for each layer
        Klp = zeros(len(n)-2).tolist() 
       #Create dummy arrays to have have the reflectance/transmittance 
        Rs = zeros(size)
        Ts = zeros(size)
        Rp = zeros(size)
        Tp = zeros(size)
        #Actual maths
        i = 0 #Becasue wl and theta are constant, all the points are the same, we just choose the first
        for j in range(len(n)): #loops over all media
            #Dynamical Matrices for s/TE waves
            Ds[j] = matrix([[1,1],[n[j][i]*cos(theta[j][i]),-n[j][i]*cos(theta[j][i])]])
            #Dynamical Matrices for p/TM waves  
            Dp[j] = matrix([[cos(theta[j][i]),cos(theta[j][i])],[n[j][i],-n[j][i]]])
            #Propagation Matrix
            P[j] = matrix([[exp(1j*phi[j][i]),0],[0,exp(-1j*phi[j][i])]]) 
        for j in range(len(n)-2): #need to exclude initial and final 
            #dummy intermediate arrays
            Kls[j] = Ds[j+1]*P[j+1]*Ds[j+1].I
            Klp[j] = Dp[j+1]*P[j+1]*Dp[j+1].I
        Ks = reduce(mul,Kls)
        Kp = reduce(mul,Klp)
        #transfer matrix
        Ms = Ds[0].I*Ks*Ds[-1]
        Mp = Dp[0].I*Kp*Dp[-1]
        #Reflectance and transmittance for s/TE waves
        Rs[i] = (abs(Ms[1,0]/Ms[0,0]))**2
        Rs = Rs[i]
        Ts[i] = real(n[-1][i]*cos(theta[-1][i]))/real(n[0][i]*cos(theta[0][i]))*(abs(1/Ms[0,0]))**2
        Ts = Ts[i]
        #Reflectance and transmittance for p/TM waves
        Rp[i] = (abs(Mp[1,0]/Mp[0,0]))**2
        Rp = Rp[i]
        Tp[i] = real(n[-1][i]*cos(np.conjugate(theta[-1][i])))/real(n[0][i]*cos(np.conjugate(theta[0][i])))*(abs(1/Mp[0,0]))**2
        Tp = Tp[i]
        if options.get("text") == "on":
            print ("For an incident angle of", theta[0][0].real*180/pi, "and a wavelength of", wl[0], "meters, the results are:")
            print ("Rs=", Rs)
            print ("Rp=", Rp)
            print ("Ts=", Ts)
            print ("Tp=", Tp)
        else:
            pass
    
    elif thetai[0] != thetai[1] and wl[0] != wl[1]: #if both are ranges
        thetai = linspace(min(thetai), max(thetai), sizetheta) #Making the angle array much smaller, otherwise code is too slow
        theta = zeros(len(thetai)).tolist()
        phi = zeros(len(thetai)).tolist()
        #Create dummy arrays to have dynamical and propagation matrices, this time they have the length of the number of angles
        Ds = zeros(len(thetai)).tolist()
        Dp = zeros(len(thetai)).tolist()
        P = zeros(len(thetai)).tolist()
        Kls = zeros(len(thetai)).tolist() #This serves as a dummy array to get D*P*D^-1 for each layer
        Klp = zeros(len(thetai)).tolist() 
       #Create dummy arrays to have have the reflectance/transmittance for each angle
        Rs = zeros(len(thetai)).tolist()
        Ts = zeros(len(thetai)).tolist()
        Rp = zeros(len(thetai)).tolist()
        Tp = zeros(len(thetai)).tolist()
        
        #To make sure constant n values don't bug
        for i in range(len(n)):
            n[i] = n[i]*ones(size) 
        # A nice message
        if options.get("text") == "on":
            print ("\n")
            print ("Please wait, I'm doing all the matrix calculations now!")
        else:
            pass
        for t in range(len(thetai)): #loop over all angles 
            #Defining variebles for each layer. Each layer needs to have a number of points with the length of "size"
            theta[t] = zeros(len(n)).tolist() #ith element gives the angle at the ith layer
            theta[t][0] = thetai[t] 
            phi[t] = zeros(len(n)).tolist() #This will have the phase of each intermediate layer
            #Making constant refractive indices into same length arrays and finding angles and phases at different layers
            for i in range(len(n)): #loop over all layers
                theta[t][i] = array(list(map(cmath.asin,(n[0]*sin(thetai[t])/n[i]))))
                phi[t][i] = 2*pi*n[i]*d[i]*cos(theta[t][i])/wl #since thickness is for initial and final is zero, their phi is also zero
            #Create dummy arrays to have dynamical and propagation matrices
            Ds[t] = zeros(len(n)).tolist()
            Dp[t] = zeros(len(n)).tolist()
            P[t] = zeros(len(n)).tolist()
            Kls[t] = zeros(len(n)-2).tolist() #This serves as a dummy array to get D*P*D^-1 for each layer
            Klp[t] = zeros(len(n)-2).tolist() 
           #Create dummy arrays to have have the reflectance/transmittance 
            Rs[t] = zeros(size)
            Ts[t] = zeros(size)
            Rp[t] = zeros(size)
            Tp[t] = zeros(size)
            #Actual maths
            for i in range(size): #loops over all points 
                for j in range(len(n)): #loops over all media
                    #Dynamical Matrices for s/TE waves
                    Ds[t][j] = matrix([[1,1],[n[j][i]*cos(theta[t][j][i]),-n[j][i]*cos(theta[t][j][i])]])
                    #Dynamical Matrices for p/TM waves  
                    Dp[t][j] = matrix([[cos(theta[t][j][i]),cos(theta[t][j][i])],[n[j][i],-n[j][i]]])
                    #Propagation Matrix
                    P[t][j] = matrix([[exp(1j*phi[t][j][i]),0],[0,exp(-1j*phi[t][j][i])]]) #it's giving "nan", might be a problem?
                for j in range(len(n)-2): #need to exclude initial and final 
                    #dummy intermediate arrays
                    Kls[t][j] = Ds[t][j+1]*P[t][j+1]*Ds[t][j+1].I
                    Klp[t][j] = Dp[t][j+1]*P[t][j+1]*Dp[t][j+1].I
                Ks = reduce(mul,Kls[t])
                Kp = reduce(mul,Klp[t])
                #transfer matrix
                Ms = Ds[t][0].I*Ks*Ds[t][-1]
                Mp = Dp[t][0].I*Kp*Dp[t][-1]
                #Reflectance and transmittance for s/TE waves
                Rs[t][i] = (abs(Ms[1,0]/Ms[0,0]))**2
                Ts[t][i] = real(n[-1][i]*cos(theta[t][-1][i]))/real(n[0][i]*cos(theta[t][0][i]))*(abs(1/Ms[0,0]))**2
                #Reflectance and transmittance for p/TM waves
                Rp[t][i] = (abs(Mp[1,0]/Mp[0,0]))**2
                Tp[t][i] = real(n[-1][i]*cos(np.conjugate(theta[t][-1][i])))/real(n[0][i]*cos(np.conjugate(theta[t][0][i])))*(abs(1/Mp[0,0]))**2
                
                """
                Here, R and T will no longer give a 1 x size array, but rather 
                a size theta x size array, which makes it much slower to plot.
                """
        "Contour Plots"
        if options.get("plot") == "on":
            if "cmap" in options: #Need to fix this
                color = "Blues"  #Type of contour plot
            else:
                color = 'viridis'#Type of contour plot
            numberofcontours = 100  #number of colours we want in the contourplot #FIX ME: The number, might need change depending on what we want
            #Rs
            plt.figure()
            plt.contourf(wl, thetai, Rs, numberofcontours, cmap = color) 
            plt.colorbar();
            plt.xlabel('${\lambda}/m$') 
            plt.ylabel('${\Theta}/rad$') 
            plt.title('Rs')
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #force exponential in x axis
            #Ts
            plt.figure()
            plt.contourf(wl, thetai, Ts, numberofcontours, cmap = color) 
            plt.colorbar();
            plt.xlabel('${\lambda}/m$') 
            plt.ylabel('${\Theta}/rad$') 
            plt.title('Ts')
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #force exponential in x axis
            #Rp
            plt.figure()
            plt.contourf(wl, thetai, Rp, numberofcontours, cmap = color) 
            plt.colorbar();
            plt.xlabel('${\lambda}/m$') 
            plt.ylabel('${\Theta}/rad$') 
            plt.title('Rp')
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #force exponential in x axis
            #Tp
            plt.figure()
            plt.contourf(wl, thetai, Tp, numberofcontours, cmap = color) 
            plt.colorbar();
            plt.xlabel('${\lambda}/m$') 
            plt.ylabel('${\Theta}/rad$') 
            plt.title('Tp')
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #force exponential in x axis
        else:
            pass
        
        "Saving arrays with results"
        if options.get("file") == "on":
            #Rs
            Rsdummy = array(zeros((len(thetai)+1,len(wl)+1))) #+1 to have the wl and theta arrays
            for j in range(len(wl)):
                Rsdummy[0][j+1] = wl[j] #First row but skip first value
            for i in range(len(thetai)):
                Rsdummy[i+1][0] = thetai[i] #First column but skip first value
                for j in range(len(wl)):
                    Rsdummy[i+1][j+1] = Rs[i][j]
            np.savetxt('Rs.txt', Rsdummy, header= " wavelength = Rs[0,1:] \n angle = Rs[1:,0] \n matrix with s reflection = Rs[1:,1:]")

            #Rp
            Rpdummy = array(zeros((len(thetai)+1,len(wl)+1))) #+1 to have the wl and theta arrays
            for j in range(len(wl)):
                Rpdummy[0][j+1] = wl[j] #First row but skip first value
            for i in range(len(thetai)):
                Rpdummy[i+1][0] = thetai[i] #First column but skip first value
                for j in range(len(wl)):
                    Rpdummy[i+1][j+1] = Rp[i][j]
            np.savetxt('Rp.txt', Rpdummy, header= " wavelength = Rp[0,1:] \n angle = Rp[1:,0] \n matrix with p reflection = Rp[1:,1:]")

            #Ts
            Tsdummy = array(zeros((len(thetai)+1,len(wl)+1))) #+1 to have the wl and theta arrays
            for j in range(len(wl)):
                Tsdummy[0][j+1] = wl[j] #First row but skip first value
            for i in range(len(thetai)):
                Tsdummy[i+1][0] = thetai[i] #First column but skip first value
                for j in range(len(wl)):
                    Tsdummy[i+1][j+1] = Ts[i][j]
            np.savetxt('Ts.txt', Tsdummy, header= " wavelength = Ts[0,1:] \n angle = Ts[1:,0] \n matrix with s transmission = Ts[1:,1:]")

            #Tp
            Tpdummy = array(zeros((len(thetai)+1,len(wl)+1))) #+1 to have the wl and theta arrays
            for j in range(len(wl)):
                Rsdummy[0][j+1] = wl[j] #First row but skip first value
            for i in range(len(thetai)):
                Rsdummy[i+1][0] = thetai[i] #First column but skip first value
                for j in range(len(wl)):
                    Tpdummy[i+1][j+1] = Tp[i][j]
            np.savetxt('Tp.txt', Tpdummy, header= " wavelength = Tp[0,1:] \n angle = Tp[1:,0] \n matrix with p transmission = Tp[1:,1:]")

        else:
            pass
    if wl[0] == wl[2]: #use index 2 just because
        wl = wl[0]
    if thetai[0] == thetai[2]:
        thetai = thetai[0]
    return Rs, Ts, Rp, Tp, wl, thetai #need to fix this to give a single wl and thetai when necessary
    
"Electric Field"""""""""""""""""""""""""""""""""""""""""""""""""""       
def efield(size, wl, thetai, md, **options):
    """
    >This function finds the electric field and magnetic field at all points of
    all layers.
    >The first part of this code is exactly like the previous.
    >size is the resolution (number of points) the user wants on each individual
    layer, not the total!!
    >wl is the (constant) wavelenght wanted
    >theta is the (constant) angle of incidence wanted
    >>in md, insert an array with 2 arrays inside, one for refractive indices,
    the other for thicknesses. All of these should be inserted in order. Each 
    one of these can be simplified by multiplying repeated layers. For example,
    an optical cavity from Si3N4 and SiO2: 
        
        size = 500
        wl = 650e-9
        thetai = 30
        n1 = [1e-06, 'Si3N4_real.dat'] #Silicon Nitride
        n2 = [1e-06, 'SiO2_real.dat'] #Silica or silicon dioxide
        #need to get the refractive index at a specific wavelength
        n1w = refractive_index(wl, size, 1e-6, "Si3N4_real.dat") #Silicon Nitride
        n2w = refractive_index(wl, size, 1e-6, "SiO2_real.dat") #Silica or silicon dioxide
        d1 = wl/(4*real(n1w))
        d2 = wl/(4*real(n2w))
        
        m = [1,[n2,n1]*10,n2,[n1,n2]*10,1]
        d = [[d2,d1]*10,4*d2,[d1,d2]*10]  
        md = [m,d]
        E, H, x, h = efield(size, wl, thetai, md, plot = "on")
        
    >plot = "on" in the options to enable plots
    >Returns E-field and H-field along spatial coordinate x. h is the total 
    thickness of all the layers combined
    """
    thetaidummy = thetai #This is just to say which angle is in being used in the plot
    wldummyplot = wl #Just to say which wavelength is used in the plot
    #First lets deal with the angle
    if isinstance(thetai, (list, tuple, np.ndarray)) is False:
        thetai = linspace(thetai,thetai,size)*pi/180
    elif len(thetai) == 1:
        thetai = linspace(thetai[0],thetai[0],size)*pi/180
    else:
        sys.exit("For the electric field, you must provide a constant wavelength and angle")
    #check if wavelength is constant
    if isinstance(wl, (list, tuple, np.ndarray)) is False or len(wl) == 1:
        pass
    else:
        sys.exit("For the electric field, you must provide a constant wavelength and angle")
    #Check that incident angle is in the right place    
    if max(thetai) > pi/2 or min(thetai) < 0:
        sys.exit("The incident angle must be between 0 and 90 degrees")
    else:
        pass
    
    #extract information from the m arguments
    m = md[0]
    d = md[1]
    #Getting the n separated per layer
    dummy = [1,2,3] #dummy array just to use in while loop
    #making constants in the list into lists to use a line below it
    while dummy[-1] > dummy[-2]: #stops when length stops increasing i.e. when it's all in a list
        for i in range(len(m)):
            if isinstance(m[i], (list, tuple, np.ndarray)) is False:
                m[i] = [m[i]]               
            else:
                if isinstance(m[i][1], str) is True:
                    m[i] = [m[i]]
                else:
                    pass
        #making list flatter
        m = [item for sublist in m for item in sublist]
        dummy.append(len(m))
    ndummy = m
    #Similarly for thickness, but this time there are no strings so it's simpler
    dummy = [1,2,3] #dummy array just to use in while loop
    #making constants in the list into lists to use a the line below it
    while dummy[-1] > dummy[-2]: #stops when length stops increasing i.e. when it's all in a list
        for i in range(len(d)):
            if isinstance(d[i], (list, tuple, np.ndarray)) is False:
                d[i] = [d[i]]               
            else:
                pass                     
        #making list flatter
        d = [item for sublist in d for item in sublist]
        dummy.append(len(d))
    d = [0] + d + [0] #adding initial and final zero thickness
    dt = array(d) #defining d as array to use in ticks
    #Now, we extract the refractive indices and wavelengths we are using
    wldummy = zeros(len(ndummy)).tolist()
    n = zeros(len(ndummy)).tolist()
    for j in range(len(ndummy)):
        if isinstance(ndummy[j], (list, tuple, np.ndarray)) is True:
            if len(ndummy[j]) == 3:
                n[j], wldummy[j] = refractive_index("", size, ndummy[j][0], ndummy[j][1], ndummy[j][2])
            else:
                n[j], wldummy[j] = refractive_index("", size, ndummy[j][0], ndummy[j][1])
        #This will give the refractive index and wavelength. The ndummy[j][n] part gives units, real, and imaginary part respectively
        else:
            n[j], wldummy[j] = ndummy[j], 1 #this gives the constant value for refractive index and an innocuous value of 1 to wl to facilitate the use later 
    ###########################################################################
    #Finding the wavelength range at which we can work
    """
    >In this part, we first check if the dummy array (wldummy) containing all
    the wavelengths has any array in it. Each one of this arrays will be a 
    constraint to the wavelenght we might want to use.
    >Each ith value of wldummy with an array will give the ith value of 
    dummydummy a value of 1, or 0 if it has no array.
    >If dummydummy has no 1's, then there are no constraints and we can just use
    the desired wavelength.
    >If there is only one constraint, then we just find the overlap between the 
    desired wavelength and the constraint
    >If there is more than one, we loop and find the overlap between the desired 
    wavelength and all constraints.

    """
    dummydummy = zeros(len(wldummy)) 
    if isinstance(wl, (list, tuple, np.ndarray)) is False: #if we care only about a constant value, need to make this to avoid error in the overlap process
        wl = [wl,wl]
    else:
        pass
    for i in range(len(wldummy)): #To get an array of 0's and 1's 
        if isinstance(wldummy[i], (list, tuple, np.ndarray)) is True:
            dummydummy[i] = 1
        else:
            dummydummy[i] = 0
    if count_nonzero(dummydummy == 1) == 0:
        wl=linspace(min(wl),max(wl),size)
    elif count_nonzero(dummydummy == 1) == 1:
        dummydummy = dummydummy.tolist() #transform array into list to be able to find position of the value=1
        wl = overlap(wl,wldummy[dummydummy.index(1)],len(wldummy[dummydummy.index(1)])) 
    else:
        indices = [i for i, x in enumerate(dummydummy) if x == 1] #gives an array with the positions of all the elements equal to one in ndummy
        for j in range(len(indices)):
            wl = overlap(wl,wldummy[indices[j]],len(wldummy[indices[j]]))
    #Need to change the n's to match in wavelength we can work in
    for i in range(len(n)):   
        if isinstance(n[i], (list, tuple, np.ndarray)) is True:
            n[i] = interp(wl,wldummy[i],n[i]) #interpolation that cuts down the n's to the wavelength we're using
        else:
            pass
    
    #Error message if user has put the wrong m, n combination
    if len(m) != len(d):
        sys.exit("Your refractive index vector must have the same length as your as your thicknesses vector plus two")
    else:
        pass

    #NEW PART##################################################################
    #Create dummy arrays, for the electric field distribution. nth element correspond to the nth layer constants
    As = zeros(len(n)).tolist() #Left to right movement
    Bs = zeros(len(n)).tolist() #Right to left movement
    Ap = zeros(len(n)).tolist() #Left to right movement
    Bp = zeros(len(n)).tolist() #Right to left movement

    k = zeros(len(n)).tolist() #Wavevector 
    theta = zeros(len(n)).tolist() #ith element gives the angle at the ith layer
    theta[0] = thetai 
    Es = zeros(len(n)).tolist() #E-field in each layer for s waves
    Ep = zeros(len(n)).tolist() #E-field in each layer for p waves
    d[0] = d[-1] = max(d) #defining thicknesses for initial and final layers
    ddummy = zeros(len(n)).tolist() #these will have the x values at each layer
    phi = zeros(len(n)).tolist() #this is needed to find the initial B
    for i in range(len(n)): #loops over all layers
        n[i] = n[i]*ones(size) #To make sure constant n values don't bug 
        ###Careful with arcsin ambiguity!! This might be a problem with complex refractive indices
        theta[i] = array(list(map(cmath.asin,(n[0]*sin(thetai)/n[i])))) #find angle at each layer
        k[i] = 2*pi*n[i]*cos(theta[i])/wl #wl gives "size" points to each layer
        ###CAREFUL HERE, IN THE BOOK IT IS DIFFERENT, BUT A PAPER SOLVING THE INSTABILITIES WRITES IT LIKE THIS
        ddummy[i] = linspace(0, d[i], size) #an x between the initial and final part of the layer: FINAL OR INITIAL LAYER MAYBE WRONG
        phi[i] = k[i]*d[i] #Need to fix initial and final phi
        ###Need to correct TMM instabilities
#        if imag(n[i][0]) < 0: #carefull, we choose [0] because we know that the wavelength is constant
#            k[i] = abs(real(k[i])) + 1j*abs(imag(k[i])) #set up Im(k) > 0 nad Re(k) > 0.
#        else:
#            pass
    phi[0] = phi[-1] = zeros(len(wl)) #This shouldn't matter anyways
    #ddummy[-1] = linspace(0, d[-1], size) #According to book this is what it should be
    ddummy[0] = linspace(-d[0], 0, size) #According to paper fixing TMM, this is how it should be

    #Create dummy arrays to have dynamical and propagation matrices
    Ds = zeros(len(n)).tolist()
    Dp = zeros(len(n)).tolist()
    P = zeros(len(n)).tolist()
    Kls = zeros(len(n)-2).tolist() #This serves as a dummy array to get D*P*D^-1 for each layer
    Klp = zeros(len(n)-2).tolist() 
    Ks = zeros(len(n)).tolist() #This will have the product of all Kls up to the ith element in the ith element
    Kp = zeros(len(n)).tolist()
    Ms = zeros(len(n)).tolist() #This is the transfer matrix
    Mp = zeros(len(n)).tolist()
    
    #Actual maths##############################################################
    #Because wl and theta are constant, all the points are the same, we just choose the first
    for j in range(len(n)): #loops over all layers
        #Dynamical Matrices for s/TE waves
        Ds[j] = matrix([[1,1],[n[j][0]*cos(theta[j][0]),-n[j][0]*cos(theta[j][0])]])
        #Dynamical Matrices for p/TM waves  
        Dp[j] = matrix([[cos(theta[j][0]),cos(theta[j][0])],[n[j][0],-n[j][0]]])
        #Propagation Matrix
        P[j] = matrix([[exp(1j*phi[j][0]),0],[0,exp(-1j*phi[j][0])]]) 
    for j in range(len(n)-2): #need to exclude initial and final 
        #dummy intermediate arrays
        Kls[j] = Ds[j+1]*P[j+1]*Ds[j+1].I #this is DjPjDj^-1
        Klp[j] = Dp[j+1]*P[j+1]*Dp[j+1].I
    ###########################################################################
    #need to find reflectivity to find B[0]
    Ksd = reduce(mul,Kls)
    Kpd = reduce(mul,Klp)
    #transfer matrix
    Msd = Ds[0].I*Ksd*Ds[-1]
    Mpd = Dp[0].I*Kpd*Dp[-1]
    #Reflectivity for s/TE waves
    rs = Msd[1,0]/Msd[0,0]
    #Reflectivity for p/TM waves
    rp = Mpd[1,0]/Mpd[0,0]
    ###########################################################################
    #starting A and B
    As[0] = 1 
    Ap[0] = 1
    Bs[0] = rs*As[0]
    Bp[0] = rp*Ap[0]
    
    for i in range(len(n)-1): #skip first one to avoid bug since they have been defined. Use [0,0] to go from matrix to constant
        As[i+1] = (Ds[i+1].I*Ds[i]*P[i].I*matrix([[As[i]],[Bs[i]]]))[0][0,0]
        Bs[i+1] = (Ds[i+1].I*Ds[i]*P[i].I*matrix([[As[i]],[Bs[i]]]))[1][0,0]
        Ap[i+1] = (Dp[i+1].I*Dp[i]*P[i].I*matrix([[Ap[i]],[Bp[i]]]))[0][0,0]
        Bp[i+1] = (Dp[i+1].I*Dp[i]*P[i].I*matrix([[Ap[i]],[Bp[i]]]))[1][0,0]
    Bp[-1] = Bs[-1] = 0 #Typically, these are zero anyways, but to make sure, I added this
    
    H = zeros(len(n)).tolist() #Need to find the magnetic field to plot. Ep waves are discontinuous along the x axis
    c = scipy.constants.c #speed of light
    #And finally, what we're actually interested in 
    for i in range(len(n)):
        Es[i] = As[i]*exp(-1j*k[i]*ddummy[ i]) + Bs[i]*exp(1j*k[i]*ddummy[i]) #[0,0] because A and B elements are matrices
        H[i] = n[i]/(c*4*pi*10**(-7)) * (Ap[i]*exp(-1j*k[i]*ddummy[i]) - Bp[i]*exp(1j*k[i]*ddummy[i])) #[0,0] because A and B elements are matrices
        #Note that this should have a cos(theta[i]) but that messes up everything, but like this is perfect, I don't know why
    
    "Plots"
    #cartesian coordinate x
    xnew = zeros(len(n)).tolist()
    for i in range(len(n)):
        xnew[i] = sum(d[:i+1]) #sum up to the ith element
    xnew[0] = linspace(0, xnew[0], size)
    for i in range(len(n)-1): #skip first term
            xnew[i+1] = linspace(xnew[i][-1], xnew[i+1], size) #need [-1] to tell it's from last element
    interface = max(xnew[0]) #going to push left every layer by this
    for i in range(len(n)): #making zero at first interface
        xnew[i] = xnew[i] - interface
        
    #give option to plot graphs
    if options.get("plot") == "on":
        #s waves###########
        plt.figure()
        for i in range(len(n)): #plot for each layer
            plot(xnew[i], abs(Es[i])**2,'r',linewidth=1.5)
        if len(dt) > 10: #If it has too many layers, things will get too messy
            for i in range(len(dt)):
                plt.axvline(dt[0:i].sum(),color='gray',linestyle='dashed',linewidth=1)
        else:
            for i in range(len(dt)):
                plt.axvline(dt[0:i].sum(),color='gray',linestyle='dashed',linewidth=1.5)
            plt.grid()  
        plt.gca().legend(("TE (s)",))
        plt.title('Electric Field E')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #force exponential in x axis
        plt.figtext(0.5, 0, "${\Theta}_i$ = %.2f$^\circ$" %thetaidummy, position=(0.137,0.925), fontsize = 12, bbox={"facecolor":"orange", "alpha":0.5, "pad":5}) 
        plt.figtext(0.5, 0, "${\lambda}_i$ = %.2e $m$" %wldummyplot, position=(0.7,0.925), fontsize = 12, bbox={"facecolor":"orange", "alpha":0.5, "pad":5}) 
        plt.ylabel('$|E_{\mathrm{TE}}|^2$')
        plt.xlabel("$x/m$")
        
#        plt.axvspan(sum(d[:17]) - interface, sum(d[:18]) - interface, facecolor='0.2', alpha=0.5) #grey area for specific x values
#        plt.axvspan(sum(d[:18]) - interface, sum(d[:19]) - interface, facecolor='y', alpha=0.5) #yellow area for specific x values
        #see if we need to for exponential in y axis too
        sd = zeros(len(n))
        for i in range(len(n)):
            sd[i] = max(abs(Es[i])**2)
        if max(sd) < 0.1:
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) #force exponential in x axis
        else:
            pass
        #p waves#############
        plt.figure()
        for i in range(len(n)): #plot for each layer
            plot(xnew[i], abs(H[i])**2,'r',linewidth=1.5)
        if len(dt) > 10: #If it has too many layers, things will get too messy
            for i in range(len(dt)):
                plt.axvline(dt[0:i].sum(),color='gray',linestyle='dashed',
                            linewidth=1)
        else:
            for i in range(len(dt)):
                plt.axvline(dt[0:i].sum(),color='gray',linestyle='dashed',
                            linewidth=1.5)
            plt.grid()  
        plt.gca().legend(("H (p)",))
        plt.title('Magnetic Field H')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) #force exponential in x axis
        plt.figtext(0.5, 0, "${\Theta}_i$ = %.2f$^\circ$" %thetaidummy, position=(0.137,0.925), fontsize = 12, bbox={"facecolor":"orange", "alpha":0.5, "pad":5}) 
        plt.figtext(0.5, 0, "${\lambda}_i$ = %.2e $m$" %wldummyplot, position=(0.7,0.925), fontsize = 12, bbox={"facecolor":"orange", "alpha":0.5, "pad":5}) 
        plt.ylabel('$|H_{\mathrm{TM}}|^2$')
        plt.xlabel('$x/m$')
        
#        plt.axvspan(sum(d[:17]) - interface, sum(d[:18]) - interface, facecolor='0.2', alpha=0.5) #grey area for specific x values
#        plt.axvspan(sum(d[:18]) - interface, sum(d[:19]) - interface, facecolor='y', alpha=0.5) #yellow area for specific x values
        #see if we need to for exponential in y axis too
        sd = zeros(len(n))
        for i in range(len(n)):
            sd[i] = max(abs(H[i])**2)
        if max(sd) < 0.1:
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) #force exponential in y axis
        else:
            pass
    else:
        pass #doesn't plot graphs
    return concatenate(Es), concatenate(H), concatenate(xnew), sum(d[:17]) - interface #concatenate joins all arrays in a list

#Need to fix the bug that changes m, d, and md after running this code