# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 15:47:24 2016

@author: Elena Koumpli


"""

import numpy as np
from scipy.optimize import newton,fsolve
#from numba import jit - this doesn't work with current python installation for some reason-should work with Anaconda
import warnings

from pvlib import pvsystem


class PVDevice(object):
    
    '''
    PV cell builder: this module is used to build a cell with 
    characteristics as extracted and/or taken
    from manufaturer's data sheets. It can also be used for 
    modules if the number of cells in series (Ns) is other than 1. 
    
    *Initially the device is built for STC, G =1000W/m2 and T =25C or 298.0 K*
    
    *NOTE*
    
    Make sure you enter the right units for temperature coefficients, 
    as they may given as absolute or percentage values
    
    References
    -----------
    De Soto, W.,Klein, S.a. Beckman, W. "Improvement and validation of a model for photovoltaic 
    array performance" Solar energy, 2006
    
    Villalva, M, Gazoli, J,Filho, E "Comprehensive Approach to Modeling and Simulation of Photovoltaic Arrays"
    IEEE Transaction on Power Electronics, 2009
    
    '''
    
    
    def __init__(self, Rs, Rsh, I0_stc, IL_stc,  n, Ki, Kv,  Isc_stc, Voc_stc, Ns =1,Tcell=298.0, 
                 Gin = 1000, abr = 0, Vbr = None, m=0, change_Rsh=False,
                 Eg0 =1.16, k =1.38E-23 , q = 1.602E-19):
        
        
        '''
        Parameters
        -----------
        
        
        NOTE: initially the device is defined at Standard Testing Conditions (STC). 
        
        
        
        Rs      :  series resistance 
        Rsh     :  shunt resistance 
        I0_stc  :  diode saturation current 
        IL_stc  :  light current (photocurrent)
        Ki      :  temperature coefficient for current in A/Celsius
        Kv      :  temperature coefficient for voltage in A/Celsius
        n       :  diode ideality factor
        change_Rsh: If true Rsh changes with irradiance 
        Isc      :  short circuit current
        Voc_stc :  open circuit voltage 
        Tcell   :  cell temperature in Kelvin (298 degrees Kelvin at STC)
        Gin     :  irradiance (1000 W/m2 at STC)
        Ns      :  number if cells in series
        Eg0     :  Energy gap (not used at the moment)
        k  = 1.38E-23    : Boltzmann constant   
        q = 1.602E-19  : electric charge
        
        
        Breakdown properties: a,b and Vbr
        
        Vbr : breakdown voltage
        
        '''
        self.Rs = Rs
        self.Rsh = Rsh
        self.Rsh_stc = Rsh
        self.I0 = I0_stc
        self.IL = IL_stc
        self.I0_stc = I0_stc
        self.IL_stc = IL_stc
        self.Ki = Ki
        self.Kv = Kv
        self.n = n
        self.Isc = Isc_stc
        self.Isc_stc = Isc_stc
        self.Voc_stc = Voc_stc
        self.Voc = Voc_stc
        self.Tcell = Tcell
        self.Gin = Gin
        self.Vbr = Vbr
        self.abr = abr
        self.m = m
        self.k = k
        self.q = q
        self.Ns = Ns
        self.Eg0 = Eg0
        self.change_Rsh = change_Rsh
        
        if Vbr is None:
            self.Vbr = (-3)*self.Voc
        else:
            self.Vbr = Vbr
            
        self.a = self.n*self.Ns*self.k*self.Tcell/self.q
        
    
    def __repr__(self):
        
        return ('A device with: Rseries=' + str(self.Rs) +
                ', Rshunt = ' + str(self.Rsh) +
                ', ideality factor = ' + str(self.n) +
                ',\n photocurrent = ' + str(self.IL) + 
                #',\n open circuit voltage (STC) = ' + str(self.Voc) +
                #',\n short circuit current  = ' + str(self.Isc) +
                ', diode current  = ' + str(self.I0) +
                ',\n under irradiance =  ' + str(self.Gin) +
                ' W/m2 and temperature = ' + str(self.Tcell)+'K and' +
                ' breakdown voltage = ' +str(self.Vbr)+'Volts and '+
                'b = ' +str(self.abr) +'and m = '+str(self.m))
    
    
    def get5params(self):
        
        
        print('Rshunt: '+str(round(self.Rsh,6))+ '\nRseries: '+str(round(self.Rs,6))+
              '\nn: '+str(round(self.n,2)) + '\nIph: '+str(round(self.IL,3)) +
              '\nIo: '+str(self.I0) +' and with number of cells: '+ 
              str(self.Ns))
    
      
    def setEnvironment(self,Tcell,Gin):
    
        # 8 Kelvin is the limit for the simulation i.e. -265 Celsius 
        
        kev=8.167E-5
        
        t = Tcell
        deltaT = t-298.0  
        
        # modified ideality factor
        self.a = self.a*t/self.Tcell # do not confuse self.Tcell with Tcell, which is the argument of this function
        
        # photocurrent
        IL_n = (self.IL_stc +self.Ki*deltaT)*Gin/1000 
        
        # short circuit current
        Isc_n = (self.Isc_stc +self.Ki*deltaT)*Gin/1000
        
        # temperature difference from STC
        if deltaT == 0:
            
            I0_n= self.I0 #no change
        
        else:
            
            I0_n = (self.Isc_stc + self.Ki*deltaT)/(np.exp((self.Voc_stc+ self.Kv*deltaT)/self.a)-1) #Villalva et al.
            # more research needs to go in this equation as when deltaT =0 it doesn't return the initial value of I0 for
            # T= 25C and G=1000W/m2
        
        # some papers suggest a change in the Rsh for more accurate prediction 
        # haven't seen this myself in c-Si though
       
        if self.change_Rsh:
            
            Rsh_n = (1000/Gin)*self.Rsh_stc 
        else:
            Rsh_n = self.Rsh_stc 
        
        #This is blocked out as weird behaviour is seen near Voc: Eg for crystalline silicon, seems not so accurate
        #=======================================================================
        # Eg = self.Eg0 - ((0.000702*t**2) / (t+1108))  
        # tr = t/298.0
        # I0_n = (self.I0*(tr**3))*np.exp(((Eg*self.q)/(kev*self.n))*((1/298.0)-(1/t)))
        #=======================================================================

    
        # Changed parameters
        
        self.IL = IL_n
        self.I0 = I0_n
        self.Isc = Isc_n
        #self.Voc = Voc_n 
        self.Gin = Gin
        self.Tcell = Tcell
        self.Rsh = Rsh_n
        
       
    def IVcurveI(self,i,v):
        
         
        func = self.IL - self.I0*(np.exp((v+i*self.Rs)/self.a)-1) - (v+i*self.Rs)/self.Rsh - i - \
                    ((v+i*self.Rs)/self.Rsh)*self.abr*((1-(v+i*self.Rs)/self.Vbr)**(-self.m))
    
        return func   
        
     
    def IVcurveV(self,v,i):
        
        
        # same as IVcurvei but with voltage as the 1st parameter
        func = self.IL - self.I0*(np.exp((v+i*self.Rs)/self.a)-1) - (v+i*self.Rs)/self.Rsh - i - \
                    ((v+i*self.Rs)/self.Rsh)*self.abr*((1-(v+i*self.Rs)/self.Vbr)**(-self.m))
    
        return func     
     
    
    # 1st derivative for Voltage(V)
    def IVprimeV(self,v,i): 
    
        funcvp  = -(self.I0/self.a)*np.exp((v+i*self.Rs)/self.a) - (1/self.Rsh) -\
                    (self.abr/self.Rsh)*(1-(v+i*self.Rs)/self.Vbr)*(1+(self.m/self.Vbr)*(v+i*self.Rs)*\
                        (1-(v+i*self.Rs)/self.Vbr)**(-self.m-2))
                    
    
        return funcvp
    
    
    # 1st derivative for current (I)
    
    def IVprimeI(self,i,v):
        
        funcip =  -(self.I0*self.Rs/self.a)*np.exp((v+i*self.Rs)/self.a) - self.Rs/self.Rsh - 1 -\
                - (self.abr*self.Rs/self.Rsh)*((1-(v+i*self.Rs)/self.Vbr)**(-self.m))*\
                (1+self.m/((self.Vbr/(v+i*self.Rs))-1))
    
        return funcip
    
    
    
    def solveForCurrent(self, guess =0.0, start=None, stop =0.0, step = -0.1, vdata=None,  method = 'newton'):
        
        # Chooses values for current and finding voltage on the curve
        
        if (start is None and vdata is None): #when no data is given and no start is given either
            start = 1.1*self.Voc
      
        if (self.Vbr is not None and  np.abs(stop) > np.abs(0.80*self.Vbr)):
            print('Warning: stopping voltage is very close to the critical voltage of Cell with properties:\n')
            print(vars(self))
            
        if vdata is None:#if data are given then no array needs to be specified
            
            vspace = np.arange(start,stop,step)

        else:
            vspace = vdata
        
        ilist = []
        vlist = []
        
        guess_value = guess
        
        for v in vspace:
            
            if method == 'newton':  
                
                i = newton(self.IVcurveI, args =(v,), x0 = guess_value, fprime=self.IVprimeI,  maxiter = 3000)
            
            else:
            
                i = fsolve(self.IVcurveI,args = (v,), x0 = guess_value)
                i = i[0]
            
            vlist.append(v)
            ilist.append(i)
            
        voltage = np.array(vlist)
        current = np.array(ilist)       
        return  voltage, current
        
    
    def solveForVoltage(self, guess=None, start=0.0, stop =None, step = 0.01, idata=None,  method = 'newton'):
        
        ''' 
        This works more efficiently when starting from 0.0 
        and incrementing towards IL where guess value is Voc
        
        '''
        
        if stop is None and idata is None:
            stop = 1.1*self.Isc
        
        if guess is None and idata is None:
            guess = self.Voc_stc

        
        ilist = []
        vlist = []
        
        if idata is None:
            
            ispace = np.arange(start,stop,step)
        else:
            ispace = idata
        
        guess_value = guess
        
        for i in ispace:
            
            if method == 'newton':  
            
                v = newton(self.IVcurveV, args =(i,), x0 = guess_value, fprime=self.IVprimeV,  maxiter = 3000)
            
            else:
            
                v = fsolve(self.IVcurveV, args = (i,), x0 = guess_value)
                v = v[0]
            
            vlist.append(v)
            ilist.append(i)
            
        
        voltage = np.array(vlist)
        current = np.array(ilist)
        
        return voltage, current
        
    
    def findIVmax(self):
    
    
    # based on fsolve find current, voltage at maximum power point on curve 
        
        def _maxpoint(p):
            
            i,v = p
            
            func1 = self.IVcurveI(i,v)
            num = self.IVprimeV(v,i)
            den = self.IVprimeI(i,v)
        
            func2 = i-v*(num/den)
            
            z = np.array([func1,func2])
            return z
            
        # based on minpack hybrid algorithms
        
        x0 = [self.Isc,self.Voc]
        
        
        max_point = fsolve(_maxpoint,x0=x0)
        i_max = max_point[0] 
        v_max = max_point[1]
        p_max = i_max*v_max
        
        return i_max, v_max, p_max

    # an estimation of the new Voc can be given with this function using the analytical IV expressions
    # another method to be added is the graphical method which is probably more accurate
    
    def findVoc(self):
        
        def _vocpoint(p):
            
            v,i = p
            
            func1 = self.IVcurveV(v,i)
        
            func2 = 0
            
            z = np.array([func1,func2])
            
            return z
        
        voc_point = fsolve(_vocpoint,x0=[self.Voc,0])
        
        
        return voc_point
    















  

 
