
"""
Created on Tue Oct 11 12:26:25 2016

@author: Elena Koumpli

Parameter extraction from manufacturer's datasheets 


Reference:

Villalva et al. 2009 "Comprehensive approach to modelling and simulation of photovoltaic arrays",
IEEE transactions in power electronics 

"""

import PVDevice
import numpy as np
import pandas as pd
import os

def extractParams(n,Vmp_n,Imp_n, Voc_n,
                  Pmax, Isc_n, Ki,Kv, Ns,
                  T =298, Gin = 1000):
    
    '''
    Parameters
    -----------
        
        
    NOTE: initially the device is defined at Standard Testing Conditions (STC). 
        
        

        n       :  guess value for the diode ideality factor 
    
        manufacturer data:
        
        Ki      :  temperature coefficient for current in A/Celsius
        Kv      :  temperature coefficient for voltage in V/Celsius
        Imp_n   :  current at maximum power point
        Vmp_n   :  voltage at maximum power point
        Isc_n   :  short circuit current
        Voc_n   :  open circuit voltage 
        T       :  module temperature in Kelvin (298 degrees Kelvin at STC)
        Gin     :  irradiance (1000 W/m2 at STC)
        Ns      :  number if cells in series
        


    Returns
    -------

    The 5 parameters of the one diode model and the differences of the calculated values for the 3
    characteristic points of the IV curve from the manufacturer's data values.
    
    '''
    
    k  = 1.38E-23
    q = 1.602E-19
    
    plist = []
    
    # manufacturer data at nominal (n) conditions

    mif = n*Ns*k*T/q
    Ki = Ki*Isc_n/100
    Kv = Kv*Voc_n/100
    
    # Pmp = Pmax #this one may not be very accurate, as manufacturer's modify these values for the binning 
    # therefore, for more accurate representation Imax and Vmax are used to calculate pmax
    
    Pmp = Imp_n*Vmp_n
    
    # thermal voltage
    
    Vt_n = Ns*k*T/q
    

    # initial values 
    
    FF_n = Pmp/(Isc_n*Voc_n) 
    Rs = 0.0
    Rsh = Vmp_n/(Isc_n-Imp_n)-(Voc_n-Vmp_n)/Imp_n
    IL_n = Isc_n 
    Io_n = Isc_n/(np.exp(Voc_n/(n*Vt_n))-1) 
    
    
    # find initial pmax 
    
    eps = 0.001
    count = 0
    pmax = 0
    voc = Voc_n
    
    while (np.abs(pmax-Pmp)>eps):  
        
        
        mymod =PVDevice.PVDevice(Rs=Rs, Rsh=Rsh, I0_stc=Io_n, IL_stc=IL_n,  n=n, Ki=Ki, Kv=Kv,  Isc_stc=Isc_n, 
                         Voc_stc = Voc_n,Tcell=T, Gin = Gin, Ns =Ns, abr = 0, Vbr = None, m=0, change_Rsh=False,
                 Eg0 =1.16, k =1.38E-23 , q = 1.602E-19)
        
        im, vm, pmax = mymod.findIVmax()

        plist.append({'n':n, 'voc': voc,'Rs':Rs, 'Rsh':Rsh, 'pmax':pmax, 'eps': np.abs(pmax - Pmp),'IL':IL_n,'Io':Io_n})
        
        
        if count>1 and plist[count]['eps']>plist[count-1]['eps']:

            break
        
        if (count >4000):
            
            break
        
        count+=1

        # Assign new values

        Rs +=0.001
        voc = (im*vm)/(FF_n*Isc_n)
        IL_n = Isc_n*((Rs+Rsh)/Rsh)
        Io_n = IL_n/(np.exp(Voc_n/(n*Vt_n))-1)
        A = np.abs(Vmp_n*IL_n-Vmp_n*Io_n*np.exp((Vmp_n+Imp_n*Rs)*q/(Ns*n*k*T))+Vmp_n*Io_n-Pmp)
        Rsh = Vmp_n*(Vmp_n+Imp_n*Rs)/A 
    
  
    dicts = plist[count-1]
    
    
    irm = (im-Imp_n)/Imp_n
    vrm = (vm-Vmp_n)/Vmp_n
    pr = (pmax-Pmp)/Pmp
    vocr = (voc-Voc_n)/Voc_n
    
    dicts.update({'Irel':irm, 'Vrel':vrm, 'Prel':pr, 'Vocr':vocr}) 
    
    
    return dicts



def optim(**kwargs):

    # keywords are taken for the extractParams
    
    f = []
    for x in np.arange(1.0,1.9,0.01):
        
        dict = extractParams(n=x,**kwargs)
        
        ir = dict['Irel']
        vr = dict['Vrel']
        pr = dict['Prel']
        vocr = dict['Vocr']
        
        
        delta = np.sqrt(pr**2)
        
        dict.update({'delta': delta})
        f.append(dict)
    
    df = pd.DataFrame(f)
    optimum = df.loc[df['delta'] == df['delta'].min()]
    
    return optimum
    



def runFromFile(panel_id, csv_outfile = None):

    '''

    Runs the extraction based on the file containing manufacturer data for different panels.
    
    '''
    os.chdir('..')
     
    manufacturer_datafile = os.path.join(os.getcwd(),'Test Files//panel_info.csv')


    
    df = pd.read_csv(manufacturer_datafile, index_col = 'panel_id')

    try: 
        dfi = df.iloc[panel_id-1] 
      
        Vmp_n = dfi['Vmp']
        Imp_n = dfi['Imp']  
        Voc_n = dfi['Voc']
        Pmax=dfi['Pmax']
        Isc_n = dfi['Isc']
        Ki = dfi['a_isc (%/K)']
        Kv = dfi['b_voc (%/K)'] 
        Ns = dfi['Cells in series']
          
            
        params = optim(Vmp_n =Vmp_n ,Imp_n=Imp_n, Voc_n=Voc_n,
                              Pmax=Pmax, Isc_n=Isc_n,Ki=Ki,  Kv=Kv, Ns=Ns)    

        if csv_outfile:
           params.to_csv(csv_outfile)
        
        print(params)

    except IndexError:

        raise Exception('Index information for id = {} does not exist in the datasheet or there is not enough information available for that module.'.format(panel_id))
    
if __name__=='__main__':

    
    runFromFile(panel_id = 10)    

    
















