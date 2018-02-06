'''

Extract modelling parameters from experimental IV curves at STC using the one-diode moodel.
Works both for PV cells and modules

Main author(s)     : Elena Koubli

Test files by Michael-Owen Bellini (c-Si mini module) and Francesco Bittau (CdTe cell)

NOTE 1: At the moment fitting is based on the particular optimisation algorithm and its inherent error criterion alone 
and no extra error criterion is added here. 

NOTE 2:The trust-region optimisation has been found to outperform Levenberg-Marquardt 
in terms of convergence rates and succession and therefore it may be used without an extra error criterion (such area).

NOTE 3: The algorithm does take the IV as positive in the first quadrant. The algorithm returns best results if the supplied IV curves 
have more points in the *positive* voltage part than at negative voltage (towards breakdown). In the next update there will be an option 
to remove the negative part.Check the test files for your information.

NOTE 4: For seaborn graph styles check some options here: https://github.com/mwaskom/seaborn/issues/672

NOTE 5: Author is grateful for any feedback/forking provided.

Happy extraction!


'''

from PVDevice import PVDevice
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from pandas.io.common import CParserError
from plotGraphs import plotDynGraphsOnXY
import seaborn as sns
from statistical_functions import returnStats
import os


# set parameters

k = 1.38E-23
q = 1.602E-19
Tm=298.0


def fitcurve(V_points, I_points, IVfunc, initial): 
    
    '''
    Curve fitting function using the trust-region function, 
    Runs the scipy.optimize.curve_fit with certain conditions
    
    '''
   
    y_data = np.zeros(len(V_points)) #array of zeros
    I_V_points = np.array([V_points, I_points])
   
   

    lower = [0,0,0,0,0]# guess limits
    upper = [10, 1E-2, 1E4, 1E4, 35.0] #guess limits
    

    guess = curve_fit(IVfunc, I_V_points, y_data, p0=initial, 
                        bounds=(lower, upper), method='trf')#trust-region function
    
    
    
    return guess 







def main(input_file,cells_in_series,columns=['V','I'], output_file=None):

    '''
    
    About this function:
    1) prints out the 5 parameters n, I0, IL, Rs , Rsh.
    2) prints out the statistics for V, I, P at maximum power point and the whole IV curve
    
    Input parameters
    ----------------
    
    input_file:  file path, either a CSV or EXCEL file 
    
    cells_in_series: integer, number of cells in series 
                    *put the right number in otherwise ideality factor is weird*
    
    columns: list, names of voltage (column[0]) and current([0])  columns in the file 
              * if these are mixed there will probably be a convergence error *      
    
    output_file: file path, default None, a CSV file to save the simulated results 
    
    
    '''
    
    
    
    
    try:
        df = pd.read_csv(input_file, index_col = None)# reads CSV
    
    except CParserError:
        
        df = pd.read_excel(input_file) #reads excel
    
    
    # function used in fitting - one diode model 
    def _I_of_V(I_V_points,Iph,I0,Rs,Rsh,n,b=0,Vbr=100,m=0):
    
    
        V,I = I_V_points
        a = n*cells_in_series*k*Tm/q
        func = Iph - I0*(np.exp((V+I*Rs)/a)-1) - ((V+I*Rs)/Rsh)*(1+b*(1-((V+I*Rs)/Vbr))**(-m))- I 
        
        return func
  
    
    V_points = df[columns[0]].values
    I_points = df[columns[1]].values
    
    max_point = np.argmax(V_points*I_points)
    
    
    vmax = V_points[max_point]
    imax = I_points[max_point]
    pmax = vmax*imax
    
    #starting point for photo-current
    indexI=np.where(np.abs(V_points)==min(np.abs(V_points)))
    I_zero = I_points[indexI[0]]
    
    #starting point for voltage
    indexV=np.where(min(np.abs(I_points))==np.abs(I_points))

    if len(indexV[0]>1.0): #need to resolve when after the absolute you get two values corresponding to minimum voltage
        V_zero = float(V_points[indexV[0][0]])

    else:
        V_zero = float(V_points[indexV[0]]) 
    
    
    # setting guess values    
    guess_values ={'Iph':I_zero, 'I0':1E-8, 'Rs':0.01, 'Rsh':200, 'n':1,'abr':0.01,'m':1,'Vbr':-3*V_zero}
    
    
        
    initial = [guess_values['Iph'],guess_values['I0'],guess_values['Rs'],
               guess_values['Rsh'],guess_values['n']]
  
    
    
    fitted_params = fitcurve(V_points, I_points, _I_of_V, initial)[0]
    
    
    guess_values['Iph']= fitted_params[0]
    guess_values['I0'] = fitted_params[1]
    guess_values['Rs'] = fitted_params[2]
    guess_values['Rsh']= fitted_params[3]
    guess_values['n']  = fitted_params[4]
    
    
    
    
    mydevice = PVDevice(Rs= guess_values['Rs'], Rsh=guess_values['Rsh'], 
                        I0_stc=guess_values['I0'], Ki= None, Kv=None, Voc_stc = V_zero, Isc_stc = I_zero,
                        IL_stc=guess_values['Iph'], n=guess_values['n'], 
                          Ns =cells_in_series,Tcell=Tm)
    
    
    # this is optional testing block for fitting at different environmental conditions

    #===========================================================================
    # mydevice = PVDevice(Rs= 0.177, Rsh=190.2, 
    #                     I0_stc=1.34E-08, Ki= 0.047, Kv=-0.32, Voc_stc = V_zero, Isc_stc = I_zero[0],
    #                     IL_stc=8.688, n=1.2, Gin=1000,
    #                       Ns =cells_in_series,Tcell=Tm)
    #
    # #change the device curve here
    # mydevice.setEnvironment(Tcell=298.0, Gin=200)
    #===========================================================================
   
    
    mydevice.get5params() #gives you the 5 params as extracted from fitting
    
    
    v,i_sim = mydevice.solveForCurrent(vdata = V_points)
    
    
    max_point_sim = np.argmax(v*i_sim)
    
    
    vmax_sim = v[max_point_sim]
    imax_sim = i_sim[max_point_sim]
    pmax_sim = imax_sim*vmax_sim
    
  
    # stats at maximum power point
    
    imax_stats = returnStats(measured=imax, modeled=imax_sim)
    vmax_stats = returnStats(measured=vmax, modeled=vmax_sim)  
    pmax_stats = returnStats(measured=pmax, modeled=pmax_sim)  
    
    #print('pmax:',pmax,'pmax_sim:', pmax_sim)
    #print('vmax:',vmax,'vmax_sim:',vmax_sim) 
    print('Pmax stats:\n ',pmax_stats)
    print('Imax stats:\n ',imax_stats)
    print('Vmax stats:\n ',vmax_stats)
    
    # stats for the whole IV curve
    stats = returnStats(measured=I_points, modeled=i_sim)
    print('IV curve stats:\n',stats)
    
    
    
    
    # saves data into a pandas dataframe and into a file if this is given
    all_data = pd.DataFrame({'voltage':v,'current_exp':I_points,'current_sim':i_sim})
    
    if output_file:
        all_data.to_csv(output_file)
        
        

    # here you can change the styles in the graph - check the top of the module
    sns.set_style({'font.family':'serif', 'font.serif':'Trebuchet MS', 'text.color': '#6B33FF', 'font.size':14})
   
    # also see plotGraphs for more info about my custom made function with lots of options that can be found there

    
    plt.plot(v,i_sim,'b--',label = 'Simulated IV')
    plt.plot(v,I_points,'gp', label = 'Experimental IV')
    plt.xlabel('Voltage(V)')
    plt.ylabel('Current (A)')
    plt.legend()
    plt.show()
            
  



#One thin film CdTe cell
os.chdir('..')
     
f1 = os.path.join(os.getcwd(),'Test Files//FBinput_iv.csv')

#mini _module with six cells

f2 = os.path.join(os.getcwd(),'Test Files//iv_compare_mini_moduleMOB.csv')


# This command runs the module
if __name__=='__main__':
    
    main(f2,cells_in_series=6, columns=['V','I'], output_file=None)
    


