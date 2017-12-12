'''
Created on 13 May 2015

@author: elek2

All my statistical analysis in terms of RMSE, MBE, MAE files go here 

'''
import StatMetrics as sm
import pandas as pd
import numpy as np




def returnStats(measured, modeled):

    rmseabs = sm.absoluteRMSE(modeled,measured) 
    maeabs = sm.absoluteMAE(modeled,measured)
    mbeabs = sm.absoluteMBE(modeled,measured)
    #rmserel = sm.relativeRMSE(modeled,measured)            
    
    rmserel = sm.relativeRMSE(modeled,measured)  
    maerel = sm.relativeMAE(modeled,measured)
    mberel = sm.relativeMBE(modeled,measured)
     
    #===========================================================================
    # print("absolute RMSE is: ",rmseabs) 
    # print("absolute MAE is: ", maeabs) 
    # print("absolute MBE is: ", mbeabs)          
    # print('\n')         
    # print("relative RMSE is: ",rmserel)
    # print("relative MAE is: ", maerel) 
    # print("relative MBE is: ", mberel)
    #===========================================================================
    
    maxdev = (np.max(np.abs(measured-modeled)))
    
    #print("maximum deviation is: ",maxdev)
   
    #print("\n")
    
    dfstats = pd.DataFrame({'RMSE':rmseabs,'MAE':maeabs, 'MBE':mbeabs, 'rRMSE':rmserel, 
                            'rMAE':maerel, 'maxdev':maxdev,'rMBE':mberel}, index =[0])
    
    return dfstats
    
 
    
    


























     