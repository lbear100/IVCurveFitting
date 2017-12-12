'''
Created on 13 May 2015

This module has the most common statistical metrics
Results agree with the Python StatsModels.tools.eval_measures


'''
import numpy as np

def absoluteRMSE(x_model,y_meas):
    
    res = (x_model-y_meas)**2
    sum1 = np.sum(res)
    size = np.size(x_model)
    result = np.sqrt(sum1/size) 
    return result
    
def absoluteMAE(x_model,y_meas):
    
    res = np.abs(x_model-y_meas)
    sum1 = np.sum(res)
    size = np.size(x_model)
    result = (sum1/size) 
    return result    

def absoluteMBE(x_model,y_meas):
    
    result = np.average((x_model-y_meas))
    
    return result    

def relativeRMSE(x_model,y_meas):
    
    num = np.sqrt(np.average((x_model-y_meas)**2))
    den  = np.average(y_meas)
    
    result = num/den *100
   
    
    return result



def relativeMAE(x_model,y_meas):
     
    num = np.sum(np.abs(x_model-y_meas))
    den  = np.sum(y_meas)
    
    result = num/den *100
    
    return result
       

def relativeMBE(x_model,y_meas):
    
    num = np.sum(x_model-y_meas)
    
    den  = np.sum(y_meas)
    
    result = num/den *100
    
    return result
    


def MAD(data):
    
    '''
    Data: is a numpy array
    
    Returns: 
    
    Median absolute deviation
    
    '''
    med = np.median(data)
    #mean absolute deviation
    mad = np.median(np.abs(data-med))
    
    return mad

#===============================================================================
# EXAMPLE
#
# x = np.array([1,2,3,4,5])
# y = np.array([0.9,2.1,3.4,4.3,5.2])
#  
# z = relativeMBE(y,x)
#  
# print(z)
#===============================================================================









