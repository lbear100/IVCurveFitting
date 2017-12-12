'''
This module contains useful functions 
for plotting simple  2D,3D graphs, contours and histograms
of higher quality than excel and exporting them at .tiff or .eps etc

Author: Elena Koumpli
Last updated: 21/07/2017

'''


from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator
from scipy.stats import linregress
from tkinter import filedialog, Tk




rcParams.update({'figure.autolayout': True}) #very important command for sizing automatically
rcParams['mathtext.default'] = 'regular' # this one removes the italics from inserted functions

savefolder = r'C:\Users\elek2_backup\Helena\DATA_WORKSPACES\matplotlib_savefolder'

# the following expression works with the 3.0 interpreter

lbl = 18 #label font size
leg= 16 # legend font size
tck = 16 # ticks font size
res =300 

#filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
# legend locations 
    #===========================================================================
    # center
    # right
    # lower right
    # upper left
    # center left
    # upper center
    # center right
    # best
    # upper right
    # lower center
    # lower left
    #===========================================================================



def plotDynGraphsOnXY(x,y, ysec=[], legends = [], axlabels = [], title='', 
                      scatter= False, marker_line = True, xticklocator = None, 
                      yticklocator = None,lower_limits=None, marker = False,
                      prefs = {'lbl':20,'leg':18, 'tck':18,'res': 300, 
                               'leg_loc':'upper right','fformat':'tiff',
                               'leg_layout':'vertical', 'bg_colour':'white'},
                      regression= None, showlegend = True, show_grid=False, 
                      set_ticks = None, show = True, yscale='linear',
                      autoformat_date = False, save_graph = False):
    
    
    '''
    A useful function which wraps the most essential functionalities of Matplotlib and keeps
    consistent styles for every graph produced, based on personal preferences.
    Draw line and scatter graphs automatically by adding Y arguments. 
    Colours and styles are re-cycled, taken from default lists
    IF scatter is False then a line graph is produced else a scatter graph is produced.
    
    NOTE: Some styles can be overriden/altered by using the python package seaborn, 
    Just type seaborn.set() prior to plotDynGraphsOnXY() to adopt the seaborn styles.
    
    Parameters
    -----------
    
    x        : float, array or list
               x-axis parameter of length M
        
    
    y        : list of floats, arrays or lists
             y-axis parameters. The size of each argument in the list, N, must be of equal length to M
             (M = N)
    
    ysec     : list of floats, arrays or lists, Default empty list
             Secondary y-axis parameters. The size of each argument in the list, N, must be of equal length to M
             (M = N) 
    
    legends  : A list of string(s), Default empty 
    
          
    axlabels : A list of string(s), default Empty, 
                labels for X and Y axes 
    
    
    title    : String,  default Empty
                The title of the graph
    
      
    scatter :  Boolean, default False
                Whether the graph is a line or scatter
    
     
    marker_line: Boolean, default True
                 Whether for the line graph, a marker line is preferred or not
    
    marker      : Boolean, default False
                Whether for the line graph, a marker  is preferred or not
                
    xticklocator :Integer, default None
                 Where the tick labels are located in X axis, 
                 if None then the default style is applied
    
                  
    yticklocator : Integer, default None
                 Where the tick labels are located in primary Y axis, 
                 if None then the default style is applied
        
    regression   : string, default None 
                   if not None, then a trendline is applied to the scatter plot 
                   currently two strings are recognised 'linear' and 'log' for linear and
                   logarithmic fit respectively 
    
                   
    showlegend   : boolean, default False,
                    whether to show legend or not. 
                    Properties of the legend will have to be hard coded 
    
    show_grid    : boolean, default False
                    Whether to show grid (True for both axes) or not (False)
                    
    set_ticks   : A list of strings, default None
                if not None, this replaces the tick labels of the x axis
                    
    lower_limits : A list of floats, size = 2, default None
                if not None it sets the lower limits of the X and (primary) Y axis where
                X limit = lower_limits[0]
                Y limit = lower_limits[1]
    
    prefs      : A dictionary default {'lbl':20,'leg':18, 'tck':18,'res': 300, 'fformat':300 }
                Common font size preferences where 
                      lbl denotes label font size
                      leg denotes legend font size
                      tck denotes tick font size
                      res denotes dpi resolution size for saved graphs
                      fformat is the file format for the produced graph
                      leg_layout horizontal or vertical legend
                      leg_loc location of the legend 
                      bg_colour this sets the background colour of the figure  
             
    show       : Boolean, default True
                 if true then the graph is shown else it is not.
    
    autoformat_date : Boolean, default False
                    if autoformat is True then dates are formatted accordingly
    
    yscale     :string, default 'linear'  
                The scale of the primary Y axis
    
    save_graph : Boolean, default False
                if True a dialogue window opens where it asks for the file to save the graph in 
    
    Returns
    -------

    A graph object which has the option to be saved in a file (default format .tiff)
    
    '''
    
    from itertools import cycle
    
    #default dictionary
    
    pprefs = {'lbl':18,'leg':16, 'tck':16,'res': 300, 
                               'leg_loc':'upper right','fformat':'tiff',
                               'leg_layout':'vertical'}
    
    lw = 1.5 #this looks good
    
    # updates based on personal preferences given as defaults in the arguments
    pprefs.update(prefs)
    
    # colours are cycled as well as styles for both scatter and line graphs
    
    fig, ax = plt.subplots() # you need this when you add second axis
    
    
    # more elements can be added in the lists below for very large number of Y data
    gen_colors = ['black','blue', 'orange','cadetblue','purple','red', 'green','mintcream','steelblue','cyan','lightgreen','orchid','pink', 
                  'yellow', 'tan', 'sandybrown', 'goldenrod', 'yellowgreen'] 
     
    gen_lines = ["-o","--^","--s","--v","-.","--x","--o",":o","-.^","-o"]
    
    gen_lines_no_marker = ['-', '--', '-.', ':']
    
    gen_scatter = ["o", "x","^", "<", ">","1","2","3", "4","8","s","p","P","*","h","H","+",".","D","d"]

    lns = []
        
    
    if (scatter or not marker_line):
        linecycler = cycle(gen_scatter) #cycles through line styles
    
    elif (marker is False and scatter is False):
        
        linecycler = cycle(gen_lines_no_marker)
              
    else:    
        linecycler = cycle(gen_lines)

    colorcycler = cycle(gen_colors)
    
    for i,yi in enumerate(y):
        
        
        if not scatter:    
            
            try:
                ln, = ax.plot(x,yi,  next(linecycler), color = next(colorcycler), 
                           label = legends[i], linewidth = lw)
                
            except IndexError:
                
                ln, = ax.plot(x,yi,  next(linecycler), color = next(colorcycler), 
                           linewidth = lw)
            
            
            lns.append(ln)

        else:
            
            colour = next(colorcycler) # this is to get same colour for both scatter and regression lines
            
            try:
                ln = ax.scatter(x,yi,30, marker = next(linecycler), color = colour, label = legends[i])
            except IndexError:
                ln = ax.scatter(x,yi,30, marker = next(linecycler), color = colour)
            
            lns.append(ln)
            
            if regression is not None:
                
                if regression =='linear':
                
                    a,b,r_value,p_value,std = linregress(x,yi)
                    yi_m = a*x + b
                    ax.plot(x,yi_m,'--',color = colour, linewidth = 1.5)
                    
                    # NOTE: text requires x and y coordinates as the first two arguments 
                    plt.figtext(0.55,0.8,'R$^2$ = '+str(r_value.round(2)),fontsize=18,backgroundcolor='w')
                    
                elif regression =='log':
                
                    #params, cov = curve_fit(lambda t,a,b: a+b*np.log(t),  x, yi)
                
                    yi_m = -0.179*np.log(x)+1.3535 #values are an example here
                    ax.plot(x,yi_m,'--',color = colour, linewidth = 1.5)
                    
    
    if ysec:   #checks for empty list

        
        ax2 = ax.twinx()
        
        for k,yk in enumerate(ysec): 
            
            try:
                ln, = ax2.plot(x,yk,  next(linecycler), color = next(colorcycler), 
                           label = legends[len(y)+k], linewidth = lw)
            except IndexError:
                
                ln, = ax2.plot(x,yk,  next(linecycler), color = next(colorcycler), 
                            linewidth = lw)
                
                
            ax2.set_ylabel(axlabels[2], fontsize = pprefs['lbl'])  
            ax2.tick_params(labelsize = pprefs['tck'])
            lns.append(ln)
           
    
    try:
        
        ax.set_xlabel(axlabels[0], fontsize = pprefs['lbl'])
        ax.set_ylabel(axlabels[1], fontsize = pprefs['lbl'])
    
    except IndexError:
        
        ax.set_xlabel('', fontsize = pprefs['lbl'])
        ax.set_ylabel('', fontsize = pprefs['lbl'])
             
    
    if autoformat_date:
        
        plt.gcf().autofmt_xdate()
    
    if showlegend:
        
        leg_loc = pprefs['leg_loc']     
        labs = [l.get_label() for l in lns]
        
        if pprefs['leg_layout']=='horizontal':
            
            plt.legend(lns, labs, loc=leg_loc,prop={'size':pprefs['leg']},ncol = len(y)+len(ysec))
        
        else:
            plt.legend(lns, labs, loc=leg_loc,prop={'size':pprefs['leg']}, facecolor = 'w')
   
    if title:
        
        plt.title(str(title),size = pprefs['tck'])
    
    
    ax.tick_params(axis='both', labelsize = pprefs['tck'])
    
    if show_grid:
        ax.grid(axis = 'both')#,color='AliceBlue') #color='grey')
    
    
    
    if lower_limits:
        
        plt.xlim(lower_limits[0])
        plt.ylim(lower_limits[1])
    
    
   
    
    if xticklocator is not None:
    
        majorLocator = MultipleLocator(xticklocator) 
        ax.xaxis.set_major_locator(majorLocator)
        
    
    
    if yticklocator is not None:
    
        majorLocator = MultipleLocator(yticklocator) 
        ax.yaxis.set_major_locator(majorLocator)
        
    
    if set_ticks is not None:
         
        if xticklocator: #set the frequency for the set_xticks according to the xticklocator
            
            x0 =  [x[i] for i in range(0,len(set_ticks),xticklocator)]
            set_ticks0 = [set_ticks[i] for i in range(0,len(set_ticks),xticklocator)]
            x = x0
            plt.xticks(x,set_ticks0)
        else:
            plt.xticks(x,set_ticks)

    
    plt.yscale(yscale)
    
    
    
    #more options here
    
    #plt.yticks([-1,0,1])
    
    #sets axes to always cros at given x,y
    #plt.axhline(y=0, color='k')
    #plt.axvline(x=0, color='k')
    
    #plt.fill_between(x,y[0],y[1],hatch="\\", facecolor='green') #fill
    #plt.text(0.6,8.7,'Shaded area', fontdict={'family':'Arial','size':lbl}) #write a text

    # adds background colour 
    # check https://www.w3schools.com/cssref/css_colors.asp
    
    ax.set_facecolor(pprefs['bg_colour'])#('lavender') # background color in the graph
    
    
    
    if save_graph: 
        
        from tkinter import filedialog
        
        fformat = pprefs['fformat']
        
        file = filedialog.asksaveasfilename(defaultextension=fformat)
        
        filename = file#+'.'+fformat
        
        plt.savefig(filename,dpi = pprefs['res'],figsize=(10,6))
        
    if show:
        
        plt.show() 
    

