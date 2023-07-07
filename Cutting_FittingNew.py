#############################################################################################################
# @Giulia Faletti                                                                                           #
# Cutting and Fitting algorithm for the luminosity evolution model considering ATLAS redifined (LUCID) data #
#############################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import LoadData as ld
from lmfit import Model
import scipy.integrate as integrate

#Font setting
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica",
  "font.size": 12
})

def Cut_Fit(year, text):
    """Function that performs the necessary cut on the current fill

    Args:
        year (int): current year
        text (str): current fill

    Returns:
        L_fit: cutted data
        T_fit_real: times in second for the cutted data fit
        Y: Luminosity evolution form the fit
        a: fitting parameter
        b: fitting parameter
        c: fitting parameter
        d: fitting parameter
        chi: reduced chi square of the fit
        L_evol: raw luminosity data
        Times: raw Unix time
    """
    Times, L_evol= ld.AtlasData(text, year, grl=True)

    #deleting the null values of the luminosity
    zero=np.where(L_evol<100)
    L_zero=np.delete(L_evol, zero)
    T_zero=np.delete(Times, zero)
        
    #check for enough points
    if len(L_zero)<10:
        zero=np.where(L_evol<5)
        L_zero=np.delete(L_evol, zero)
        T_zero=np.delete(Times, zero)

    #defining the derivative 
    dy = np.zeros(L_zero.shape)
    dy[0:-1] = np.diff(L_zero)/np.diff(T_zero)


    #start to slim down the fit interval       
    L_tofit=[]
    T_tofit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1.5:
            L_tofit.append(L_zero[idx])
            T_tofit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1.5:
            continue     
        
    #evaluating the differences between two subsequent points
    diff=np.diff(L_tofit)
        
    #deleting the discrepancies
    thr=np.max(abs(diff))*0.05
    idx_diff= np.where(abs(diff)>thr)[0]+1
        
    #new slim down of data
    L_tofit2=np.delete(L_tofit, idx_diff)
    T_tofit2=np.delete(T_tofit, idx_diff)
        
    #check for enough points
    if len(L_tofit2) < 30:
        L_tofit2=L_tofit
        T_tofit2=T_tofit
        
    L_fit=L_tofit2
    T_fit=T_tofit2     

    L_fit=np.array(L_fit)
    T_fit=np.array(T_fit) 

    #normalization of the fit interval    
    norm_T_fit=[]
    norm_T_fit=np.array(norm_T_fit)
    for element in T_fit:
        z=(element-np.amin(T_fit))/(np.amax(T_fit)-np.amin(T_fit))
        norm_T_fit=np.append(norm_T_fit, z)


    #defining the fit function
    def fit(x, a, b, c, d):
        return (a*np.exp((-b)*x))+(c*np.exp((-d)*x))

    model=Model(fit)      

    #performing fit 
    if year==16:
        model.set_param_hint('b', value=0.2, min=0, max=80)
        model.set_param_hint('d', value=0.2, min=0, max=80)
        model.set_param_hint('a', value=15, min=1, max=11500)
        model.set_param_hint('c', value=15, min=1, max=11500)
        fit_result=model.fit(L_fit, x=norm_T_fit, a=10, b=0.2, c=10, d=0.2)
    elif year==17:
        model.set_param_hint('b', value=0.2, min=0, max=70)
        model.set_param_hint('d', value=0.2, min=0, max=70)
        model.set_param_hint('a', value=1000, min=1, max=19500)
        model.set_param_hint('c', value=1000, min=1, max=19500)
        fit_result=model.fit(L_fit, x=norm_T_fit, a=1000, b=0.2, c=1000, d=0.2)
    elif year==18:
        model.set_param_hint('b', value=0.2, min=0, max=30)
        model.set_param_hint('d', value=0.2, min=0, max=30)
        model.set_param_hint('a', value=10, min=1, max=16500)
        model.set_param_hint('c', value=10, min=1, max=16500)
        fit_result=model.fit(L_fit, x=norm_T_fit, a=1000, b=0.2, c=1000, d=0.2)

    #transforming the times from unix in seconds
    T_fit_real=T_fit    
    Y=fit(T_fit_real, fit_result.params['a'].value, (fit_result.params['b'].value)/(np.amax(T_fit)-np.amin(T_fit)), fit_result.params['c'].value, fit_result.params['d'].value/(np.amax(T_fit)-np.amin(T_fit)))
    
    return  L_fit, T_fit_real,Y, fit_result.params['a'].value, (fit_result.params['b'].value/(np.amax(T_fit)-np.amin(T_fit))), fit_result.params['c'].value, fit_result.params['d'].value/(np.amax(T_fit)-np.amin(T_fit)), fit_result.redchi, L_evol, Times

def Partial_Cut_Fit(year, text, h):
    """Function that performs the necessary cut on the current fill

    Args:
        year (int): current year
        text (str): current fill
        h (int): number of hours for the evaluation of the fitting parameters

    Returns:
        a: fitting parameter
        b: fitting parameter
        c: fitting parameter
        d: fitting parameter
    """
    Times, L_evol= ld.AtlasData(text, year, grl=True)
    
    
    #deleting the null values of the luminosity
    zero=np.where(L_evol<100)
    L_zero=np.delete(L_evol, zero)
    T_zero=np.delete(Times, zero)
        
    #check for enough points
    if len(L_zero)<10:
        zero=np.where(L_evol<5)
        L_zero=np.delete(L_evol, zero)
        T_zero=np.delete(Times, zero)

    #defining the derivative 
    dy = np.zeros(L_zero.shape)
    dy[0:-1] = np.diff(L_zero)/np.diff(T_zero)

     #start to slim down the fit interval       
    L_tofit=[]
    T_tofit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1.5:
            L_tofit.append(L_zero[idx])
            T_tofit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1.5:
            continue     
        
    #evaluating the differences between two subsequent points
    diff=np.diff(L_tofit)
        
    #deleting the discrepancies
    thr=np.max(abs(diff))*0.05
    idx_diff= np.where(abs(diff)>thr)[0]+1
        
    #new slim down of data
    L_tofit2=np.delete(L_tofit, idx_diff)
    T_tofit2=np.delete(T_tofit, idx_diff)
        
    #check for enough points
    if len(L_tofit2) < 30:
        L_tofit2=L_tofit
        T_tofit2=T_tofit
        
    L_fit=L_tofit2
    T_fit=T_tofit2   
    
    L_fit=np.array(L_fit)
    T_fit=np.array(T_fit)
      
    #selecting only the designed hour for the partial knowledge
    if h==1:
        indi=np.where(T_fit<=(np.amin(T_fit)+3600))
        if len(T_fit[indi])<5:
            indiz=np.where(T_fit<=(np.amin(T_fit)+7200))
            L_fit=L_fit[indiz]
            T_fit=T_fit[indiz]
            
        elif len(T_fit[indi])>=5:
            L_fit=L_fit[indi]
            T_fit=T_fit[indi]
    elif h!=1:
        indi=np.where(T_fit<=(np.amin(T_fit)+(h*3600)))
        L_fit=L_fit[indi]
        T_fit=T_fit[indi]  
      
    #normalization of the fit interval    
    norm_T_fit=[]
    norm_T_fit=np.array(norm_T_fit)
    for element in T_fit:
        z=(element-np.amin(T_fit))/(np.amax(T_fit)-np.amin(T_fit))
        norm_T_fit=np.append(norm_T_fit, z)
    
      
    #defining the fit function
    def fit(x, a, b, c, d):
        return (a*np.exp((-b)*x))+(c*np.exp((-d)*x))

    model=Model(fit)    

    #performing fit 
    if year==16:
        model.set_param_hint('b', value=0.2, min=0, max=50)
        model.set_param_hint('d', value=0.2, min=0, max=50)
        model.set_param_hint('a', value=10, min=1, max=8500)
        model.set_param_hint('c', value=10, min=1, max=8500)
        fit_result=model.fit(L_fit, x=norm_T_fit, a=10, b=0.2, c=10, d=0.2)
    elif year==17:
        model.set_param_hint('b', value=0.2, min=0, max=50)
        model.set_param_hint('d', value=0.2, min=0, max=50)
        model.set_param_hint('a', value=1000, min=1, max=19000)
        model.set_param_hint('c', value=1000, min=1, max=19000)
        fit_result=model.fit(L_fit, x=norm_T_fit, a=1000, b=0.2, c=1000, d=0.2)
    elif year==18:
        model.set_param_hint('b', value=0.2, min=0, max=50)
        model.set_param_hint('d', value=0.2, min=0, max=50)
        model.set_param_hint('a', value=10, min=1, max=18500)
        model.set_param_hint('c', value=10, min=1, max=18500)
        fit_result=model.fit(L_fit, x=norm_T_fit, a=10, b=0.2, c=10, d=0.2)

    return fit_result.params['a'].value, (fit_result.params['b'].value/(np.amax(T_fit)-np.amin(T_fit))), fit_result.params['c'].value, fit_result.params['d'].value/(np.amax(T_fit)-np.amin(T_fit))
  
def PlotLumiEvol(L_fit, T, L, a,b,c,d, text, chi, year, L_ev, Tim):
    """Function that plots all the needed graphs.

    Args:
        L_fit: cutted data
        T_fit_real: times in second for the cutted data fit
        Y: Luminosity evolution form the fit
        a: fitting parameter
        b: fitting parameter
        c: fitting parameter
        d: fitting parameter
        chi: reduced chi square of the fit
        L_evol: raw luminosity data
        Times: raw Unix time
    """
    fig, ax= plt.subplots()
    ax.plot(T/3600, L_fit, "b.", markersize=4)
    ax.plot(T/3600, L, 'r-', label='Fit')
    ax.plot([], [], 'kx ', label=r'$\tilde{\chi}^2$='+'{:.2f}'.format(chi))
    ax.set_title('{}'.format(text))
    ax.set_xlabel('Times [h]')
    ax.set_ylabel('Luminosity evolution [$\mathrm{Hz}/\mathrm{\mu b}$]')
    plt.legend(loc='best')
    plt.savefig('Cutting_FittingNew/20{}_Graphs/{}_LuminosityEvolutionCutFit.pdf'.format(str(year),text))
    
    fig1, ax1= plt.subplots()
    ax1.plot(Tim, L_ev, "k-")
    ax1.set_title('{}'.format(text))
    ax1.set_xlabel('Times [s]')
    ax1.set_ylabel('Luminosity evolution [$\mathrm{Hz}/\mathrm{\mu b}$]')
    plt.savefig('Cutting_FittingNew/20{}_Graphs/init/{}_LuminosityEvolution.pdf'.format(str(year),text))
    
    
    fig2, ax2= plt.subplots()
    ax2.plot(T/3600, L_fit, "b.", markersize=4)
    ax2.plot(T/3600, L, 'r-', label='Fit')
    ax2.plot([], [], 'kx ', label=r'$\tilde{\chi}^2$='+'{:.2f}'.format(chi))
    ax2.set_yscale('log')
    ax2.set_title('{}'.format(text))
    ax2.set_xlabel('Times [h]')
    ax2.set_ylabel('Luminosity evolution [$\mathrm{Hz}/\mathrm{\mu b}$]')
    plt.legend(loc='best')
    plt.savefig('Cutting_FittingNew/20{}_Graphs/log/{}_LuminosityEvolutionCutFit_log.pdf'.format(str(year),text))



#selecting current year
year=18

#defining the knowledge level
PartialKnowledge=False
h=15 #number of hours for the evaluation of the fitting parameters

#Plotting Luminosity Evolution Graphs
plot=True

#loading fill number
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()
FillNumber17=np.delete(FillNumber17, np.where(FillNumber17==6160)[0])



#load turnaround times and fill times 
data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 = ld.loadFill()
data_ta16_sec = data_ta16*3600 
data_tf16_sec = data_tf16*3600  
data_ta17_sec = data_ta17*3600 
data_tf17_sec = data_tf17*3600
data_ta18_sec = data_ta18*3600 
data_tf18_sec = data_tf18*3600

if year==16:
    FillNumber=FillNumber16
    ta=data_ta16_sec
    tf=data_tf16_sec
elif year==17:
    FillNumber=FillNumber17
    FillNumber_Prev=FillNumber16
    previous_year=16
    ta=data_ta17_sec
    tf=data_tf17_sec
elif year==18:
    FillNumber=FillNumber18
    FillNumber_Prev=FillNumber17
    previous_year=17
    ta=data_ta18_sec
    tf=data_tf18_sec
     

if PartialKnowledge==True:
    with open('Cutting_FittingNew/FitCoefficients{}_{}h.txt'.format(str(year), str(h)), 'w') as f:
        f.write('')
        f.close()    
        
elif PartialKnowledge==False:
    with open('Cutting_FittingNew/FitCoefficients{}.txt'.format(str(year)), 'w') as f:
        f.write('')
        f.close()     
    
    for i in range(len(FillNumber)):
        text = str(int(FillNumber[i])) #number of current fill
        with open('Cutting_FittingNew/20{}/{}.txt'.format(str(year),text), 'w') as f:
            f.write('')
            f.close()  

    
for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill
    
    if PartialKnowledge==True:
        
        #performing the cutting and fitting algorithm
        A,B,C,D = Partial_Cut_Fit(year, text, h)
    
        #saving the results
        with open('Cutting_FittingNew/FitCoefficients{}_{}h.txt'.format(str(year), str(h)), 'a') as f:
            f.write(text)
            f.write(' ')
            f.write(str(A))
            f.write(' ')
            f.write(str(B))
            f.write(' ')
            f.write(str(C))
            f.write(' ')
            f.write(str(D))
            f.write('\n')
            
    elif PartialKnowledge==False:
        
        #performing the cutting and fitting algorithm
        L_fit, T, L, a,b,c,d, chi, L_ev, Tim= Cut_Fit(year, text)
    
        if plot==True:
            PlotLumiEvol(L_fit, T, L, a,b,c,d, text, chi, year, L_ev, Tim)
        
        #saving results
        for x in range(len(T)):
            with open('Cutting_FittingNew/20{}/{}.txt'.format(str(year),text), 'a') as f:
                f.write(str(T[x]))
                f.write(' ')
                f.write(str(L_fit[x]))
                f.write('\n')
                
        with open('Cutting_FittingNew/FitCoefficients{}.txt'.format(str(year)), 'a') as f:
            f.write(text)
            f.write(' ')
            f.write(str(a))
            f.write(' ')
            f.write(str(b))
            f.write(' ')
            f.write(str(c))
            f.write(' ')
            f.write(str(d))
            f.write('\n')
            