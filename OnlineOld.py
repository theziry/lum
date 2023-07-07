#########################################################
#      @Giulia Faletti                                  #
# Online optimisation strategy for analysis of 2017     #
#########################################################
from scipy.optimize import least_squares
import scipy.integrate as integrate
import scipy.interpolate
import math
from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt
import LoadData as ld
from MPL import *
#import ConfrontoFitMPL

#font settings
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica",
  "font.size": 12
})
#Selecting Current Year
year=17


#loading fill number
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()

#load turnaround times and fill times - cheking
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

#loading fill evolution fit coefficients  
f=open('FitCoefficients/a{}.txt'.format(str(year)),"r")
lines=f.readlines()
a=[]
for x in lines:
    a.append(float(x.split(' ')[0]))

f=open('FitCoefficients/b{}.txt'.format(str(year)),"r")
lines=f.readlines()
b=[]
for x in lines:
    b.append(float(x.split(' ')[0]))
    
f=open('FitCoefficients/c{}.txt'.format(str(year)),"r")
lines=f.readlines()
c=[]
for x in lines:
    c.append(float(x.split(' ')[0]))
    
f=open('FitCoefficients/d{}.txt'.format(str(year)),"r")
lines=f.readlines()
d=[]
for x in lines:
    d.append(float(x.split(' ')[0]))
    
a=np.array(a)
b=np.array(b) 
c=np.array(c)
d=np.array(d) 

#defining the fit function
def fit(x, aa, bb, cc, dd):
    return (aa*np.exp((-bb)*x))+(cc*np.exp((-dd)*x))

#loading numerical results - for checking
f=open('NumericalOptimization/res_opt_20{}.txt'.format(str(year)),"r")
lines=f.readlines()
t_opt_num=[]
for x in lines:
    t_opt_num.append(float(x.split(' ')[0]))
    
t_opt_num=np.array(t_opt_num)

#cleaning file for saving the results
with open('OnlineStrategyComplete/FutureFill{}.txt'.format(str(year)), 'w') as f:
    f.write('')
    f.close() 
with open('OnlineStrategyComplete/OptimalFill{}.txt'.format(str(year)), 'w') as f:
    f.write('')
    f.close() 



#defining the mostprobable sampling
sampling=np.arange(0, 30*3600, 1)


#defining the phyisics time
Constraint=np.sum(ta)+np.sum(tf)   

t_future=[]
t_optimal=[]
for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill

    #Current Fill Peak Luminosity determination
    L_peak=PeakLumi(year, text)

    #MostProbableLuminosity evaluation
    Sampling(year, L_peak)


    #Checking the model of most probable luminosity
    #fig1, ax1=plt.subplots()
    #ax1.set_xlabel('Times [s]')
    #####ax1.set_ylabel('Luminosity evolution [$\mathrm{Hz}/\mathrm{\mu b}$]')
    ####ax1.set_title('2016')
    ###for ione in range(len(FillNumber_Prev)):
        ##text2 = str(int(FillNumber_Prev[ione])) #number of current fill
        #x, y = ConfrontoFitMPL.Checking_MPL_16(text2)
        #ax1.plot(x, y, '--', linewidth=2, alpha=0.6)
    
    #for ionz in range(len(FillNumber[:i+1])):
        #text3=str(int(FillNumber[ionz])) #number of current fill
        #xx, yy = ConfrontoFitMPL.Checking_MPL_17(text3)
        #ax1.plot(xx, yy, '--', color='b', linewidth=2, alpha=0.6)
 
    
    #check the constraint
    test=np.sum(ta[:i+1])+ sum(t_optimal)
    text1= str(int(FillNumber[i-1]))
    if test>=Constraint:
    #if math.isclose(test,Constraint, rel_tol=0.015):
        print("________________________________________________")
        print("| End of Available Time!")
        print("| The last fill is the fill numered:", text1)
        print("| ", i, "/", len(FillNumber), "fills have been performed.")
        print("________________________________________________")
        break
    
    
    #defining average of the turnaround times
    tau=(np.sum(ta[:i+1]))/(i+1)

    T_fit_real, Y=OnlineSampling(year, text, L_peak)
    #ax1.plot(T_fit_real, Y, linewidth=3, color='k', label="Current Fill Lumi Evol")
    T_mp, L_mp=Mode(year)

    #print(T_mp, '\n', L_mp)
    """
    #Most probable luminosity interpolation
    fL_mp=scipy.interpolate.interp1d(T_mp, L_mp)

    
    #defining the first equation of the optimization
    f1=lambda x1: fL_mp(x1)-((1/(tau+x1))*(integrate.quad(fL_mp, 0, x1)[0]))
    
    #solving the first equation
    x0=15*3600
    res1=least_squares(f1, x0, bounds=(1800, 28*3600))
    t_future.append(res1.x[0])
    

    x2=np.linspace(0, res1.x[0], 150)
    y2=fL_mp(x2)
    #ax1.plot(x2, y2, linewidth=3, color='r', label='Current Fill Most Probable Luminosity')
    #plt.legend()
    #plt.savefig("OnlineStrategyComplete/CheckMPL/CheckingMPL_{}.pdf".format(text))
    
    #save future fill optimal time
    with open('OnlineStrategyComplete/FutureFill{}.txt'.format(year), 'a') as f:
        f.write(str(res1.x[0]))
        f.write('\n') 
        
    #defining the second equation of the optimization
    f2=lambda x2: (a[i]*np.exp(-b[i]*x2)+c[i]*np.exp(-d[i]*x2))-((1/(tau+res1.x[0]))*(integrate.quad(fL_mp, 0, res1.x[0])[0]))
    
    #solving the second equation
    x00=data_tf17_sec[i]
    res2=least_squares(f2, x00, bounds=(1800, 28*3600))
    t_optimal.append(res2.x[0])
    
    #save current fill optimal time
    with open('OnlineStrategyComplete/OptimalFill{}.txt'.format(year), 'a') as f:
        f.write(str(res2.x[0]))
        f.write('\n')   
        
    #print(x00, res2.x)


t_optimal=np.array(t_optimal)

tta_opt=sum(ta[:len(t_optimal)])

print('______________________________________________________________________________________________________')
print("|____________________________________________Check the Constraint____________________________________|")
print("| T_real=", (np.sum(tf)/3600)+(np.sum(ta)/3600), "[h]", "T_optimal=", (np.sum(t_optimal)/3600)+(tta_opt/3600), "[h]", "T_num=", (np.sum(t_opt_num)/3600)+(np.sum(ta)/3600), "[h]")
print("| T_real=", np.sum(tf)/(3600*24)+(np.sum(ta)/(3600*24)), "[d]", "T_optimal=", (np.sum(t_optimal)/(3600*24))+(tta_opt/(3600*24)), "[d]", "T_num=", (np.sum(t_opt_num)/(24*3600))+(np.sum(ta)/(3600*24)), "[d]")
print('|_____________________________________________________________________________________________________')

if len(t_optimal)<len(FillNumber):
    FillNumber=FillNumber[:len(t_optimal)]

fig3, ax3= plt.subplots()
ax3.plot(FillNumber, t_optimal/3600, "b.", markersize=4)
ax3.axhline(1800/3600, color="r", linestyle='-')
ax3.axhline(28, color="r", linestyle='-')
ax3.set_title('Optimization Boundaries')
ax3.set_xlabel('Fill Number')
ax3.set_ylabel('Optimized Times [h]')
plt.savefig('OnlineStrategyComplete/OptimizationBoundaries{}.pdf'.format(year))
plt.show()
"""