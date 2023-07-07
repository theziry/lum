##############################################################################
#      @Giulia Faletti                                                       #
# Online optimisation strategy for analysis of Run2 with data from LUCID     #
##############################################################################

from scipy.optimize import least_squares
import scipy.integrate as integrate
import scipy.interpolate
import math
from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt
import LoadData as ld
from MPL import *
#import time as t

#start=t.time()

#font settings
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica",
  "font.size": 12
})

#Selecting Current Year 17 or 18
year=18

#plotting 
plot=True

#defining the knowledge level
PartialKnowledge=False
h=15 #number of hours for the evaluation of the fitting parameters

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

#parameters setting
if year==16:
    FillNumber=FillNumber16
    ta=data_ta16_sec
    tf=data_tf16_sec
elif year==17:
    FillNumber=np.delete(FillNumber17, np.where(FillNumber17==6160)[0])
    FillNumber_Prev=FillNumber16
    previous_year=16
    ta=np.delete(data_ta17_sec, np.where(FillNumber17==6160)[0])
    tf=np.delete(data_tf17_sec, np.where(FillNumber17==6160)[0])
elif year==18:
    FillNumber=FillNumber18
    FillNumber_Prev=FillNumber17
    previous_year=17
    ta=data_ta18_sec
    tf=data_tf18_sec


#loading fill evolution fit coefficients  
if PartialKnowledge==False:
    f=open('Cutting_FittingNew/FitCoefficients{}.txt'.format(str(year)),"r")
    lines=f.readlines()
    a=[]
    b=[]
    c=[]
    d=[]
    for x in lines:
        a.append(float(x.split(' ')[1]))
        b.append(float(x.split(' ')[2]))
        c.append(float(x.split(' ')[3]))
        d.append(float(x.split(' ')[4]))  
    
elif PartialKnowledge==True:
    f=open('Cutting_FittingNew/FitCoefficients{}_{}h.txt'.format(str(year), str(h)),"r")
    lines=f.readlines()
    a=[]
    b=[]
    c=[]
    d=[]
    for x in lines:
        a.append(float(x.split(' ')[1]))
        b.append(float(x.split(' ')[2]))
        c.append(float(x.split(' ')[3]))
        d.append(float(x.split(' ')[4]))  
 
f.close()   
a=np.array(a)
b=np.array(b) 
c=np.array(c)
d=np.array(d) 

#defining the fit function
def fit(x, aa, bb, cc, dd):
    return (aa*np.exp((-bb)*x))+(cc*np.exp((-dd)*x))

#loading numerical results - for checking
f=open('NumericalOptimizationNew/res_opt_20{}.txt'.format(str(year)),"r")
lines=f.readlines()
t_opt_num=[]
for x in lines:
    t_opt_num.append(float(x.split(' ')[0]))
    
t_opt_num=np.array(t_opt_num)

#cleaning file for saving the results
with open('OnlineNew/FutureFill{}.txt'.format(str(year)), 'w') as f:
    f.write('')
    f.close() 
with open('OnlineNew/OptimalFill{}.txt'.format(str(year)), 'w') as f:
    f.write('')
    f.close() 


#defining the phyisics time
Constraint=np.sum(ta)+np.sum(tf)   

#performing the online optimisation
t_future=[]
t_optimal=[]
for i in range(len(FillNumber)): #len(FillNumber)
    fill=FillNumber[i]
    text = str(int(fill)) #number of current fill
    
    #check the constraint
    test=np.sum(ta[:i+1])+ sum(t_optimal)
    text1= str(int(FillNumber[i-1]))
    if test>=Constraint:
        print("________________________________________________")
        print("| End of Available Time!")
        print("| The last fill is the fill numered:", text1)
        print("| ", i, "/", len(FillNumber), "fills have been performed.")
        print("________________________________________________")
        break
    
    #Current Fill Peak Luminosity determination
    L_peak=PeakLumiNew(year, text)

    #MostProbableLuminosity evaluation
    #Tmp, Lmp=NormMPLNew(L_peak, year, fill)
    Tmp, Lmp=MPLNew(year, fill)
    #Tmp2, Lmp2=MPL(year, fill)
    
    
    if plot==True:
        #mostProbableLuminosity Plots
        fig, ax=plt.subplots()
        ax.plot(Tmp, Lmp, "b.")
        ax.set_xlabel("Tmp")
        ax.set_ylabel("Lmp")
        ax.set_title("{}".format(text))
        #plt.show()
        plt.savefig("MPL/Graphs{}/{}.pdf".format(year, text))
    
        #fig, ax=plt.subplots()
        #ax.plot(Tmp, Lmp, "b.")
        #ax.set_xlabel("Tmp")
        #ax.set_ylabel("Lmp")
        #ax.set_title("{} - without normalisation".format(text))
        #plt.savefig("MPL/Graphs{}/{}withoutNormalisation.pdf".format(year, text))

    #defining average of the turnaround times
    tau=(np.sum(ta[:i+1]))/(i+1)
     
    #Most probable luminosity interpolation
    fLmp=scipy.interpolate.interp1d(Tmp, Lmp)
    F=lambda x : fLmp(x) 
    
    if plot==True:
        #plotting the interpolation
        fig, ax=plt.subplots()
        ax.plot(Tmp, Lmp, "b.")
        ax.plot(Tmp, fLmp(Tmp), "r--")
        ax.set_xlabel("Tmp")
        ax.set_ylabel("Lmp")
        ax.set_title("{}".format(text))
        plt.savefig("MPL/Graphs{}/Interp_{}.pdf".format(year, text))
        #plt.show()
    
    T_mp=np.array(Tmp)
    L_mp=np.array(Lmp)

    #defining the first equation of the optimization
    f1=lambda x1: fLmp(x1)-((1/(tau+x1))*(integrate.quad(F, 0, x1)[0]))
    
    #solving the first equation
    x0=15*3600
    res1=least_squares(f1, x0, bounds=(1800, 28*3600))
    t_future.append(res1.x[0])
    
    #save future fill optimal time
    with open('OnlineNew/FutureFill{}.txt'.format(year), 'a') as f:
        f.write(str(res1.x[0]))
        f.write('\n') 
   
    #defining the second equation of the optimization
    f2=lambda x2: (a[i]*np.exp(-b[i]*x2)+c[i]*np.exp(-d[i]*x2))-((1/(tau+res1.x[0]))*(integrate.quad(F, 0, res1.x[0])[0]))
    
    #solving the second equation
    x00=tf[i]
    res2=least_squares(f2, x00, bounds=(1800, 28*3600))
    t_optimal.append(res2.x[0])
    
    #save current fill optimal time
    with open('OnlineNew/OptimalFill{}.txt'.format(year), 'a') as f:
        f.write(str(res2.x[0]))
        f.write('\n') 
    
t_optimal=np.array(t_optimal)
tta_opt=sum(ta[:len(t_optimal)])
#print(t_optimal/3600)

#Check the constraint in hours and days
print('__________________________________________________________________________20{}___________________________________________________________________'.format(year))
print("|______________________________________________________________Check the Constraint_________________________________________________________|")
print("| T_real=", (np.sum(tf)/3600)+(np.sum(ta)/3600), "[h]", "T_optimal=", (np.sum(t_optimal)/3600)+(tta_opt/3600), "[h]", "T_num=", (np.sum(t_opt_num)/3600)+(np.sum(ta)/3600), "[h]")
print("| T_real=", np.sum(tf)/(3600*24)+(np.sum(ta)/(3600*24)), "[d]", "T_optimal=", (np.sum(t_optimal)/(3600*24))+(tta_opt/(3600*24)), "[d]", "T_num=", (np.sum(t_opt_num)/(24*3600))+(np.sum(ta)/(3600*24)), "[d]")
print('|____________________________________________________________________________________________________________________________________________')

#Plotting optimisation boundaries
if len(t_optimal)<len(FillNumber):
    FillNumber=FillNumber[:len(t_optimal)]

fig3, ax3= plt.subplots()
ax3.plot(FillNumber, t_optimal/3600, "b.", markersize=4)
ax3.axhline(1800/3600, color="r", linestyle='-')
ax3.axhline(28, color="r", linestyle='-')
ax3.set_title('20{} - Optimization Boundaries'.format(year))
ax3.set_xlabel('Fill Number')
ax3.set_ylabel('Optimized Times [h]')
plt.savefig('OnlineNew/OptimizationBoundaries{}.pdf'.format(year))
#plt.show()


#stop=t.time()
#print('time=', stop-start)
