import numpy as np
import matplotlib.pyplot as plt
import LoadData as ld
from matplotlib.ticker import MaxNLocator
import scipy.integrate as integrate

def set_axis_style(ax, labels):
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1), labels)
    #ax.set_xlim(0.25, len(labels) + 0.75)
    #ax.set_xlabel('Sample name')

plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica",
  "font.size": 16
})

#current year
year=18

#Load Data
#loading fill number
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()

#load turnaround times and fill times --- REAL
data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 = ld.loadFill()
data_ta16_sec = data_ta16*3600 
data_tf16_sec = data_tf16*3600  
data_ta17_sec = data_ta17*3600 
data_tf17_sec = data_tf17*3600
data_ta18_sec = data_ta18*3600 
data_tf18_sec = data_tf18*3600

Lmes16, Lmes17, Lmes18=ld.MeasuredLuminosity() #femtobarn^-1

#current year parameters setting
if year==16:
    FillNumber=FillNumber16
    ta=data_ta16_sec
    tf=data_tf16_sec
    Lmes=Lmes16
elif year==17:
    FillNumber=FillNumber17
    FillNumber_Prev=FillNumber16
    previous_year=16
    ta=data_ta17_sec
    tf=data_tf17_sec
    Lmes=Lmes17
elif year==18:
    FillNumber=FillNumber18
    FillNumber_Prev=FillNumber17
    previous_year=17
    ta=data_ta18_sec
    tf=data_tf18_sec
    Lmes=Lmes18

#loading fill evolution fit coefficients  
f=open('Cutting_Fitting/FitCoefficients{}.txt'.format(str(year)),"r")
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

#loading data --- ONLINE
f=open('Online/OptimalFill{}.txt'.format(year),"r")
lines=f.readlines()
t_online=[]
for x in lines:
    t_online.append(float(x.split(' ')[0]))

t_online=np.array(t_online)  

#loading data --- NUMERICAL
f=open('NumericalOptimization/res_opt_20{}.txt'.format(year),"r")
lines=f.readlines()
t_num=[]
for x in lines:
    t_num.append(float(x.split(' ')[0]))
    
t_num=np.array(t_num)

#optimal turn around time
tta_opt=sum(ta[:len(t_online)])

def TimePlots():
    """Function that produces the violin plots for the times of the actual case and online and numerical ones.
    """
    #############VIOLIN PLOT
    pos = [0, 1, 2]
    data = [tf/3600, t_online/3600, t_num/3600]

    fig, axs=plt.subplots()
    violinz=axs.violinplot(data, pos, points=20, widths=0.3, showmeans=True)
    # set style for the axes
    axs.set_xticks([0, 1, 2])
    axs.set_xticklabels(labels = ['LHC', 'Online', 'Numerical'])
        
    # Make all the violin statistics marks red:
    for partname in ('cbars','cmins','cmaxes','cmeans'):
        vp = violinz[partname]
        vp.set_edgecolor('blue')
        vp.set_linewidth(1)

    # Make the violin body blue with a red border:
    for vi in violinz['bodies']:
        vi.set_facecolor('steelblue')
        vi.set_edgecolor('blue')
        vi.set_linewidth(0.25)
        vi.set_alpha(0.3)

    #tta_opt=sum(ta[:len(t_online)])
    axs.plot([],[], "k.", label="$T_\mathrm{LHC}=$"+"{:.1f} days".format(np.sum(tf)/(3600*24)+(np.sum(ta)/(3600*24))))
    axs.plot([],[], "k.", label="$T_\mathrm{Online}=$"+"{:.1f} days".format(np.sum(t_online)/(3600*24)+(np.sum(tta_opt)/(3600*24))))
    axs.plot([],[], "k.", label="$T_\mathrm{Numerical}=$"+"{:.1f} days".format(np.sum(t_num)/(3600*24)+(np.sum(ta)/(3600*24))))
    axs.set_xlabel('DataSet', fontsize=16)
    axs.set_ylabel(r'Fill Times [h]', fontsize=16)
    axs.set_title('20{} - Comparison between Optimisations'.format(year))
    #plt.legend(loc='best')

    plt.savefig('Online/ViolinPlot{}.pdf'.format(year))
    #plt.show()
    

def LumiEval():
    """Function that evaluates the needed integrated or total luminosities.

    Returns:
        Lint_online: Integrated online luminosity
        Ltot_online: Total online Luminosity
        Lint_num: Integrated numerical luminosity
        Ltot_num: Total Numerical luminosity"""
    #evaluating the integrated Luminosities and the total luminosity for 2017
    Lint_online=[]
    for i in range(len(t_online)):
        Linst=lambda x: fit(x, a[i], b[i], c[i], d[i])
        L_i=integrate.quad(Linst, 0, t_online[i])[0]
        Lint_online.append(L_i)

    Lint_online=np.array(Lint_online)
    Ltot_online=np.sum(Lint_online)

    Lint_num=[]
    for i in range(len(FillNumber)):
        Linst=lambda x: fit(x, a[i], b[i], c[i], d[i])
        L_i=integrate.quad(Linst, 0, t_num[i])[0]
        Lint_num.append(L_i)

    Lint_num=np.array(Lint_num)
    Ltot_num=np.sum(Lint_num)
    return Lint_online, Ltot_online, Lint_num, Ltot_num #microbarn^-1

def LumiPlots(Lint_online, Ltot_online, Lint_num, Ltot_num):
    """Function that produces the violin plots for the luminosities of the actual case and online and numerical ones.
    """
    
    #defining phyisics time
    T_num=sum(t_num)/(3600*24)+np.sum(ta)/(3600*24)
    T_real=np.sum(tf)/(3600*24)+(np.sum(ta)/(3600*24))
    T_online=np.sum(t_online)/(3600*24)+(np.sum(tta_opt)/(3600*24))

    #Evaluating total Real Luminosities
    Lmes_tot=np.sum(Lmes)


    #############VIOLIN PLOTs
    pos = [0, 1, 2]
    data = [Lmes, Lint_online/1e9 , Lint_num/1e9]

    fig, axs=plt.subplots()
    violinz=axs.violinplot(data, pos, points=20, widths=0.3, showmeans=True)
    # set style for the axes
    axs.set_xticks([0, 1, 2])
    axs.set_xticklabels(labels = ['LHC', 'Online', 'Numerical'])
        
    # Make all the violin statistics marks red:
    for partname in ('cbars','cmins','cmaxes','cmeans'):
        vp = violinz[partname]
        vp.set_edgecolor('blue')
        vp.set_linewidth(1)

    # Make the violin body blue with a red border:
    for vi in violinz['bodies']:
        vi.set_facecolor('steelblue')
        vi.set_edgecolor('blue')
        vi.set_linewidth(0.25)
        vi.set_alpha(0.3)

    axs.plot([],[], "k.", label="$\mathcal{L}^\mathrm{LHC}_\mathrm{Tot}/day=$"+"{:.4f}".format((Lmes_tot/1e9)/T_real)+"[$\mathrm{fb}^{-1}$]")
    axs.plot([],[], "k.", label="$\mathcal{L}^\mathrm{Online}_\mathrm{Tot}/day=$"+"{:.4f}".format((Ltot_online/1e9)/T_online)+"[$\mathrm{fb}^{-1}$]")
    axs.plot([],[], "k.", label="$\mathcal{L}^\mathrm{Num.}_\mathrm{Tot}/day=$"+"{:.4f}".format((Ltot_num/1e9)/T_num)+"[$\mathrm{fb}^{-1}$]")
    axs.set_xlabel('DataSet', fontsize=16)
    axs.set_ylabel(r'Integrated Luminosities [$\mathrm{fb}^{-1}$]', fontsize=16)
    axs.set_title('20{} Comparison between Optimisations'.format(year))
    #plt.legend(loc='best')

    plt.savefig('Online/ViolinPlot{}_Luminosity.pdf'.format(year))
    #plt.show()

    print('_________________________20{}___________________________________'.format(year))
    print('LHC Total Luminosity=', (Lmes_tot), '[fb^-1]')
    print('Online Total Luminosity=', Ltot_online/1e9, '[fb^-1]')
    print('Numerical Total Luminosity=', Ltot_num/1e9, '[fb^-1]')
    print('________________________________________________________________')


#MACRO
TimePlots()
Lint_online, Ltot_online, Lint_num, Ltot_num=LumiEval()
LumiPlots(Lint_online, Ltot_online, Lint_num, Ltot_num)
#print(np.sum(Lmes), Ltot_online, Ltot_num)