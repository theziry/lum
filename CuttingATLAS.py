#########################################################
#      @Giulia Faletti                                  #
#   From ATLAS ROOT data to python objects              #
#########################################################

#from collections import OrderedDict
#from sortedcollections import SortedDict
import numpy as np
import csv
import ROOT


def FillNumber():
     """Creates the arrays that contain the number of fills' list reading data from 3 .txt files.

     Returns:
         FillNumber16, FillNumber17, FillNumber18: array
     """
     #defining lists to save data
     FillNumber16=[]
     FillNumber17=[]
     FillNumber18=[]
     #read the tvs file and saving data in the lists
     f=open('FillData16.txt',"r")
     lines=f.readlines()
     for x in lines:
          FillNumber16.append(int(float(x.split(' ')[0])))
          
     f1=open('FillData17.txt',"r")
     lines1=f1.readlines()
     for x1 in lines1:
          FillNumber17.append(int(float((x1.split(' ')[0]))))
          
     f2=open('FillData18.txt',"r")
     lines2=f2.readlines()
     for x2 in lines2:
          FillNumber18.append(int(float(x2.split(' ')[0])))
          
     return FillNumber16, FillNumber17, FillNumber18



#setting the current year
year=18

#setting the preselected-fill numbers
FillNumber16, FillNumber17, FillNumber18=FillNumber()

if year==16:
    MyFillNumber=FillNumber16
elif year==17:
    MyFillNumber=FillNumber17
elif year==18:
    MyFillNumber=FillNumber18


#defining lists to save data
ATLAS_Data=[]
FillNumberATLAS=[]

#read the tvs file and saving data in the lists
with open('20{}/lumi_20{}.tsv'.format(year, year), "r") as file:
    tsv_file = csv.reader(file, delimiter="\t")
    for row in tsv_file:
        if row[0].isalpha():
            continue
        ATLAS_Data.append(int(row[0]))
        FillNumberATLAS.append(int(row[4]))

#from lists to array
MyFillNumber=np.array(MyFillNumber)
FillNumberATLAS=np.array(FillNumberATLAS) 
         
#Choosing only the preselected fills
newKeys=[]
newValues=[]
for el in MyFillNumber:
    try:
        ind=np.where(FillNumberATLAS==el)[0][0]
        #print(ind)
        newKeys.append(ATLAS_Data[ind])
        newValues.append(FillNumberATLAS[ind])
    except:
        continue

#from lists to array
newKeys=np.array(newKeys)
newValues=np.array(newValues)



#Loading empty data files
for i in range(len(newValues)):
    fill=newValues[i]
    with open('20{}Data/{}_Data.txt'.format(year, fill), 'w') as f:
        f.write('')
        f.close() 
        
#saving data from the .ROOT files into .txt file whose title contains the fill number
for i in range(len(newValues)):
    fill=newValues[i]
    runAtlas=newKeys[i]
    infile=ROOT.TFile("20{}/PerLB_{}.root".format(year, runAtlas), "READ")
    t=infile.Get("t")
    N=t.GetEntries()
    #print(N)
    for i in range(N):
        t.GetEntry(i)
        #saving selected data in txt files: lb, dt, nbx, grl, Algorithm for Lumi Measurements, Algorithm errors
        #lb-> Luminosity block number
        #dt-> Luminosity block duration in seconds
        #nbx-> Number of colliding BCIDs
        #grl-> Good run list flag to reject bad ATLAS lumi block
        with open('20{}Data/{}_Data.txt'.format(year, fill), 'a') as f:
            f.write(str(t.lb))
            f.write('\t')
            f.write(str(t.dt))
            f.write('\t')
            f.write(str(t.nbx))
            f.write('\t')
            f.write(str(t.grl))
            f.write('\t')
            #select the correct algorithm according to different years
            if year==16 or year==17:
                f.write(str(t.LUCID_HITOR_BI))
                f.write('\t')
                f.write(str(t.LUCID_HITOR_BI_err))
                f.write('\t')
            elif year==18:
                f.write(str(t.LUCID_C12_Bi2))
                f.write('\t')
                f.write(str(t.LUCID_C12_Bi2_err))
                f.write('\t')
            f.write('\n')
            f.close() 


#Pro
#checking for multi-fill ATLAS runs
#for el in newValues:
   #index=np.where(newValues==el)[0]
    #print(index)

#creating a dictionary to take in to account the relation between newKeys and newValues 
#relation=dict(zip(newKeys, newValues))
#relation1=OrderedDict(SortedDict(relation.items()))
#print(relation)