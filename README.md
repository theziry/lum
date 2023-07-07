# How does ONLINEOPTIMISER work?
## @GiuliaFaletti
The OnlineOptimiser allows the user to find the optimal values of fill times for the LHC.
The folder contains:
* **ATLAS**: Folder with Massi files data (istantaneous luminosity and integrated ones) and graphs of luminosity evolutions;
* **Cutting_Fitting**: Folder that stores the results of the Cutting_Fitting.py script;
* **MPL**: Folder that stores the most probable luminosity model of each year;
* **NumericalOptimization**: Folder that stores the results of the numerical optimiser;
* **Online**: Folder that stores the results of the online optimisations (Future fill times and optimal fill times);
* **Cutting_Fitting.py**: Script that performs the necessary cuts and fits for all selected fills also considering partial knowledge of the fitting parameters;
* **LoadData.py**: Script that loads the needed data from FillData.xlsx;
* **MPL.py**: Most Probable luminosity model;
* **NumericalOptimisation.py**: Numerical optimisation of the Run 2 control room operations;
* **Online.py**: Script that performs the online optimisation;
* **FillData.xlsx**: Stores the preselected fills number, times and turnaround times for the studied years (2016-2017-2018);
* **OptimisationGraphs.py**: Scripts that produces the violin plots for times and luminosities for the studied years (which in this case ca be only 2017 and 2018).

The same files and folders with the denomination _New_ represent the same files as above, but, in this case, we considered the offline-processed data coming from the LUCID detector. The only differences are:
* **LumiData**: folder in which are store all the LUCID data and some reports on the problematic fills, that must be solved;
* **OnlineNewMC.py**: Same algorithm as OnlineNew.py but with the presence of a MonteCarlo function that can be eventually called to simulate missing fills.


# How to run a simulation?

1. Run ***Cutting_Fitting.py*** or ***Cutting_FittingNew.py*** for 2016, 2017 and 2018 changing the year parameter.
2. Run ***NumericalOptimisation.py*** or ***NumericalOptimisationNew.py*** choosing the year and the plotting desired configuration.
3. Run ***Online.py*** or ***OnlineNew.py*** (or ***OnlineNewMC.py***) for the desired year (i.e. in this case only 2017 or 2018).
4. Run ***OptimisationGraphs.py*** or ***OptimisationGraphsNew.py*** for the desired year (i.e. in this case only 2017 or 2018).



# A more detailed look
## More on LoadData.py
1. ***Create_DataSet():*** Reads the Excel file and creates the corret datasets. 
2. ***DataToLists(data1, data2, data3):*** Transfroms Pandas Dataframes into Lists.
3. ***FromListsToArrays(data_16, data_17, data_18):*** Transfroms Lists into arrays.
4. ***TotalDataSet(data_16, data_17, data_18):*** Creates the total dataset and transforms the inital list into a dataframe and an array.
5. ***PartialDataSets(data_16, data_17, data_18):*** Creates the partial dataseta and transforms the inital lists into dataframea and arrays.
6. ***Data():*** Generate the whole set of data sample needed.
7. ***loadFill():*** Reads the Excel file and creates the corret dataset lists.   
8. ***Data_sec(array16, data_ta16, data_tf16, array17, data_ta17, data_tf17, array18, data_ta18, data_tf18):*** Transform the data from hours to seconds.
9. ***FillNumber():*** Creates the arrays that contain the number of fills' list.

10. ***MeasuredLuminosity():*** Evaluates the measured luminosity from ATLAS data.
11. ***savingFillNumber_txt(year, FillNumber):*** saving fill number arrays into txt files.

12. ***AtlasData(text, year, grl=True):*** Function that sets the data extracted from Atlas ROOT files for them to be usable in optimization algorithms, allowing all data or only those with positive flags to be selected in the good atlas run list.
13. ***CuttedData(year, text):*** Function that read the saved cutted and fitted data from txt files.
14. ***CuttedNewData(year, text):*** Function that read the saved cutted and fitted new data from txt files.
15. ***OldCuttedData(year, text):*** Function that read the saved cutted and fitted data from txt files.

16. ***AtlasCuttedData(year, text):*** Function that read the saved cutted and fitted new data from txt files.



## More on NumericalOptimisation.py
Select the year to be analize setting the year variable (i.e. 16,17,18). At this point, the script will produce an optimisation of what actually happened during that year, saving the optimal fill values to file. If you then set the plot variable to True, then the script will also produce two graphs. The first graph will represent a comparison between the histogram of the actual values and that of the optimised values. The second graph will always produce similar histograms but for the case of fills' integrated luminosities.

## More on MPL.py

1. ***CuttingAlgorithm*** (current year, current fill): Function that performs the necessary cut on the current fill returning the fit times in seconds, the istantaneous luminosities and the fitting parameters.
2. ***PeakLumi*** (current year, current fill): Function that returns the peak luminosity of the current fill.
3. ***NormMPL*** (peak lumi, current year, current fill): Function that evaluates the Most Probable Luminosity with normalization to the current fill's peak luminosity, returning the most probable times and the most probable luminoisties.
4. ***MPL*** (current year, current fill): Function that evaluates the Most Probable Luminosity, returning the most probable times and the most probable luminoisties.

