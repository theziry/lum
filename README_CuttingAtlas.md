# How does CuttingAtlas work?
#### @GiuliaFaletti
The CuttinAtlas software allows the user to transform ATLAS ROOT lumi data into txt files, which store (or should store if you are downloading this project): 
1. **20{}**: ATLAS data of the considered year({16-17-18}) stored in ROOT files (named after the specific ATLAS run);
2. **20{}Data**: Atlas data (***lb***-> luminosity block number; ***dt***-> luminosity block duration in seconds; ***nbx***-> number of colliding BCIDs; ***grl***-> good run list flag to reject bad ATLAS lumi block; ***\<current year algorithm name>***-> integrated luminoisty value evaluated with a specific algorithm for the LUCID detector; ***\<current year algorithm name>_err***-> errors on the previous measurements.) saved in txt files, named after the fill number to which they refer; 
3. **OldData_and_Utilities**: Folder that stores old trials and data for the task;
4. **FillData{}.txt**: txt files which store the preselected fills (selection made considering ...);
5. **CuttingATLAS.py**: Script to extract data from the ROOT files (**20{}**) and write them down on txt files (**20{}Data**).

## What is needed?:
To compile this script a ROOT and PyROOT enviroment is necessary. To do so with conda:

$ ```conda create -n my_root_env root -c conda-forge```

To enter the environment:
$ ```conda activate my_root_env```

The first time you enter the environment, you should add the conda-forge channel to the search list (otherwise, you will have to add -c conda-forge every time you install or update something):
$ ```conda config --env --add channels conda-forge```

To leave the environment:
$ ```conda deactivate```

## How does CuttingATLAS work?
First, load and compare ATLAS fills with preselected fills. Once the intersection between the two sets is obtained, it opens the ROOT files corresponding to the chosen fills and saves the necessary data in text files of the type:
|lb|dt|nbx|grl|algorithm|algorithm_err|
|---|---|---|---|---|---|
| ... | ... | ... | ... | ... | ... | 

