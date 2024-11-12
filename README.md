The Perovskite encapsulation comparison python code takes the data downloaded from the Material Zone Perovskite Project and compares the stability of different encapsulations.

No pip install is needed - in order to get the out-put files all you need to do is run the code.
There are 2 input csv files needed to run this code: 'Perovsite database query.csv' (contains the data as downloaded from the Material Zone project, 'Manual additions.csv' (contains corrections\additions to any of the lines in the Material Zone project.  
The out-put files are:
'Experiments_with_T80_greater_than_100.csv', 
'Experiments_with_T80_for_supp.csv', 
'Experiments_with_Encap_and_stability.csv', 
'frequency table logT80.csv', 
6 csv t-test info files (1 file for all experiments and a seperate file for each of the 5 sub-devisions of the experiments), 
6 box plot files (1 graph for all experiments and a seperate graph for each of the 5 sub-devisions of the experiments)

The 5 sub-devisions are: Ambient experiments, experiments at 65 degrees, experiments in dark conditions, experiments in light conditions, experiements at MPPT

Download the newest data from the Peroskite Project, the file should contain only encapsulated experiments and it's name should be 'Perovsite database query.csv' 
(which is the default name of the downloaded data)

What the code does:
1) Estimates the T80 if it is not in the downloaded file (the value in column 'calculated T80').
2) Classifies each experiments into 9 encapsulation groups: Other, Glass + epoxy, Glass + UV glue, Glass + polymer, Glass + butyl rubber, Glass +no sealant, Foils and Polymers.
3) Compares the T80 of 7 groups (Glass + epoxy, Glass + UV glue, Glass + polymer, Glass + butyl rubber, Foils and Polymers) by plotting the groups T80 and box plots and preforms t-test for each and every combination of 2 groups.

