The Perovskite encapsulation comparison python code takes the data downloaded from the Material Zone Perovskite Project and compares the stability of different encapsulations.

No pip install is needed - in order to get the output analysis files all you need to do is run the code.
There are 2 input csv files needed: 'Perovsite database query.csv' (contains the data as downloaded from the Material Zone project) and 'Manual additions.csv' (contains corrections\additions to any of the lines in the MaterialZone project).  
The output files are:
'Experiments_with_T80_greater_than_100.csv', 
'Experiments_with_T80_for_supp.csv', 
'Experiments_with_Encap_and_stability.csv', 
'frequency table logT80.csv', 
6 csv t-test info files (1 file for all experiments and a seperate file for each of the 5 sub-devisions of the experiments), 
6 box plot files (1 graph for all experiments and a seperate graph for each of the 5 sub-devisions of the experiments)

The 5 sub-devisions are: experiments at ambient temperatures, experiments at elevated tempuratures, experiments in dark conditions, experiments in light conditions and experiements at MPPT

Download the newest data from the Peroskite Project, the file should contain only encapsulated experiments and it's name should be 'Perovsite database query.csv' 
(which is the default name of the downloaded data)

What the code does:
1) Fixes encapsulation and PCE data in the database according to the 'Manual additions' file. Cleans out all experiments that were not preformd in the Air.
2) Estimates the T80 of the remaining experiments if the T80 was not in the downloaded file (the value in column 'calculated T80' contains the T80 from the database and if it was not in the database than incontains the estimated T80).
3) Classifies each experiments into one of the following 8 encapsulation groups: Other, Glass + epoxy, Glass + UV glue, Glass + polymer, Glass + butyl rubber, Foils and Polymers.
4) Compares the T80 of 7 groups (Glass + epoxy, Glass + UV glue, Glass + polymer, Glass + butyl rubber, Foils and Polymers) by plotting the groups T80 and box plots and preforming a t-test for every possible combination of 2 groups.

