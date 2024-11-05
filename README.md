The Perovskite encapsulation comparison python code takes the data downloaded from the Material Zone Perovskite Project and compares the stability of different encapsulations

No pip install is needed - in order to get the out-put files all you need to do is run the code.
There are 3 input csv files needed to run this code: 'Perovsite database query.csv', 'additions_to_Encapsulation_info.csv' and 'additions_to_PCE_info.csv'.  
The out-put files are:
'Experiments_with_T80_greater_than_100.csv', 
'Experiments_with_T80_for_supp.csv', 
'Experiments_with_Encap_and_stability.csv', 
'frequency table logT80.csv', 
The t-test p-values for all 6 sub-devisions (there is a seperate csv file for each of the 6 sub devisions), 
Box plot graphs for each 6 sub-devision of the experiments.

The sub-devisions are: ALL, Ambient, 65 deg., Dark, Light, MPPT

Make sure you have the newest data downloaded from the Peroskite Project, the file should contain only encapsulated experiments (It's name should be 'Perovsite database query.csv' 
which is the default name of the downloaded data file)

What the code does:
Estimates the T80 (the value in column 'calculated Ts80').
Classifies each experiment into 9 encapsulation groups (Other, Glass + epoxy, Glass + UV glue, Glass + polymer, Glass + butyl rubber, Glass +no sealant, Foils and Polymers.
Compares the T80 of each classification by plotting the groups values and box plots. 

