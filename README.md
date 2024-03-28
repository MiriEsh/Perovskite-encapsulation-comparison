The Perovskite encapsulation comparison python code takes the data downloaded from the Material Zone Perovskite Project and compares the stability of different encapsulations

No pip install is needed - all you need to do is run the code and you will get the outputs of:
Experiments_with_Ts80_for_supp.csv
Experiments_with_Encap_and_stability.csv
frequency table Ts80.csv
Box plot graphs for each 6 sub-devision of the experiments (ALL, Ambient, 65 deg., Dark, Light, MPPT)
The t-test p-values for all 6 sub-devisions (seperate csv for each one)

Make sure you have the newest data downloaded from the Peroskite Project, the file should contain only encapsulated experiments (It's name should be 'Perovsite database query.csv' 
which is the default name of the downloaded data file)

What the code does:
Estimates the Ts80 (the value in column 'calculated Ts80').
Classifies each experiment into one of the 7 encapsulation groups (Glass + epoxy, Glass + UV glue, Glass + polymer, Glass + butyl rubber, Foils and Polymers.
Compares the Ts80 of each classification by plotting the groups values and box plots. 

