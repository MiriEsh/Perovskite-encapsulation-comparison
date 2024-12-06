The Perovskite encapsulation comparison Python code takes the data downloaded
from the Material Zone Perovskite Database Project and compares the stability of different
encapsulations. 
No pip install is needed - in order to get the output analysis files: download the newest data from the Peroskite Project (the file should contain only encapsulated experiments and it's name should be 'Perovsite database query.csv') 
and run the code. There are 2 required input files: "Perovsite database query.csv" 
(contains the data as downloaded from the Perovskite Database Material Zone project) and
"Manual additions.csv" (contains corrections\additions to some of the rows in the MaterialZone
project, this file is required but all rows after the column names can be empty).
The output files are: "Experiments_with_T80_greater_than_100.csv",
"Experiments_with_T80_for_supp.csv", "Experiments_with_Encap_and_stability.csv",
"frequency table logT80.csv", 6 csv t-test info files (1 file for all experiments and a seperate
file for each of the 5 sub-devisions of the experiments), 6 box plot files (1 graph for all
experiments and a seperate graph for each of the 5 sub-devisions of the experiments).

What the code does:
1. Fixes encapsulation and lifetime data in the database according to the "Manual
additions" file. Cleans out all data from stability experiments that were not preformed
in the air.
2. Estimates the T80 of the remaining experiments if the T80 was not indicated in the
downloaded file data (the value in column "calculated T80" contains the T80 from the
database and, if it was not indicated in the data, than it contains the estimated
T80 values).
3. Classifies each experiment into one of the following 8 encapsulation groups: Other,
Glass + epoxy, Glass + UV glue, Glass + polymer, Glass + butyl rubber, Al2O3, Foils
and Polymers.
4. Compares the T80 of 7 groups (excluding “Other”) by plotting the group’s
T80 in box plots and performing a t-test for every possible combination of 2
groups.

