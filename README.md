# Project-HT22
GitHub Page for the Project "Exploring the effects of replication on chromosomal features and orientation of genes" in the course Applied Bioinformatics HT 22
# How to run the pipeline
Before any python script can be run for the validation of co-orientation of the rRNA process the Biopython package needs to be installed using for example pip. If not installed, running the script will result in an error in the log file and the run is interrupted. 

The python script used for the validation of co-orientation of the rRNA process is started from the command line by having the folder containing start.py as the current directory and starting python and calling the script with the appropriate arguments.  

Example of command: 
python start.py <path to FilteredDataFile.csv> <firstname.lastname@gmail.com> <7b4a5e9841f79495be77491+223ad485fda08> <D:/>
All of these arguments are mandatory.
The first argument is the path to what the filtered data is named or is going to be named
The second argument is the mail used for the user's NCBI account
The third argument is the API-key from NCBI for that specific user.
The fourth argument is the path to local storage with at least 120 GB of free storage. 

There are also optional commands for start.py
The flag -v is used to toggle verbose logging which will increase the amount of information in the log file. This will cause it to display most information about the run which can be useful when debugging.
The -l flag will allow the user to run the script using specified accession numbers. 

Before retrieving any data or running any functions the script will first create a log file with the date and time and check if Biopython is installed. The first function will retrieve data from SkewDB unless a filtered version already exists locally on the computer.
