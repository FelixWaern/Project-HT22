# Project-HT22
GitHub Page for the Project _"Exploring the effects of replication on chromosomal features and orientation of genes"_ in the course Applied Bioinformatics HT 22

## Project description
Bacteria (and some archaea) are replicated bidirectionally, from the origin of replication (Ori) to the terminus (Ter). This property structures the whole genome, impacting – among others – the nucleotide composition and the location and orientation of genes, emphasizing the importance of these two locations.

Since the beginning of the genomic era, researchers have started to use these properties to identify Ori and Ter. The most reliable indicator is the GC-skew, i.e. the difference between the number of Gs and Cs on one strand, but other sources of evidence have been used (some marker sequences, gene orientation, etc.). A newly published database, SkewDB (https://skewdb.org/), proposes a robust GC-skew-based method to predict the origin and terminus of replication in most chromosomes. This opens up whole new avenues of research, allowing to correlate chromosomal features with replication.

## Aims
In the course of this project we aimed firstly to test if the hypothesis that rRNA genes are co-oriented with replication is still true (Guy & Roten 2004). The second task is to investigate the cases where non-colocalization of dnaA with Ori is related to taxonomic distribution.

## How to run the pipeline
Before any python script can be run for the validation of co-orientation of the rRNA process the Biopython package needs to be installed using for example pip. If not installed, running the script will result in an error in the log file and the run is interrupted. 

The python script used for the validation of co-orientation of the rRNA process is started from the command line by having the folder containing start.py as the current directory and starting python and calling the script with the appropriate arguments.  

Example of command: 

_python start.py <path to FilteredDataFile.csv> <firstname.lastname@gmail.com> <7b4a5e9841f79495be77491+223ad485fda08> <D:/>_
  
All of the above arguments are mandatory.
  
- The first argument is the path to what the filtered data is named or is going to be named.
- The second argument is the mail used for the user's NCBI account.
- The third argument is the API-key from NCBI for that specific user.
- The fourth argument is the path to local storage with at least 120 GB of free storage. 

There are also optional commands for start.py:
- The flag -v is used to toggle verbose logging which will increase the amount of information in the log file. This will cause it to display most information about the run which can be useful when debugging.
- The -l flag will allow the user to run the script using specified accession numbers. 

Before retrieving any data or running any functions the script will first create a log file with the date and time and check if Biopython is installed. The first function will retrieve data from SkewDB unless a filtered version already exists locally on the computer.
