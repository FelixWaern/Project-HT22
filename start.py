# The top script which calls and checks 
#TODO Email NCBI about the faulty files, KEEP UNTIL LATER SO WE KNOW WE HAVE ONE FAULTY.
#TODO Perform for larger dataset, investigate faulty or strange records

#TODO The new dnaA investigation things 
#TODO Fix comments. Can add more if possible. 
#TODO Run the entire thing
#TODO Report
#TODO Won't open. 

import importlib.util
import pandas as pd
import argparse
import sys
import os
import logging
import datetime
from combined_scripts import get_rRNA_intervals as get_rrna
sys.path.insert(0, '/skewDB/')
from skewDB import dowloading_filtered_csvFile as download_filtered
from rrna_leading_lagging import rrna_lead_lag as rll
from download_gcfits import run_download_gcfits as rdg
from plots import plotting_graphs


parser = argparse.ArgumentParser(description='The start.py script downloads chromosomes data from the SkewDB by Bert Hubert. It filters it, calculates leading and lagging strand, and passes it '
                                "to to a subscript which fetches the NCBI records corresponding to the remaining chromosomes. "
                                "The NCBI record contains rRNA intervals which are compared to the leading and lagging strand intervals. "
                                "If the primary rRNA intervals do not overlap with leading they are stored in a logfile and reported. And vice versa for complementary rRNA and the lagging strand.")

# 4 positional arguments which are mandatory
parser.add_argument('csv_path', metavar='csv_path',
                    help='The csv path to what the filtered data from the SkewDB will be name. Example: C:/Users/Felix/Documents/FilteredDataFile.csv')

parser.add_argument('email', metavar='email',
                    help='The email used for the NCBI account of the user. Must be the account as the API-key is for. Example: firstname.lastname@gmail.com')

parser.add_argument('api_key', metavar='api_key',
                    help='The API-key for the NCBI account. API-key must be used and must match the email for the account. Example: 1y4a5e5641h73645fg73759384ad485lot05')

parser.add_argument('local_storage_path', metavar='local_storage_path',
                    help='The local storage path to a where the NCBI records will be stored, must have sufficient space. Approximately 70 GB. Example: D:/ ')

# Optional arguments
parser.add_argument('-v', '--verbose',
                    action='store_true',
                    help = 'Toggle verbose logging. Changes the amount of information in the output log file.' )  # on/off flag

parser.add_argument('-l', '--a_list',
                    nargs='+',
                    help = 'Option for running the script for a specific list of records. Example: -l NC_034600.1 NC_005780.1' )  # Get a accession nr list  

parser.add_argument('-t', '--test_set',
                    action='store_true',
                    help = 'Test the script by running it on a artificial dataset. Specify test_chromosome_set.csv location instead of csv_path' )  # Runs the script on a test data set. 

#Parsing all arguments
args = parser.parse_args()



def start(csv_path, email, api_key, local_storage_path, verbose=False, a_list=[], test_set=False):
    """ Starting the script which matches the strands for the chromosomes from SkewDB against the 16S rRNA locations from
    NCBI records corresponding to the chromosomes. The script is started from the command line.
    Use python start.py -help for instructions on setup and running the script.
    """

    print("--------------Starting SkewDB rRNA interval match--------------")

    # Start logging
    if __name__ == "__main__":
        # Starting main log file
        now = datetime.datetime.now()
        start_datetime = now.strftime('%Y-%m-%d__%H%M')
        log_file_name = str("logfile_" + start_datetime + ".log")
        logging.basicConfig(level=logging.DEBUG, filename=log_file_name, filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info(f"Logging verbose set to {verbose}") 

    #Check if biopython is installed. 
    name = 'Bio'
    if name in sys.modules: #Check if already loaded. 
        if verbose == True:
            logging.debug(f"\n --------Biopython already loaded---------")
        print(f"{name!r} already in sys.modules") 
    elif (spec := importlib.util.find_spec(name)) is not None: #Check if package exists on machine
        if verbose == True:
            logging.debug(f"\n --------Biopython is installed--------- ")
        print(f"{name!r} is installed")
    else:
        print(f"can't find the {name!r} module, not installed")
        logging.warning("BIOPYTHON NEEDS TO BE INSTALLED BEFORE RUNNING START.PY")
        sys.exit()

    # Checks if the processes are to be run using the test set available. 
    if test_set == True:
        print("Running script on test set.")

    #Check if csv is downloaded & filtered 
    # Skip if running on test set
    if test_set == True:
        pass
    else: 
        if not os.path.isfile(csv_path):
            if verbose == True:
                logging.debug(f"\n --------Filtered csv file not found--------- \n run_downloaded_filtered_csvfile input: {csv_path} ")
            download_filtered.run_download_filtered_csvfile(csv_path)
            print("csv filtered downloaded")


    # Get the rRNA interval dict
    # Skip if running on test set
    if test_set == True:
        pass
    else:
        if verbose == True:
            logging.debug(f"\n --------parameters into rrna_dict()--------- \n csv_path = {csv_path} \n email = {email} \n api_key = {api_key} \n local_storage_path = {local_storage_path}")
        rrna_dict = get_rrna(csv_path, email, api_key, local_storage_path, a_list, verbose) #rna dict is list. 0 is the position and 1 is the locus_tag. And change inputs.
    
    
    if test_set == True:
        # Loading both test rrna csv and locus csv and loading it as a dataframe
        # Starting the matching process using a dictionary of all chromosomes and rRNAs, and strand data
        temp_rna_path = "test_rna_dict.csv"
        df = pd.read_csv(temp_rna_path)
        skip = df.iloc[:,0]
        temp_rrna_dict = {}
        locus_dict = {}
        for i in range(0, len(df)):
            rna_row = []
            locus_tag_row =[]
            for x in df.loc[i]:
                if x == skip[i]:
                    pass
                else:
                    if (x[0]) == "l": # Check if locus or rRNA position. All locus for test set should start with l.
                        locus_tag_row.append(x)
                    else:
                        rna_row.append(x)
            #rna_row = rna_row[1:]
            temp_rrna_dict[i + 1] = rna_row
            locus_dict[i + 1] = locus_tag_row
        rrna_dict = [temp_rrna_dict, locus_dict]
        rll(csv_path, rrna_dict) # Need to change the rll so that True uses only test set. Add if statment in rll. 

    else:
        # Starting the matching process using a dictionary of all chromosomes and rRNAs, and strand data
        if verbose == True:
            logging.debug(f"\n parameters into rrna_dict: \n csv_path = {csv_path} \n rrna_dict = rrna_dict, too long to display")
        rll(csv_path, rrna_dict) # Need to specicfy which in the list now. 


    # Creating graphical representation. 
    if test_set == True:
        logging.debug(f"\n --------Graphical representation is not available for for test set-------------- ")
        pass
    else: 
        file_path = csv_path.rstrip("FilteredDataFile.csv")
        gcfit_path = file_path + "gcfits"
        #Check if gc_fits.csv is downloaded
        if not os.path.isdir(gcfit_path):
            print("Downloading gcfits")
            if verbose == True:
                logging.debug(f"\n --------Filtered csv file not found--------- \n run_download_gcfits input: {gcfit_path} ")
            #download_filtered.run_download_filtered_csvfile(csv_path)
            rdg(file_path)
        # Run the graphical representation script
        print("Plotting graphs")
        plotting_graphs(csv_path)

    print("Everything is done")
start(args.csv_path, args.email, args.api_key, args.local_storage_path, args.verbose, args.a_list, args.test_set)


"""
csv_path = 'C:/Users/Felix/Documents/FilteredDataFile.csv'
email = "Felix.wae@gmail.com"
api_key = "7b4a5e9841f79495be73767323ad485fda08"
local_storage_path = 'D:/'
start(csv_path, email, api_key, local_storage_path)
"""