# The top script which calls and checks 
#TODO Email NCBI about the faulty files, KEEP UNTIL LATER SO WE KNOW WE HAVE ONE FAULTY.
#TODO Perform for larger dataset, investigate faulty or strange records
#TODO Fix the log file gets the date - DONE
#TODO Make it compatible for bash commands - DONE
#TODO Add optional verbosity as bash commands - 
#TODO Write the description of the program for the terminal parser. 
#TODO The new dnaA investigation things - 
#TODO Need a check if biopython is downloaded?

import argparse
import sys
import os
import logging
import datetime
from combined_scripts import get_rRNA_intervals as get_rrna
sys.path.insert(0, '/skewDB/')
from skewDB import dowloading_filtered_csvFile as download_filtered
from rrna_leading_lagging import rrna_lead_lag as rll

"""
parser = argparse.ArgumentParser(description='What the process does.')

# 4 positional arguments which are mandatory
parser.add_argument('csv_path', metavar='csv_path',
                    help='The csv path to what the filtered data from the SkewDB will be name. Example: C:/Users/Felix/Documents/FilteredDataFile.csv)

parser.add_argument('email', metavar='email',
                    help='The email used for the NCBI account of the user. Must be the account as the API-key is for. Example: firstname.lastname@gmail.com')

parser.add_argument('api_key', metavar='api_key',
                    help='The api key for the NCBI account. API-key must be used and must match the email for the account.' Example: 1y4a5e5641h73645fg73759384ad485lot05)

parser.add_argument('local_storage_path', metavar='local_storage_path',
                    help='The local storage path to a where the NCBI records will be stored, must have sufficient space. Approximately 70 GB' Example: D:/)

# Optional arguments

#Parsing all arguments
args = parser.parse_args()
"""

def start(csv_path, email, api_key, local_storage_path):
    # Start logging
    if __name__ == "__main__":
        # Starting main log file
        now = datetime.datetime.now()
        start_datetime = now.strftime('%Y-%m-%d__%H%M')
        log_file_name = str("logfile_" + start_datetime + ".log")
        print(log_file_name)
        print(type(log_file_name))
        logging.basicConfig(level=logging.DEBUG, filename=log_file_name, filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
        
    

    #Check if csv is downloaded & filtered 
    if not os.path.isfile(csv_path):
        download_filtered.run_download_filtered_csvfile(csv_path)
        print("csv filtered downloaded")


    #Get the rRNA interval dict
    rrna_dict = get_rrna(csv_path, email, api_key, local_storage_path)    


    # Send to rrna leading lagging script
    rll(csv_path, rrna_dict)

    print("Everything is done")
#start(args.csv_path, args.email, args.api_key, args.local_storage_path)



csv_path = 'C:/Users/Felix/Documents/FilteredDataFile.csv'
email = "Felix.wae@gmail.com"
api_key = "7b4a5e9841f79495be73767323ad485fda08"
local_storage_path = 'D:/'
start(csv_path, email, api_key, local_storage_path)
