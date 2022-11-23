# The top script which calls and checks 
#TODO Email NCBI about the faulty files, perform for larger dataset, investigate faulty or 
# strange record, fix the log file gets the date, make it compatible for bash commands 
#TODO The new dnaA investigation things

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
                    help='The csv path')
parser.add_argument('email', metavar='email',
                    help='The email')
parser.add_argument('api_key', metavar='api_key',
                    help='The api key')
parser.add_argument('local_storage_path', metavar='local_storage_path',
                    help='The local storage path')

# Optional arguments

#Parsing all arguments
args = parser.parse_args()
"""

def start(csv_path, email, api_key, local_storage_path):
    # Start logging
    if __name__ == "__main__":
        now = datetime.datetime.now()
        start_datetime = now.strftime('%Y-%m-%d %H:%M')
        logging.basicConfig(level=logging.DEBUG, filename="logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
                        
        """Function that searches for 16S rRNA genes"""
        # Get start date and time of analysis
        
        

        # Create log file
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)
        
    

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
