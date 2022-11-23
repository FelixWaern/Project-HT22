# The top script which calls and checks 
#TODO Email NCBI about the faulty files, fix so only one log file with date


import sys
import os
import logging
import datetime
from combined_scripts import get_rRNA_intervals as get_rrna
sys.path.insert(0, '/skewDB/')
from skewDB import dowloading_filtered_csvFile as download_filtered
from rrna_leading_lagging import rrna_lead_lag as rll


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
    
csv_path = '/Users/saralindberg/Documents/Applied_bioinformatics/Code/FilteredDataFile.csv'
email = "lindberg.sara2@gmail.com"
api_key = "bddd70a82f6fd9c645410095b4a732aa7c08"
local_storage_path = '/Volumes/LIONELGUY/'
start(csv_path, email, api_key, local_storage_path)