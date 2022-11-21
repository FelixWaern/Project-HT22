# The top script which calls and checks 
#TODO Email NCBI about the faulty files, 

# create function which takes different settings as input
# setting, could be location fo different files and paths and so forth

#download location for csv
import sys
import os
import logging
sys.path.insert(0, '/skewDB/')
from skewDB import dowloading_filtered_csvFile as download_filtered
from combined_scripts import get_rRNA_intervals as get_rrna

def start(csv_path, email, api_key, local_storage_path):
    # Start logging
    if __name__ == "__main__":
        logging.basicConfig(level=logging.DEBUG, filename="faulty_NCBI_records", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    

    #Check if csv is downloaded & filtered 
    if not os.path.isfile(csv_path):
        download_filtered(csv_path)
        print("csv filtered downloaded")

    #Get the rRNA interval dict
    rrna_dict = get_rrna(email, api_key, local_storage_path)    

    # Send to Saras script


csv_path = 'C:/Users/Felix/Documents/FilteredDataFile.csv'
email = "Felix.wae@gmail.com"
api_key = "7b4a5e9841f79495be73767323ad485fda08"
local_storage_path = 'D:/'
start(csv_path, email, api_key, local_storage_path)
