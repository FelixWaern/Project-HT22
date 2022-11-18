from Bio import SeqIO
import datetime
import logging


def find_rrna():
    """Function that searches for 16S rRNA genes"""
    # Get start date and time of analysis
    now = datetime.datetime.now()
    start_datetime = now.strftime('%Y-%m-%d %H:%M')

    # Create log file
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(f'warnings_rrna_{start_datetime}.log', 'w', 'utf-8')
    root_logger.addHandler(handler)

    # Local variables
    rrna_16s = []
    rrna_other = []

    with open("GCF_000006765.1_ASM676v1_genomic.gbff") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == "tRNA":
                    for product in feature.qualifiers.get("product"):
                        if "16S" in product:
                            rrna_16s.append(str(feature.location))       
                        elif "RNA" in product:
                            rrna_other.append(product)   
                         
    # Print warning or info to log file
    if len(rrna_16s) == 0:
        if len(rrna_other) == 0:
            logging.warning("---------- WARNING ----------")
            logging.warning(f" \nNo rRNA was found for {record.id}") 
        else:
            logging.warning("---------- WARNING ----------")
            logging.warning(f" \nNo 16S rRNA genes were found for {record.id}, but these products were found:") 
            for e in rrna_other:
                logging.warning(f"{e}") 
    else:
        logging.info("---------- INFO ----------")
        logging.info("All good!")
        
    # Print to command line
    print("---------- SUMMARY ----------")
    print("Search for rRNA genes is done, for details see 'warnings_rrna_YYYY-MM-DD HH:MM.log' file")


find_rrna()