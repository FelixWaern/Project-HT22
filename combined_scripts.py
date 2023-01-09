# Combining the fetching data & csv_filtered
def get_rRNA_intervals(csv_path, email, api_key, local_storage_path, a_list, verbose=False):
    """ Iterating over all chromosomes and retreiving rRNAs as a dictionary with accession number as keys and rRNA locations as values.
        Also fetching the locus tags for all rRNAs as a dictionary using accession number as keys and rRNA locus tags as values. 
    """
    import time
    from collections import defaultdict
    import sys
    import logging
    sys.path.insert(0, '/skewDB/')
    from skewDB import fetching_data as fd
    sys.path.insert(0, '/NCBI_DATA_FETCH/')
    from NCBI_DATA_FETCH import main_script as ms
    import platform

    # Fetching the SkewDB data as a dataframe
    org_df = fd.fetch_csv_as_df(csv_path) 
    
    if a_list != None:
        df = org_df.loc[org_df['name'].isin(a_list)]
    else:
        df = org_df
    
    i = 0
    j = 1
    batch = []
    faulty = []
    dict = {}
    locus = {}
    t_tot = []
    no_16s = []
    t_fin_1 = time.time()
    # Fetching the rRNA data in batches of 10 from NCBI since 10 is the maximum number 
    # of fetches per second using an ordinary API-key from NCBI
    for index, row in org_df.iterrows():
        if i == 9:
            batch.append(row["name"])
            if verbose == True:
                logging.debug(f"\n --------Batch {j} sent to batch_operator--------- \n Contains NCBI records: {batch} ")
            t0 = time.time()

            # Sending a batch of 10 accession number to the fetching script
            res = ms.batch_operator(batch, faulty, email, api_key, local_storage_path, no_16s ,verbose )
            dict.update(res[0])
            locus.update(res[1])
            i = 0
            batch = []
            t1 = time.time()
            total = t1-t0
            t_tot.append(total)  

            print("")
            print("Batch:",j, "done!")
            print("Batch:",j, "took:", total, "seconds!")
            print("Estimated time left: ", (((sum(t_tot)/len(t_tot))*2771)-((sum(t_tot)/len(t_tot))*j))/60 ,"minutes")
            print("Estimated mean total time: ", ((sum(t_tot)/len(t_tot))*2771)/60,"minutes")
            print("Estimated percentage done", ((1-(((sum(t_tot)/len(t_tot))*2771)-((sum(t_tot)/len(t_tot))*j))/((sum(t_tot)/len(t_tot))*2771))*100), "%")
            print("Current faulty records: ", faulty)
            
            j += 1
        else:
            i += 1
            batch.append(row["name"])

    # Adding the records of unfinished last batch
    print("Number of chromosomes left: ",len(batch))
    if len(batch) > 0:
        if verbose == True:
                logging.debug(f"\n --------Remaining NCBI records sent to batch_operator--------- \n Contains NCBI records: {batch} ")
        res = ms.batch_operator(batch, faulty, email, api_key, local_storage_path, no_16s)
        dict.update(res[0])
        locus.update(res[1])


    if faulty == []:
        logging.debug("OK! No faulty records from NCBI recorded!")
    else:
        # Retry fetching faulty records t0 subdictionary
        print("")
        print("Retrying fetch for faulty records.")
        batch = []
        no_16s = []
        retry_list = faulty.copy()
        faulty = []
        if platform.system() == 'Windows':
            new_path = local_storage_path + "second\\"
        else:
            new_path = local_storage_path + "second/"
        i = 0
        j = 1
        for x in retry_list:
            if i == 9:
                batch.append(x)
                res = ms.batch_operator(batch, faulty, email, api_key, new_path, no_16s ,verbose )

                print("Batch", j, "of retry")
                print("Number of retrieved retries for batch: ",len(res[0]))
                print("")

                dict.update(res[0])
                locus.update(res[1])
                i = 0
                batch = []
                j += 1
            else:
                i += 1
                batch.append(x)
        if len(batch) > 0:
            res = ms.batch_operator(batch, faulty, email, api_key, new_path, no_16s)
            dict.update(res[0])
            locus.update(res[1])
        logging.warning("Faulty records from NCBI: "+ str(faulty))

    # Recording the records without 16S rRNAs
    if no_16s != []:
            string = ""
            for i in range(len(no_16s)):
                string = string + no_16s[i]
            logging.warning(f" \n -------- No rRNA was found for: -------- {string} \n -------------------------------------------------------------------------------")


    print("-----All chromosmes with corresponding rRNA intervals should be in dict now-----")
    t_fin_2 = time.time()
    print("Total time :",t_fin_2-t_fin_1)    
    print("")
    print("------------------rRNA fetch done----------------")
    return([dict, locus])
