# Combining the fetching data & csv_filtered
def get_rRNA_intervals(csv_path, email, api_key, local_storage_path, a_list, verbose=False):
    import time
    import sys
    import logging
    sys.path.insert(0, '/skewDB/')
    from skewDB import fetching_data as fd
    sys.path.insert(0, '/NCBI_DATA_FETCH/')
    from NCBI_DATA_FETCH import main_script as ms
    
    org_df = fd.fetch_csv_as_df(csv_path) 
    #Ta från [22000:23000]
    #test_df = org_df.iloc[900:1000] # [700:800] Did not work :/ 
    test_df = org_df.head(12000)

    if a_list != None:
        df = org_df.loc[org_df['name'].isin(a_list)]
    else:
        df = org_df
    
    i = 0
    j = 1
    batch = []
    faulty = []
    dict = {}
    t_tot = []
    no_16s = []
    t_fin_1 = time.time()
    for index, row in test_df.iterrows():
        if i == 9:
            batch.append(row["name"])
            if verbose == True:
                logging.debug(f"\n --------Batch {j} sent to batch_operator--------- \n Contains NCBI records: {batch} ")
            t0 = time.time()
            res = ms.batch_operator(batch, faulty, email, api_key, local_storage_path, no_16s ,verbose )
            dict.update(res)
            i = 0
            batch = []
            t1 = time.time()
            total = t1-t0
            t_tot.append(total)
            
            print("")
            print("Batch:",j, "done!")
            print("Batch:",j, "took:", total, "seconds!")
            print("Estimated time left: ", ((sum(t_tot)/len(t_tot))*2700)-((sum(t_tot)/len(t_tot))*j) ,"seconds")
            print("Estimated mean total time: ", ((sum(t_tot)/len(t_tot))*2700),"seconds")
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
        dict.update(res)
    
    if no_16s != []:
            string = ""
            for i in range(len(no_16s)):
                string = string + no_16s[i]
            logging.warning(f" \n -------- No rRNA was found for: -------- {string} \n -------------------------------------------------------------------------------")    

    print("-----All chromosmes with corresponding rRNA intervals should be in dict now-----")
    t_fin_2 = time.time()
    print("Total time :",t_fin_2-t_fin_1)

    print("")
    print("Faulty records from NCBI: ", faulty)
    if faulty == []:
        logging.debug("OK! No faulty records from NCBI recorded!")
    else:
        logging.warning("Faulty records from NCBI: "+ str(faulty))
    print("")
    print("------------------rRNA fetch done----------------")
    return(dict)

