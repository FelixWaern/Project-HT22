import pandas as pd


def fetch_csv_as_df(csv_path):
    """ Function that reads the filtered csv file and creates a DataFrame,
        then calculates ori and ter for all sequences and writes to the DataFrame """

    # Reading filtered csv file
    df = pd.read_csv(csv_path)

    # Finding ori and terminus for all chromosomes
    ter = []
    ori = []

    # Method defined for calculating the origin of replication based on the positive/negative value of shift
    def shift_value(shift, siz):
        if shift > 0:
            origin = shift + 1
        elif shift < 0:
            origin = siz + shift + 1
        else:
            origin = 1
        return origin

    # Method defined for calculate the size after adding the shift value 
    def shift_value_terminus(shift, siz):
        ter = (siz * df.loc[i, "div"]) + shift + 1
        if ter > siz:
            new_ter = ter - siz
            return new_ter
        elif ter < 0:
            new_ter = (siz + ter)
            return new_ter
        else:
            return ter

    # For loop for iterating over all chromosomes
    for i in range(len(df)):
        shift = df.loc[i, "shift"]
        siz = df.loc[i, "siz"]
        origin = shift_value(shift, siz)
        ori.append(origin)
        terminus = shift_value_terminus(shift, siz)
        rounded_terminus = round(terminus)
        ter.append(rounded_terminus)
        
    # Checking the range for ori and ter and assigning those to dataframe
    if any(item <= siz and item > 0 for item in ori):
        print("Ori is within range")
        df["Ori"] = ori
        print("Ori has been added successfully")
    else:
        print("Ori is not within the range")

    if any(item <= siz and item > 0 for item in ter):
        print("Terminus is within range")
        df["Ter"] = ter
        print("Terminus has been added successfully")
    else:
        print("terminus is not within the range")
    
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    return(df)
