import pandas as pd
import numpy as np
#import babys_first_python.py
#print(track_order)

meta = "Unfinished_170315_Tracking.csv"

#returns specified column values as a list from a metadata csv

def import_meta(filename, columnname):
    '''
    filename (str) - the metadata file to import from
    columnname (str) - the column name to import data from

    returns: a list containing entries of specified column from the metadata file
    '''
    df = pd.read_csv(filename)
    return list(df[columnname].values)



