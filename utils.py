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


def rotate_embryo(tracks_matrix, angle, axis):
    '''
    Function to rotate tracks.
    Starts with initial track points and rotates by a defined angle, around a defined axis.

    tracks_matrix: matrix to rotate, likely matrix 'embryo' (tracks, coordinates(x, y, z), points)
    angle: angle to rotate by in degrees (float)
    axis: axis to rotate around, one of x, y, or z
    
    output is matrix 'rotated_embryo' (tracks, coordinates(x, y, z), points)

    '''

    rot_matrix = _generate_rotation_matrix(angle, axis)
    rotated_embryo = np.matmul(rot_matrix, tracks_matrix)
    return rotated_embryo


def _generate_rotation_matrix(angle, axis):
    rot_matrix = np.eye(3) 
    if axis == 'z':
        rot_matrix[0][0] = np.cos(np.radians(angle))
        rot_matrix[1][0] = np.sin(np.radians(angle))
        rot_matrix[1][1] = np.cos(np.radians(angle))
        rot_matrix[0][1] = -(np.sin(np.radians(angle)))
    elif axis == 'y':
        rot_matrix[0][0] = np.cos(np.radians(angle))
        rot_matrix[0][2] = np.sin(np.radians(angle))
        rot_matrix[2][2] = np.cos(np.radians(angle))
        rot_matrix[2][0] = -(np.sin(np.radians(angle)))
    else:
        rot_matrix[1][1] = np.cos(np.radians(angle))
        rot_matrix[1][2] = np.sin(np.radians(angle))
        rot_matrix[2][2] = np.cos(np.radians(angle))
        rot_matrix[1][2] = -(np.sin(np.radians(angle)))
    return rot_matrix
