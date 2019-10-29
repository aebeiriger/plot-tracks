import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



#----------------------------------------------------------------------------------
#plotting functions


def plot_tracks(tracks_matrix, view=[0, 1], color=None, square=True, smoothing=1):
    '''
    Function to plot tracks.
    with arguments to specify color schemes, create plots with equal axes aka 'square',
    and smooth tracks by averaging over a designated number of timepoints

    '''
    ntracks = tracks_matrix.shape[0]
    if color is None:
        color = range(ntracks)
    elif not len(color) == ntracks:
        print('wrong size color')
        exit(1)

    max_color = max(color)
    color = [a/max_color for a in color]
    color = plt.cm.rainbow(color)
    for i in range(smoothing):
        tracks_matrix[:,:,1:-1] = (tracks_matrix[:,:,0:-2] + tracks_matrix[:,:,1:-1] + tracks_matrix[:,:,2:])/3

    for i in range(ntracks):
        plt.plot(tracks_matrix[i, view[0], :].T,
             tracks_matrix[i, view[1], :].T,c=color[i])
    if square:  
        plt.axis('equal')


def plot_birthplaces(birthplaces, projection=[0,1], color=None, square=True):
    '''
    function to plot cell birthplaces, with arguments to specify color schemes, create plots with equal axes aka 'square',
    and project in dorsal, transverse, or lateral view

    '''
    nprogenitors = len(birthplaces)
    if color is None:
        color = list(range(nprogenitors))
    elif not len(color) == nprogenitors:
        print('wrong size color')
        exit(1)

    max_color = max(color)
    color = [a/max_color for a in color]
    color = plt.cm.rainbow(color)
    
    assert len(projection) == 2

    for i in range(nprogenitors):
        if birthplaces[i] != [-1,-1,-1]:
            plt.plot(birthplaces[i][projection[0]],
                     birthplaces[i][projection[1]],c=color[i],marker='o',markersize=12)
    if square:  
        plt.axis('equal')


#----------------------------------------------------------------------------------
#generate and calibrate embryo matrix from tracked points

    
def read_track_data(filename, notochord, len_cutoff, verbose=False):
    tracks=open(filename)
    track_order = []
    track_lengths = []
    track_duration = []
    count = 0

    for line in tracks:
        #print(line)
        points = line.split()
        #print(points)
        if points[0] == "Track":
            if count > len_cutoff:
                track_order.append(track_num)
                track_lengths.append(count)
                track_duration.append([min_time, max_time])
            track_num = int(points[1])
            count = 0
            max_time = 0
            min_time = 9999
        elif points[0] == "Point":
            count = count + 1
            if int(float(points[5])) > max_time:
                max_time = int(float(points[5]))
            if int(float(points[5])) < min_time:
                min_time = int(float(points[5]))

    track_order.append(track_num)
    track_lengths.append(count)
    track_duration.append([min_time, max_time])

    noto_index = track_order.index(notochord)

    if verbose:

        print('filename', filename)
        print('track_order',track_order)
        print('track_lengths',track_lengths)
        print('track_duration',track_duration)
        print('nr_tracks',len(track_order))

        print('noto_index',noto_index)

    tracks.close()
    return track_order, track_lengths, track_duration, noto_index


def create_embryo(filename, track_order, track_lengths, track_duration, noto_index):
    embryo = np.zeros((len(track_order), 3, track_duration[noto_index][1]), dtype=int)
    tracks = open(filename)
    t = -1

    for line in tracks:
        #print(line)
        points = line.split()
        #print(points)
        if points[0] == "Track":
            if int(points[1]) in track_order:
                append = True
                t = t + 1
            else:
                append = False
        elif points[0] == "Point":
            if append:
                p = int(float(points[5])) - 1
                embryo[t, :, p] = list(map(int,list((map(float,points[2:5])))))

    return embryo


def calibrate_to_noto(embryo, noto_index, track_duration):
    calibrated_embryo = embryo[:, :, :] - embryo[noto_index, :, :]
    for i in range(0, calibrated_embryo.shape[0]):
        start = max(track_duration[i][0], track_duration[noto_index][0])
        end = min(track_duration[i][1], track_duration[noto_index][1])
        calibrated_embryo[i, :, 0:start] = calibrated_embryo[i, :, start:start+1]
        calibrated_embryo[i, :, end:] = calibrated_embryo[i, :, end-1:end]

    return calibrated_embryo


def interpolate_missing_timepoints(calibrated_embryo, verbose=False):
    for i in range(0, calibrated_embryo.shape[0]):
        for j in range(0, calibrated_embryo.shape[1]):
            for k in range(1, calibrated_embryo.shape[2]):
                if (abs(calibrated_embryo[i, j, k] - calibrated_embryo[i, j, k-1]) > 40):
                    calibrated_embryo[i, j, k] = (calibrated_embryo[i, j, k-1] + calibrated_embryo[i, j, k+1])/2
                    if verbose:
                        print(i, j, k)

    if verbose:
        print('cal_noto', calibrated_embryo[noto_index, :, :])

    return calibrated_embryo


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



#----------------------------------------------------------------------------------
#metadata functions

def import_meta(filename, columnname):
    '''
    filename (str) - the metadata file to import from
    columnname (str) - the column name to import data from

    returns: a list containing entries of specified column from the metadata file
    '''
    df = pd.read_csv(filename)
    return list(df[columnname].values)


def parse_metadata(meta, filename, track_order, mode=None):
    '''
    imports metadata from csv provided in 'meta' variable at top
    requires columns 'Track Number', 'Sister Cell Track', 'Identity of Tracked Cell', 'Cell Division Timing'
    creates matrices for color-coding plots by lineage or fate

    '''
    track_label = import_meta(meta, 'Track Number')

    if len(track_order) != len(track_label):
        print('metadata and track data do not match')
        exit(1)

    '''
    track_order (raw data) and track_label (meta data) might not follow the same order
    'lineage' and 'fate' are built in the same sequence as track_order
    track_label is only used as a reference

    '''

    if mode is 'lineage':
        sister = import_meta(meta, 'Sister Cell Track')
        lineage = [-1 for item in track_order]
        k = 0

        for i in range(len(lineage)):
            if lineage[i] == -1:
                lineage[i] = k
                j = track_order[i]
                j = track_label.index(j)
                if sister[j] != -1:
                    l = sister[j]
                    l = track_order.index(l)
                    lineage[l] = k
                k=k+1

        return lineage

    elif mode is 'fate':
        cell_fate = import_meta(meta, 'Identity of Tracked Cell')

        if len(track_order) != len(cell_fate):
            print('fate entry missing')
            exit(1)

        fate = []
        for track in track_order:
            j = track_label.index(track)
            fate.append(cell_fate[j])

        fate = [_fate_by_num(a) for a in fate]

        return fate

    return track_label


def _fate_by_num(string):
    '''        
    converts cell fate data to integers
    creates color matrix 'fate' based on cell identities

    '''
    if string == 'NOTO':
        return 0
    elif string == 'FBMN/REN':
        return 1
    elif string == 'CEN':
        return 2
    elif string == 'Isl1(-)':
        return 3
    #Error handling
    else:
        return -1
        print('incorrect fate')


def get_birthplaces(meta, filename, tracks_matrix, track_order, track_label, track_duration):
    '''
    builds a list of cell birthplace coordinates [x,y,z]
    references point "t" in rotated_embryo matrix to get coordinates
    if no birthplace is found, [-1,-1,-1] added so list is compatible with color matrices

    '''
    
    progenitor = import_meta(meta, 'Cell Division Timing')
    birthplaces = []

    for i, time in enumerate(progenitor):
        try:
            time = int(time)
        except:
            print(i, 'division time is not a number')
        if time != -1:
            t = time
            j = track_label[i]
            j = track_order.index(j)
            if track_duration[j][0] < t:
                birthplaces.append(list(tracks_matrix[j, :, t]))
            else:
                pass
                #birthplaces.append([-1,-1,-1])    
        elif time == -1:
            pass
            #birthplaces.append([-1,-1,-1])

    return birthplaces



