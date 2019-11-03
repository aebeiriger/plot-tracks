import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



#----------------------------------------------------------------------------------
#plotting functions


def plot_tracks(tracks_matrix, projection=[0, 1], color=None, square=True, smoothing=1, scaling=[1,1,1], limits=None, birth_times=None, lineages=None):
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
        if birth_times is None:
            plt.plot(tracks_matrix[i, projection[0], :].T*scaling[projection[0]],
                 tracks_matrix[i, projection[1], :].T*scaling[projection[1]],c=color[i])
        elif birth_times is not None:
                if birth_times[i] != -1:
                    plt.plot(tracks_matrix[i, projection[0], birth_times[i]-2].T*scaling[projection[0]],
                         tracks_matrix[i, projection[1], birth_times[i]-2].T*scaling[projection[1]],c=color[i],marker='o',markersize=8)

             
    if square:  
        plt.gca().set_aspect('equal')

    if limits:
        plt.xlim(limits[0][0],limits[0][1])
        plt.ylim(limits[1][0],limits[1][1])


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


def calibrate_to_noto(embryo, noto_index, track_duration, normalize_to_noto=True, flip_perspective=True):
    if flip_perspective:
        embryo[:,1,:] = -embryo[:, 1, :]
        
    if normalize_to_noto:
        calibrated_embryo = embryo[:, :, :] - embryo[noto_index, :, :]
    else:
        calibrated_embryo = embryo
    

    
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


def rotate_embryo(tracks_matrix, angles, axes):
    '''
    Function to rotate tracks.
    Starts with initial track points and rotates by a defined angle, around a defined axis.

    tracks_matrix: matrix to rotate, likely matrix 'embryo' (tracks, coordinates(x, y, z), points)
    angle: angle to rotate by in degrees (float)
    axis: axis to rotate around, one of x, y, or z
    
    output is matrix 'rotated_embryo' (tracks, coordinates(x, y, z), points)

    '''
    for angle, axis in zip(angles, axes): 
        rot_matrix = _generate_rotation_matrix(angle, axis)
        tracks_matrix = np.matmul(rot_matrix, tracks_matrix)
    return tracks_matrix


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

def import_metadata(filename, index='Track Number'):
    metadata = pd.read_csv(filename)
    metadata = metadata.set_index(index)
    return metadata
    
def get_column_from_metadata(metadata, column_name, track_order):
    column_to_return = []
    for track in track_order:
        column_to_return.append(metadata.loc[track,column_name])
    
    return column_to_return


def get_lineages(metadata, track_order):
    '''
    def recursive_lineage(lineages, track_order, sister_tracks, index, nr_to_set):
        lineages[index] = nr_to_set
        if sister_tracks[index] == -1:
            return
        
        else:
            recursive_lineage(lineages, track_order, sister_tracks, track_order.index(sister_tracks[index]), nr_to_set)
            
    lineages = [-1 for item in track_order]
    sister_tracks = get_column_from_metadata(metadata, 'Sister Cell Track', track_order)
    for i in range(len(track_order)):
        recursive_lineage(lineages, track_order, sister_tracks, i, i)
        
    print(lineages)
    return(lineages)
   
    '''         
    lineages = [-1 for item in track_order]
    nr_lineages = -1
    for i, track in enumerate(track_order):
        if lineages[i] == -1:
            nr_lineages += 1
            lineages[i] = nr_lineages
        sister_cell = metadata.loc[track, 'Sister Cell Track']
        if sister_cell != -1:
            lineages[track_order.index(sister_cell)] = lineages[i]
    return lineages
    


def _fate_to_num(fate, fate_to_num_dict):
    '''        
    converts cell fate data to integers
    creates color matrix 'fate' based on cell identities

    '''
    return [fate_to_num_dict[a] for a in fate]


def average_birthplaces(tracks_matrix, lineages, birth_times):
    averaged_matrix = tracks_matrix
    for i in range(len(lineages)):
        for j in range(len(lineages)):
            if lineages[i] == lineages[j] and birth_times[i] == birth_times[j]:
                averaged_matrix[i, :, birth_times[i]-2] = (tracks_matrix[i, :, birth_times[i]-2] + tracks_matrix[j, :, birth_times[i]-2])/2
                averaged_matrix[j, :, birth_times[i]-2] = averaged_matrix[i, :, birth_times[i]-2]

    return averaged_matrix
    