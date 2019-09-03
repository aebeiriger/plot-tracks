import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_from_metadata import import_meta 


#update to direct towards proper tracks file (mTrackJ format)
#designate metadata file provided in csv
#designate which track is the notocord
file = "170315_tracks"
meta = "Unfinished_170315_Tracking.csv"
notocord = 34

#designate angle and axis of rotation to register tracks anatomically
angle = 60
axis = 'z'

#tracks this long or shorter will not be included
len_cutoff = 10



#function to rotate tracks, starts with initial track points and rotates by angle defined above

def rotate(rot_matrix, tracks_matrix):
    '''
    input angle in degrees (float) and initial matrix 'embryo' (tracks, coordinates(x, y, z), points)
    output is matrix 'rotated_embryo' (tracks, coordinates(x, y, z), points)

    '''

    rotated_embryo = np.matmul(rot_matrix, tracks_matrix)
    return rotated_embryo

#function to define rotation matrix, takes inputs of angle and axis of rotation

def rot_matrix(angle, axis):
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

#function to plot tracks, with arguments to specify color schemes, create plots with equal axes aka 'square',
#and smooth tracks by averaging over a designated number of timepoints

def plot_tracks(tracks_matrix,name,view=[0, 1],color=None,square=True,smoothing=1):
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
    plt.savefig(file+name)
    plt.clf()

#function to plot cell birthplaces, with arguments to specify color schemes, create plots with equal axes aka 'square',
#and project in dorsal, transverse, or lateral view

def plot_birthplaces(birthplaces,name,projection=[0,1],color=None,square=True):
    nprogenitors = len(birthplaces)
    if color is None:
        color = range(nprogenitors)
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
    plt.savefig(file+name)
    plt.clf()

#potato = rot_matrix(30, 'z')
#print(potato)

#exit()

#-----------------------------

#script to open file with tracking data, assumes format is similar to output from mTrackJ macro in FIJI.
#track_order is a 1 by x matrix storing track numbers
#track_length is a 1 by x matrix storing number of points in each track
#track_duration is a 1 by x matrix storing [min, max] times analyzed in a given track

tracks=open(file)
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

print('track_order',track_order)
print('track_lengths',track_lengths)
print('track_duration',track_duration)
print('nr_tracks',len(track_order))

noto_index = track_order.index(notocord)
print('noto_index',noto_index)



tracks.close()
#exit()


#-----------------------------


#generates matrix 'embryo' with dimensions (tracks, coordinates(x, y, z), points)

embryo = np.zeros((len(track_order), 3, track_duration[noto_index][1]), dtype=int)
tracks = open(file)
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

#print('pre-cal',embryo[1])


#-----------------------------

#creates new matrix 'calibrated_embryo' in which all values are registered to the notocord track
#duplicates first and last points to fill duration of dataset in case track is short - wil keep t lined up with 'points' dimension

calibrated_embryo = embryo[:, :, :] - embryo[noto_index, :, :]
for i in range(0, calibrated_embryo.shape[0]):
    start = max(track_duration[i][0], track_duration[noto_index][0])
    end = min(track_duration[i][1], track_duration[noto_index][1])
    calibrated_embryo[i, :, 0:start] = calibrated_embryo[i, :, start:start+1]
    calibrated_embryo[i, :, end:] = calibrated_embryo[i, :, end-1:end]

#in case of missing coordinates from a timepoint, prints missing timepoint
#takes coordinates at previous and subsequent timepoints and interpolates to generate missing point coordinates

for i in range(0, calibrated_embryo.shape[0]):
    for j in range(0, calibrated_embryo.shape[1]):
        for k in range(1, calibrated_embryo.shape[2]):
            if (abs(calibrated_embryo[i, j, k] - calibrated_embryo[i, j, k-1]) > 40):
                calibrated_embryo[i, j, k] = (calibrated_embryo[i, j, k-1] + calibrated_embryo[i, j, k+1])/2
                print(i, j, k)

#print('post-cal',calibrated_embryo[1])
print('cal_noto', calibrated_embryo[noto_index, :, :])

#-----------------------------


rotated_embryo = rotate(rot_matrix(angle, axis), calibrated_embryo)
#rotated_embryo = rotate(rot_matrix(30, 'x'), rotated_embryo)

#-----------------------------


'''
imports metadata from csv provided in 'metadata' variable at top
requires columns 'Track Number', 'Sister Cell Track', 'Identity of Tracked Cell', 'Cell Division Timing'
creates list of birthplaces and color matrices for plotting tracks/points from lineage/fate

'''

track_label = import_meta(meta, 'Track Number')
sister = import_meta(meta, 'Sister Cell Track')

if len(track_order) != len(track_label):
    print('metadata and track data do not match')
    exit(1)

lineage = [-1 for item in track_order]
k = 0

'''
track_order and track_label might not follow the same order
lineage is built in the same sequence as track_order
track_label is only used as a reference to find sister cells

'''

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
        
'''        
converts cell fate data to integers
creates color matrix 'fate' based on cell identities

'''        

def str_to_num(string):
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


cell_fate = import_meta(meta, 'Identity of Tracked Cell')

if len(track_order) != len(cell_fate):
    print('fate entry missing')
    exit(1)

fate = []
for track in track_order:
    j = track_label.index(track)
    fate.append(cell_fate[j])

fate = [str_to_num(a) for a in fate]


'''
builds a list of cell birthplace coordinates [x,y,z]
references point "t" in rotated_embryo matrix to get coordinates
if no birthplace is found, [-1,-1,-1] added so list is compatible with color matrices

'''

progenitor = import_meta(meta, 'Cell Division Timing')
birthplaces = []

for i, time in enumerate(progenitor):
    try:
        int(time)
    except:
        print(i, 'division time is not a number')
    if time != -1:
        t = time
        j = track_label[i]
        j = track_order.index(j)
        if track_duration[j][0] < t:
            birthplaces.append(list(rotated_embryo[j, :, t]))
        else:
            birthplaces.append([-1,-1,-1])    
    elif time == -1:
        birthplaces.append([-1,-1,-1])

#print(birthplaces)
        

#--------------------------------------------


plot_tracks(rotated_embryo,'_rot2')
plot_tracks(rotated_embryo,'_YZ', view=[1,2])
plot_tracks(rotated_embryo,'_rot_supersmooth2',smoothing=3)
plot_tracks(rotated_embryo,'_rot_supersmooth2',smoothing=3,color=[1,100,1,5,1,5,1,5,1,5,1,5,1,1])


plot_tracks(rotated_embryo,'_rot_lineage',smoothing=3,color=lineage) 
plot_tracks(rotated_embryo,'_rot_fate_NEW',smoothing=3,color=fate)
plot_birthplaces(birthplaces,'birthplaces')


#TO FIX: last color in array does not show up
#clean up error messages on import_meta 
#color birthplaces by fate or division type (do same way by making color matrices)
#put projections option into plot_tracks > DONE
#combine tracks from multiple embryos


exit()


