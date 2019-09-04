import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import utils as ut


#update to direct towards proper tracks file (mTrackJ format)
#designate metadata file provided in csv
#designate which track is the notochord
files = ["data/170315_tracks",
        #"181113_tracks3",
        ]

meta = ["data/Unfinished_170315_Tracking.csv",
        #"Unfinished_181113_Tracking.csv",
        ]

notochord = [34,
            #5,
            ]

#designate angle and axis of rotation to register tracks anatomically
angle = [60,
        #30,
         ]
         
axis = ['z',
        #'z',
        ]

#tracks this long or shorter will not be included
len_cutoff = 10



#-----------------------------

#script to open file with tracking data, assumes format is similar to output from mTrackJ macro in FIJI.
#track_order is a 1 by x matrix storing track numbers
#track_length is a 1 by x matrix storing number of points in each track
#track_duration is a 1 by x matrix storing [min, max] times analyzed in a given track

track_order, track_lengths, track_duration, noto_index = ut.read_track_data(files[0], notochord[0], len_cutoff, True)


    
#exit()


#-----------------------------


#generates matrix 'embryo' with dimensions (tracks, coordinates(x, y, z), points)

embryo = ut.create_embryo(files[0], track_order, track_lengths, track_duration, noto_index)


#-----------------------------

#creates new matrix 'calibrated_embryo' in which all values are registered to the notochord track
#duplicates first and last points to fill duration of dataset in case track is short - wil keep t lined up with 'points' dimension

calibrated_embryo = ut.calibrate_to_noto(embryo, noto_index, track_duration)

#in case of missing coordinates from a timepoint, prints missing timepoint
#takes coordinates at previous and subsequent timepoints and interpolates to generate missing point coordinates

calibrated_embryo = ut.interpolate_missing_timepoints(calibrated_embryo)

#print(calibrated_embryo[1, :, :])

#-----------------------------


rotated_embryo = ut.rotate_embryo(calibrated_embryo, angle[0], axis[0])
#rotated_embryo = rotate(rot_matrix(30, 'x'), rotated_embryo)

#-----------------------------


'''
imports metadata from csv provided in 'metadata' variable at top
requires columns 'Track Number', 'Sister Cell Track', 'Identity of Tracked Cell', 'Cell Division Timing'
creates list of birthplaces and color matrices for plotting tracks/points from lineage/fate

'''

track_label = ut.import_meta(meta[0], 'Track Number')
sister = ut.import_meta(meta[0], 'Sister Cell Track')

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


cell_fate = ut.import_meta(meta[0], 'Identity of Tracked Cell')

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

progenitor = ut.import_meta(meta[0], 'Cell Division Timing')
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


ut.plot_tracks(rotated_embryo, files[0]+'_rotated', smoothing=3)
#plot_tracks(rotated_embryo,'_YZ', view=[1,2])
#plot_tracks(rotated_embryo,'_rot_supersmooth2',smoothing=3)
#plot_tracks(rotated_embryo,'_rot_supersmooth2',smoothing=3,color=[1,100,1,5,1,5,1,5,1,5,1,5,1,1])


#plot_tracks(rotated_embryo,'_rot_lineage',smoothing=3,color=lineage) 
#plot_tracks(rotated_embryo,'_rot_fate_NEW',smoothing=3,color=fate)
#plot_birthplaces(birthplaces,'birthplaces')



#MAKE SURE DURATION IS RIGHT - do times go into proper slots or does everything start at first embryo entry?
#they go into proper slots but why???
#TO FIX: last color in array does not show up
#clean up error messages on import_meta 
#color birthplaces by fate or division type (do same way by making color matrices)
#put projections option into plot_tracks > DONE
#combine tracks from multiple embryos


exit()


