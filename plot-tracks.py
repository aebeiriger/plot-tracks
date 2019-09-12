import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import utils as ut


#update to direct towards proper tracks file (mTrackJ format)
#designate metadata file provided in csv
#designate which track is the notochord
files = [#"data/170315_tracks",
        "data/181113_tracks3",
        ]

metas = [#"data/Unfinished_170315_Tracking.csv",
        "data/Unfinished_181113_Tracking.csv",
        ]

notochords = [#34,
            5,
            ]

#designate angle and axis of rotation to register tracks anatomically
angles = [#60,
        30,
         ]
         
axes = [#'z',
        'z',
        ]

#tracks this long or shorter will not be included
len_cutoff = 10



#-----------------------------

#script to open file with tracking data, assumes format is similar to output from mTrackJ macro in FIJI.
#track_order is a 1 by x matrix storing track numbers
#track_length is a 1 by x matrix storing number of points in each track
#track_duration is a 1 by x matrix storing [min, max] times analyzed in a given track

for file, meta, notochord, angle, axis in zip(files, metas, notochords, angles, axes):
    track_order, track_lengths, track_duration, noto_index = ut.read_track_data(file, notochord, len_cutoff, True)

    #generates matrix 'embryo' with dimensions (tracks, coordinates(x, y, z), points)
    embryo = ut.create_embryo(file, track_order, track_lengths, track_duration, noto_index)

    #creates new matrix 'calibrated_embryo' in which all values are registered to the notochord track
    #duplicates first and last points to fill duration of dataset in case track is short - wil keep t lined up with 'points' dimension
    calibrated_embryo = ut.calibrate_to_noto(embryo, noto_index, track_duration)

    #in case of missing coordinates from a timepoint, prints missing timepoint
    #takes coordinates at previous and subsequent timepoints and interpolates to generate missing point coordinates
    calibrated_embryo = ut.interpolate_missing_timepoints(calibrated_embryo)

    rotated_embryo = ut.rotate_embryo(calibrated_embryo, angle, axis)

    track_label = ut.import_meta(meta, 'Track Number')

    fate = ut.parse_metadata(meta, file, track_order, mode='fate')

    lineage = ut.parse_metadata(meta, file, track_order, mode='lineage')

    birthplaces = ut.get_birthplaces(meta, file, rotated_embryo, track_order, track_label, track_duration)

    ut.plot_tracks(rotated_embryo, smoothing=3, color=lineage)
    plt.savefig('compiled_tracks')
    plt.clf()
    #plot_tracks(rotated_embryo,'_YZ', view=[1,2])
    #plot_tracks(rotated_embryo,'_rot_supersmooth2',smoothing=3)
    #plot_tracks(rotated_embryo,'_rot_supersmooth2',smoothing=3,color=[1,100,1,5,1,5,1,5,1,5,1,5,1,1]
    #plot_tracks(rotated_embryo,'_rot_lineage',smoothing=3,color=lineage) 
    #plot_tracks(rotated_embryo,'_rot_fate_NEW',smoothing=3,color=fate)
    ut.plot_birthplaces(birthplaces)
    plt.savefig('compiled_birthplaces')
    plt.clf()



#MAKE SURE DURATION IS RIGHT - do times go into proper slots or does everything start at first embryo entry?
#they go into proper slots but why???
#TO FIX: last color in array does not show up
#clean up error messages on import_meta 
#color birthplaces by fate or division type (do same way by making color matrices)
#put projections option into plot_tracks > DONE
#combine tracks from multiple embryos


exit()


