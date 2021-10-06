import os.path as op
import argparse, os
from utils import filter_streamlines, load_tracks, plot_streamline_outliers, save_tracks, merge_tracks
from find_outlier_streamlines import find_outlier_streamlines
from config import path_subj, path_track
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('subject', metavar='1###', help='The subject to process')
parser.add_argument('assignment', metavar='1###', help='Assignment file for streamlines (start/end ROIs)')
args = parser.parse_args()
subject = args.subject
assignment_name=args.assignment

print('Trimming tracks for subject:', subject, ' by ', assignment_name)

# subject = '1001'
path_track_subj = path_track(subject)
track_name = 'tracks_10m_wb'
assignment_name = 'assignments_SC.csv'

num_reg = 379
suffix = '_trim'

out_track_sf = 'edge_'

cutoff_member = 3
theta = 20# mm
total_num_tracks = 0
num_outlier_tracks = 0
num_trim_tracks = 0
cluster_size_track = []
for node1 in range(1, num_reg + 1):  # node number starts from 1
    for node2 in range(node1 + 1, num_reg + 1):  # account for symmetry
        out_track = "{}{}-{}".format(out_track_sf, node1, node2)
        print(out_track)
        # load with nibabel
        tracks = load_tracks(path_track_subj, out_track)
        # if there are no streamlines between these nodes, skip
        total_num_tracks += len(tracks.streamlines)
        if len(tracks.streamlines) < cutoff_member:
            num_outlier_tracks += len(tracks.streamlines)
            os.remove(op.join(path_track_subj, out_track + '.tck'))
            continue
        # trim streamlines
        streamlines_clean, index_tracks_clean, \
        streamlines_outlier, cluster_centroids, cluster_size \
            = find_outlier_streamlines(tracks, theta=theta, cutoff_member=cutoff_member)
        # if all streamlines are regarded as outlier
        #   due to less number of streamlines (5), skip
        num_outlier_tracks += len(streamlines_outlier)
        # add edge name
        cluster_size_track.append(np.append((node1,node2),cluster_size))
        os.remove(op.join(path_track_subj, out_track + '.tck'))
        if len(streamlines_clean) == 0:
            continue
        # To visualize:
        # plot_streamline_outliers(streamlines_clean[:10], streamlines_outlier[:10])
        # save trimmed streamlines
        out_track_trim = out_track + suffix + '.tck'
        save_tracks(streamlines_clean, tracks, out_track_trim, path_track_subj)
        num_trim_tracks += len(streamlines_clean)

# write stats in report file
file = open(op.join(path_track_subj, 'report_' + suffix + '.txt'), 'w')
file.write('Total number of valid tracks: {}\n'.format(total_num_tracks))
file.write('Number of outlier tracks: {}\n'.format(num_outlier_tracks))
file.write('Number of cleaned tracks: {}\n'.format(num_trim_tracks))
file.write('Percentage of outlier tracks: {}\n'.format(100*num_outlier_tracks/total_num_tracks))
file.close()
# write number of members in each cluster per edge
# First two columns are node1 and node2
clusters = pd.DataFrame(cluster_size_track)
clusters.to_csv(op.join(path_track_subj, 'report_cluster_size' + suffix + '.csv'))
