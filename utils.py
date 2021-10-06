import os.path as op
import subprocess
from matplotlib import pyplot as plt
import nibabel
from nibabel.streamlines.tractogram import Tractogram
import numpy as np


def load_tracks(path_track, out_track):
    return nibabel.streamlines.tck.TckFile.load(op.join(path_track, out_track + '.tck'))


def save_tracks(streamlines, tracks, out_track_trim, path_track):
    # tracks: original track file that will be used as a reference
    if not ('.tck' in out_track_trim):
        out_track_trim = out_track_trim + '.tck'
    new_header = tracks.header.copy()
    # change relevant information
    new_header['nb_streamlines'] = len(streamlines)
    new_header['count'] = str(len(streamlines))
    new_header['command_history'] = new_header['command_history'] + ' streamline trimming'
    new_tractogram = Tractogram(streamlines, affine_to_rasmm=np.eye(4))
    new_track = nibabel.streamlines.tck.TckFile(new_tractogram, new_header)
    new_track.save(op.join(path_track, out_track_trim))


def merge_tracks(path_track,track_name):
    cmd = ['tckedit', '*_trimmed.tck',
           op.join(path_track, track_name+'_trimmed.tck')]
    subprocess.run(cmd)


def filter_streamlines(path_track, track_name, assignment_name, node1=1,node2=2):
    # Extract the streamlines connecting node 'node' to all other nodes
    # in the parcellation, with one track file for each edge:
    out_track = 'edge_'+str(node1)+'_'+str(node2)+'.tck'
    cmd = ['connectome2tck',
           op.join(path_track, track_name + '.tck'),
           op.join(path_track, assignment_name),
           op.join(path_track,out_track),
           '-nodes', str(node1)+','+str(node2),
           '-exclusive','-files','single']
    subprocess.call(cmd)
    return out_track

#connectome2tck tracks.tck assignments.txt tracks_1_2.tck -nodes 1,2 -exclusive -files single
def plot_streamline_outliers(streamlines, outliers):
    # outliers will be plotted in different color
    # extract streamlines from tractogram
    def check_type(T):
        if isinstance(T, nibabel.streamlines.tractogram.Tractogram):
            temp = [tt.streamline for tt in T]
            return temp
        elif isinstance(T, nibabel.streamlines.tck.TckFile):
            temp = T.streamlines
            return temp
        else:
            return T

    streamlines = check_type(streamlines)
    outliers = check_type(outliers)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for i, str1 in enumerate(streamlines):
        if i == 0:
            label = 'streamline'
        else:
            label = None
        ax.plot(str1[:, 0], str1[:, 1], str1[:, 2],
                c='blue', alpha=0.8, label=label)
    for i, str1 in enumerate(outliers):
        if i == 0:
            label = 'outlier'
        else:
            label = None
        ax.plot(str1[:, 0], str1[:, 1], str1[:, 2],
                c='pink', alpha=0.8, label=label)
    fig.tight_layout()
    ax.legend()
    plt.show()
