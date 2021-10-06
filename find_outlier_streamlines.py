import numpy as np
from interpcurve import interpcurve


def find_outlier_streamlines(tracks, theta=5, cutoff_member=5):
    # Find outlier streamlines
    # :param: tracks
    # :param: theta        : distance for a streamline to be assigned to a cluster centroid
    # :param: cutoff_member: minimum number of members a centroid should have not to be classified as an outlier
    # :return:

    num_stream = len(tracks.tractogram)

    # find the median length of the track
    stream_length = [np.sqrt(np.sum(np.diff(tt.streamline, axis=0) ** 2,
                                    axis=1)).sum()
                     for tt in tracks.tractogram]
#    track_length = np.floor(np.median(stream_length))
    #use a predetermined track length
    track_length = np.min((12,np.floor(np.median(stream_length))))

    # resample streamlines so that they have equal number of points
    # first cluster is equal to the first tract
    str1 = tracks.tractogram[0].streamline
    str1_r = interpcurve(int(track_length), str1[:, 0], str1[:, 1], str1[:, 2])
    cluster_centroids = [str1_r]
    assign_track = np.empty((num_stream, 1))
    assign_track[0] = 0
    nc = 1
    for ii in range(1, num_stream):
        streamline = tracks.tractogram[ii].streamline
        # resample streamline
        streamline_r = interpcurve(int(track_length),
                                   streamline[:, 0],
                                   streamline[:, 1],
                                   streamline[:, 2])
        # calculate distance to cluster centre
        d = np.zeros((nc, 1))
        for ic in range(nc):
            d[ic] = np.sqrt(np.sum((streamline_r - cluster_centroids[ic]) ** 2,
                                   axis=1)).sum()
        d /= track_length
        # assign it to the min and < 8, if not create a new cluster
        if d[d < theta].shape[0] == 0:
            # new cluster centre
            nc = nc+1
            assign_track[ii] = nc-1
            cluster_centroids.append(streamline_r)
        else:
            I = np.argmin(d)
            # assign to closest cluster
            nmemcc = np.sum(assign_track == I)
            assign_track[ii] = I
            # re-calculate new cluster centre by re-weighting
            cluster_centroids[I] = (nmemcc * cluster_centroids[I] + streamline_r) / (nmemcc + 1)

    # number of elements in each cluster
    num_member = np.zeros(nc)
    for ic in range(nc):
        num_member[ic] = np.sum(assign_track == ic)
        # eliminate cluster with a few members
        if num_member[ic] <= cutoff_member:
            assign_track[assign_track == ic] = np.nan
            #num_member[ic] = 0

    # remove outlier tracks
    tracks_clean = tracks.tractogram[~np.isnan(assign_track).flatten()].streamlines
    index_tracks_clean = 1 * (~np.isnan(assign_track).flatten())
    outlier_streamlines = tracks.tractogram[np.isnan(assign_track).flatten()].streamlines

    return tracks_clean, index_tracks_clean, outlier_streamlines, cluster_centroids, num_member
