import numpy as np

def interpcurve(N,pX,pY,pZ):
    # Interpolate points along a curve in 3 dimensions
    # This code is based on
    # John D 'Errico (2021). interparc
    # (https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc),
    # MATLAB Central File Exchange. Retrieved June 10, 2021.
    # and python version of the same function
    # https://stackoverflow.com/questions/18244305/how-to-redistribute-points-evenly-over-a-curve
    # This version is implemented for linear and 3-dimensional interpolation
    # N: number of points that curve will be uniformly resampled
    # pX:
    # pY:
    # equally spaced in arclength
    N = np.transpose(np.linspace(0,1,N))

    #how many points will be uniformly interpolated?
    nt = N.size
    ndim = 3
    #number of points on the curve
    n = pX.size
    pxy = np.array((pX, pY, pZ)).T
#    p1=pxy[0,:]
#    pend=pxy[-1,:]
#    last_segment= np.linalg.norm(np.subtract(p1,pend))
#    epsilon= 10*np.finfo(float).eps

    #IF the two end points are not close enough lets close the curve
#    if last_segment > epsilon*np.linalg.norm(np.amax(abs(pxy),axis=0)):
#        pxy=np.vstack((pxy,p1))
#        nt = nt + 1
#    else:
#        print('Contour already closed')

    pt=np.zeros((nt,ndim))

    #Compute the chordal arclength of each segment.
    chordlen = (np.sum(np.diff(pxy,axis=0)**2,axis=1))**(1/2)
    #Normalize the arclengths to a unit total
    chordlen = chordlen/np.sum(chordlen)
    #cumulative arclength
    cumarc = np.append(0,np.cumsum(chordlen))
    # which interval did each point fall in, in
    # terms of t?
    tbins= np.digitize(N,cumarc) # bin index in which each N is in

    #catch any problems at the ends
    tbins[np.where(np.bitwise_or(tbins<=0, N<=0))] = 1
    tbins[np.where(np.bitwise_or(tbins >= n , N >= 1))] = n - 1
    # interpolate
    s = np.divide((N - cumarc[tbins-1]),chordlen[tbins-1])
    pt = pxy[tbins-1,:] + np.multiply((pxy[tbins,:] - pxy[tbins-1,:]),(np.vstack([s]*ndim)).T)

    return pt

