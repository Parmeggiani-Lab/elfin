#!/usr/bin/env python3

"""
Calculate RMSD between two XYZ files

by: Jimmy Charnley Kromann <jimmy@charnley.dk> and Lars Andersen Bratholm <larsbratholm@gmail.com>
project: https://github.com/charnley/rmsd
license: https://github.com/charnley/rmsd/blob/master/LICENSE

"""

import numpy as np
import re


def kabsch_rmsd(P, Q):
    """Rotate matrix P unto Q and calculate the RMSD."""
    P = kabsch_rotate(P, Q)
    return rmsd(P, Q)


def kabsch_rotate(P, Q):
    """Rotate matrix P unto matrix Q using Kabsch algorithm."""
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P


def run_kabsch(P, Q):
    """The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


def quaternion_rmsd(P, Q):
    """Rotate matrix P unto Q and calculate the RMSD.

    Based on doi:10.1016/1049-9660(91)90036-O
    """
    rot = quaternion_rotate(P, Q)
    P = np.dot(P,rot)
    return rmsd(P, Q)


def quaternion_transform(r):
    """Get optimal rotation.

    Note: translation will be zero when the centroids of each molecule are the
    same
    """
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3,:3]
    return rot


def makeW(r1,r2,r3,r4=0):
    """Matrix involved in quaternion rotation."""
    W = np.asarray([
             [r4, r3, -r2, r1],
             [-r3, r4, r1, r2],
             [r2, -r1, r4, r3],
             [-r1, -r2, -r3, r4] ])
    return W


def makeQ(r1,r2,r3,r4=0):
    """Matrix involved in quaternion rotation."""
    Q = np.asarray([
             [r4, -r3, r2, r1],
             [r3, r4, -r1, r2],
             [-r2, r1, r4, r3],
             [-r1, -r2, -r3, r4] ])
    return Q


def quaternion_rotate(X, Y):
    """Calculate rotation."""
    N = X.shape[0]
    W = np.asarray([makeW(*Y[k]) for k in range(N)])
    Q = np.asarray([makeQ(*X[k]) for k in range(N)])
    Qt_dot_W = np.asarray([np.dot(Q[k].T,W[k]) for k in range(N)])
    W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
    C1 = -np.sum(Qt_dot_W,axis=0)
    C2 = 0.5*N
    C3 = np.sum(W_minus_Q,axis=0)
    A = np.dot(C3.T,C3)*C2-C1
    eigen = np.linalg.eigh(A)
    r = eigen[1][:,eigen[0].argmax()]
    rot = quaternion_transform(r)
    return rot


def centroid(X):
    """Calculate the centroid from a vectorset X."""
    C = sum(X)/len(X)
    return C


def rmsd(V, W):
    """Calculate Root-mean-square deviation from two sets of vectors V and
    W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)

def main():
    """main"""
    raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
    main()