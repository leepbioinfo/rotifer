

import math
import numpy as np

def NE_dn(w, t):
    w = w * w.shape[0]
    D = sum(abs(w))

    if t == 'ave':
        D = np.divide(1, D)
        wn = D*w

    elif type == 'gph':
        D = np.divide(1, np.sqrt(D))
        wn = D * ( w * D )

    return wn

def dominateset(aff_matrix, NR_OF_KNN):


    pass

def TransitionFields(W):
    zeroindex = sum(W) == 0

    W = W * W.shape[0]
    W = NE_dn(W, 'ave')
    w = np.sqrt(sum(abs(W)))

    W = np.multiply(W, np.repeat(w,w.shape[0] 1)) # modify here

    W = W * np.transpose(W)

    Wnew = W
    Wnew[zeroindex,:] = 0
    Wnew[,zeroindex] = 0

    return Wnew

def Network_Enhancement(W_in, order = 2, K = '', alpha=0.9):
    '''
    W_in the input network size N x N
    K number of neighbors
    alpha the regularization parameter
    order determines the extent of diffusion, typical values: 0.5, 1, 2
    '''

    K = min(20, math.ceil(W_in.shape[0]/10))

    # Element wise multiplication, identity matrix
    W_in1 = np.multiply(W_in,(1-np.eye(W_in.shape[0])))

    # Column higher 0
    zeroindex = sum(abs(m))> 0

    W0 = W_in[zeroindex]

    # Function NE_dn
    W = NE_dn(W0, 'ave')
    W = (W+W.transpose())/2

    DD = sum(abs(W0))

    # Base case
    if len(np.unique(W)) == 2:
        P = W

    else:
        P = np.multiply(dominateset(abs(W), min(K, W.shape[0]-1)),
                        np.sign(W))

    P = P + np.eye(P.shape[0] + np.diag(np.diag(sum(abs(P.transpose()))))
                   )

    P = TransitionFields(P)

    D, U = np.linalg.eig(P) # D vectors

    d = np.diag(D-eps).real # define eps here

    d = np.dot(1-alpha, d) / np.multiply(1-alpha, np.power(d, order))

    D = np.diag(d.real)

    W = np.multiply(np.multiply(U, D), np.transpose(U))

    W = np.multiply(W, (1-np.eye(W.shape[0]))) / np.repeat(1-np.diag(W), W.shape[0])

    W = D * W

    W[W<0] = 0

    W = (W + np.transpose(W))/2

    W_out = zeros(W_in.shape)

    W_out(W_in.shape[0], W_in.shape[0]) =  W

    return W_out
