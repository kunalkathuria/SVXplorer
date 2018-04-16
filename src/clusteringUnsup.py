import sys
import numpy as np
from sklearn.neighbors import NearestNeighbors

if __name__=="__main__":

    f=open(sys.argv[1],"r")
    A = np.loadtxt(f,dtype=float)
    print A
    neigh = NearestNeighbors(n_neighbors=1)
    neigh.fit(A)
    print neigh.kneighbors(A, return_distance=False)
    neigh.predict(A)

