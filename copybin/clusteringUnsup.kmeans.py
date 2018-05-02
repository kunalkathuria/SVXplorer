import sys
import numpy as np
from sklearn.cluster import KMeans

if __name__=="__main__":

    f=open(sys.argv[1],"r")
    A = np.loadtxt(f,dtype=int)
    print A
    kmeans = KMeans(n_clusters=int(sys.argv[2]))
    kmeans = kmeans.fit(A)
    label = kmeans.predict(A)
    centroids = kmeans.cluster_centers_
    print label

