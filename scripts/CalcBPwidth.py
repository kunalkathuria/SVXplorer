import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

mu, sigma, rdl = 500, 50, 100
x = mu + sigma*np.random.randn(10000)

if __name__ == "__main__":

        vals = []

        f = open("../results/text/All_Clusters.txt","r")

        for line in f:

                temp = line.split()
		vals.append(400 -int(temp[-1]))		

        vals = np.asarray(vals)

        n, bins, patches = plt.hist(vals,250,normed=1,facecolor='green')

        #print vals[0:100]
        plt.xlabel('Cluster Width')
        plt.ylabel('Frequency')
        plt.axis([0, 600, 0, .3])
        plt.show()
