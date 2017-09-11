#### Generating a Cluster Configuration 

**Goals/Definitions**

We again start with a country divided into *n* regions.  The paper is interested in determining if clusters of elevated or lowered disease risk exist within this country.  Thus we need a method to divide the *n* regions into distinct clusters.

To form a clustering configuration, select without replacement *k* regions from the set of *n* regions.  Preserve the order in which this selection occurs.  These selected regions are cluster centers, and they completely determine a cluster configuration.  A cluster configuration is formed by assigning each region that is not a cluster center to the cluster center that it is closest to.  

The distance between a non-center and two separate cluster centers could be equal.  Such regions are assigned to the cluster center that appears first in the list of cluster centers.

**Method**

Create a matrix in which the entry *i*, *j* is the distance between non-center *i* and cluster center *j*.  So the rows correspond to non-centers, while the columns correspond to the ordered list of cluster centers.

The first instance of the minimum value of each row is the cluster center that each corresponding non-center will join.  Using this information, create a list containing each cluster center and the corresponding non-centers.
