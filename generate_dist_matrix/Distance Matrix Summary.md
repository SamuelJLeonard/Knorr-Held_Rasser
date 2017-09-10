#### Generating a Distance Matrix from an Adjacency Matrix

**Goal/Definitions**

Define the distance between regions *i* and *j* as the minimum number of boundaries one must cross to travel from the interior
of region *i* to the interior of region *j*.

To implement Knorr-Held and Rasser's MCMC, we will need a matrix in which entry *i, j* is the distance between regions *i* and 
*j*.

Contained in this folder is a function which generates this distance matrix from a given adjacency matrix.  An adjacency matrix
is a matrix in which entry *i, j* is 1 if regions *i* and *j* are neighbors, and 0 otherwise.  In this case let the diagonal of 
the adjacency matrix be *0*s.

**Method**

Let *A* be the adjacency matrix of some geographic region.  Consider the entry *i, j* of the matrix *A* * *A* where *i* is not 
equal to *j*.

*i, j$ will be 0 
