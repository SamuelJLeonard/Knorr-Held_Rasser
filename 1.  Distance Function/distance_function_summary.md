### Generating a Distance Matrix from an Adjacency Matrix

**Goal/Definitions**

Consider a country divided into *n* regions, each region numbered *1,...,n*.  Let *i* and *j* be distinct 
regions in this country.  Knorr-Held and Rasser define the distance between regions *i* and *j* as the minimum number of 
boundaries one must cross to travel from the interior of region *i* to the interior of region *j*.

To implement their MCMC we will need the distances between each pair of regions, organized in a matrix for convenience.

Contained in this folder is a function which generates this distance matrix from a given adjacency matrix.  An adjacency 
matrix is a matrix in which entry *i, j* is *1* if regions *i* and *j* are neighbors, and *0* otherwise.  In this case let 
the diagonal of the adjacency matrix be *0*s.  (The diagonal could just as easily be *1*, we are not concerned with distances
of regions to themselves at all.  Using *0* makes the math slightly faster later on).

**Method**

Our goal is a matrix containing the distances between each region.  Let *A* be our adjacency matrix.  And consider the list of 
*n* matrices *A*,  (*A* * *A*),  (*A* * *A* * *A*),  .... ,  (*A^n*).

The key insight here is that the distance between arbitrary regions *i* and *j* is the location in the above list of first matrix 
for which entry *i*, *j* is nonzero.

For example, if regions *i* and *j* are adjacent then their distance is obviously *1*, the first matrix in the list.

If regions *i* and *j* have distance *2* they are not adjacent but are both adjacent at to least one of the
same regions.  Call one such region *k*.

In the matrix *(A * A)*, entry *i*, *j* is the dot product of the *i*th row and *j*th column of *A*.  Because both are adjacent
to *k*, we know that the *k*th entry of the row vector corresponding to *i* is *1*, and the *k*th entry of the column vector 
corresponding to *j* is also *1*.  Hence their dot product will be greater than or equal to *1*.

In fact, the entry *i*, *j* of the matrix *(A * A)* will be the number of regions that both *i* and *j* are adjacent to.  This
number might be useful in another context, but we don't need it.  All that's important is that the first nonzero entry in the 
list of matrices described above occurs in *(A * A)*, the second element of the list.

It's not hard to see that this logic holds for further powers of *A*.  

The R function in this folder duplicates this logic.  Create a list of all possible powers of *A*, then find for each pair *i*
and *j* the first nonzero entry.
