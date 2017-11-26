This repository contains R code I wrote for implementing the disease cluster detection algorithm described "Bayesian Detection of Clusters and Discontinuities in Disease Maps" by Knorr-Held and Rasser (Biometrics 2000).  Their paper is a Bayesian, nonparametric algorithm for detecting spatial clusters of elevated disease risk.  

You can run this algorithm against oral cancer rates in the districts of Germany circa 1990.  This is the same data set Drs. Knorr-Held and Rasser used in their paper.

You'll need to download three files to run the algorithm yourself: 'Knorr Held MCMC.R', 'Knorr Held Functions and Data.R', and 'district_distances.txt' 

Once you have these files on your hard drive, open 'Knorr Held MCMC.R' with R Studio.  A few more instructions are included as comments in this file.

If you'd rather just view my results, 'Knorr-Held Rasser Model Implementation.pdf' is a summary of my implementation of the Knorr-Held and Rasser algorithm.  There you'll find plots of the results of my implementation, along with a brief summary of how the model works.














