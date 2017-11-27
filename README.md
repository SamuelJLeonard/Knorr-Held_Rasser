#### What's here?

This repository contains R code I wrote for implementing the disease cluster detection algorithm described "Bayesian Detection of Clusters and Discontinuities in Disease Maps" by Knorr-Held and Rasser (Biometrics 2000).  Their paper is a Bayesian, nonparametric algorithm for detecting spatial clusters of elevated disease risk.  

#### Can I run this algorithm on my own machine?

Yes!  You can run this algorithm against oral cancer rates in the districts of Germany circa 1990.  This is the same data set Drs. Knorr-Held and Rasser analyze in their paper.

To do so download three files included here: 'Knorr Held MCMC.R', 'Knorr Held Functions and Data.R', and 'district_distances.txt' 

Once you have these files on your hard drive, open 'Knorr Held MCMC.R' with R Studio.  A few more instructions are included as comments in this file.

#### How does this algorithm work?  What do the results mean?  What are the mathematical details and where can I learn more?

Further details are in the PDF 'Knorr-Held Model Implementation.pdf'.  You can also read Knorr-Held and Rasser's fantastic paper for a more thorough treatment.

















