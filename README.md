# 200101-treeMerging
A simple method to merge two timed phylogenies, using a (fuzzy) MRCA estimate and a fixed/random (posterior distribution) phylogeny/ies.

What do you need?
- Two separate trees of fixed topology OR
- (preferable) a nexus file with a sample of random tree topologies

AND

- an age estimate of where these two trees meet (age of their MRCA) - this can be a fixed number, or a distribution defined e.g. by its central parameter and spread.
