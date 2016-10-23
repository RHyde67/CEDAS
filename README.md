# CEDAS
Clustering of Evolving Data-streams into Arbitrary Shapes

This branch contains the first version of CEDAS for Matlab pre-2015b.

This version does not use the Matlab 'graph' function included with Matlab 2015b and on. The code runs faster, but does not easily enable visualization of the graph structure.

The data file M-G_3Dx2.csv contains a data stream consisting of 2 Mackey-Glass data streams in 3D, 2 data co-ords and a time parameter.

Run_CEDAS_MG.m reads the Mackey-Glass data and applies the CEDA|S algorithm on every sample. The micor-clusters are displayed as sphere's, coloured by their macro-clusters. The macro-clusters can be seen to marge and separate.

distinguishable_colours is not mine, and is available for downlaod from Matlab's file exchange.
