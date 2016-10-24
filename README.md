# CEDAS
CEDAS: Clustering of Evolving Data-streams into Arbitrary Shapes

Clustering of Evolving Data-streams into Arbitrary Shapes

This branch contains two versions of CEDAS for Matlab.

## CEDAS.m
This is the full algorithm implementation and requires Matlab 2015b or later as it uses Matlab's built in 'graph' function.

The file CEDAS_IS_Demo.m creates a rando data stream and applies the CEDAS algorithm in an online manner, producing complete cluster results with each data sample. By using the Matlab 'graph' function the full graph structure is available and also displayed.

## CEDAS_Pre2015b.m
This version does not use the Matlab 'graph' function included with Matlab 2015b and on. The code runs faster, but does not easily enable visualization of the graph structure.

The data file M-G_3Dx2.csv contains a data stream consisting of 2 Mackey-Glass data streams in 3D - 2 data co-ords and a time parameter.

Run_CEDAS_MG.m reads the Mackey-Glass data and applies the CEDAS algorithm on every sample. The micro-clusters are displayed as sphere's, coloured by their macro-clusters. The macro-clusters can be seen to marge and separate.

distinguishable_colours is not the work of the CEDAS author. For full details see the headers in the file. It is available for downlaod from Matlab's file exchange.
