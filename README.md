# Simple-analysis-codes
this repository contains some simple analysis codes

## count_lifetime.c
this code reads a list of integer numbers from a file, and calculate its probability histogram and the average
Usuage: ./count_lifetime.exe integer-number-file

## get_stat_histogram.c
this code is an updated one of count_lifetime.c. this code can calculate the statistic histogram from a file (its one column), with chosen range and bin-size
How to use it:
  ./get_stat_histogram.exe input-file n_skip column start end bin-size

## degree-distribution.c
this code reads a bounch of GraphGeod files (from ChemNetwork), and calculate the degree distribution from that. Before running the code, it is required to check the name of GraphGeod to be read (the source code should match the actual name), and the total snapshot number (goes from 1 to Nsnaps). This code has been transfered to lifetime-correction-code.

## mygrab-waterid.c
this code reads a file (containing all the water molecules), and grabs several water molecules based on their indexes (the water molecule index goes from 1 to Nwaters). this is used to pick several water clusters for visualizing their local structural pattern. All the coordinates are shifted with respect to the periodic boundary condition.

## matrix-inverse.c
this code is to calculate the inverse matrix using Gauss-Jordan method. it reads the matrix from a file.
  matrix-inverse.exe file-of-original-matrix

## omp-call-chemnetworks.c
this c-code is parallely calling chemnetworks for a bounch of snapshots using OpenMP library (-fopenmp)

## my_ft.f90
this is a simple code of discrete fourier transform from time-domain in unit of ps to frequency-domain in unit of cm^-1, this is naive fourier transform, not fast fourier transform.




