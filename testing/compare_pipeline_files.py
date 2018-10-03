#!/usr/bin/env python

import os
import argparse



parser = argparse.ArgumentParser(description='Compare two HDF5 files from SELF-Fluids model output using the hdiff tools.')

parser.add_argument('dir1', metavar='dir1', type=str,
                     help='First directory containing pipeline files for comparison')

parser.add_argument('dir2', metavar='dir2', type=str,
                     help='Second directory containing pipeline files for comparison')

parser.add_argument('pindex', metavar='pindex', type=str,
                     help='Second directory containing pipeline files for comparison')

args = parser.parse_args()
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundarySolution") 
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundary_flux" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundary_gradient_flux_1" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundary_gradient_flux_2" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundary_gradient_flux_3" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/flux_1" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/flux_2" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/flux_3" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/flux_divergence" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/solution" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/solution_gradient_1" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/solution_gradient_2" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/solution_gradient_3" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/source" )
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/tendency") 
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 external/prescribedState") 
os.system("h5diff  -p 1.0E-6 -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 external/externalState") 

#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundarySolution") 
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundary_flux" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundary_gradient_flux_1" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundary_gradient_flux_2" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 boundary/boundary_gradient_flux_3" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/flux_1" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/flux_2" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/flux_3" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/flux_divergence" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/solution" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/solution_gradient_1" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/solution_gradient_2" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/solution_gradient_3" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/source" )
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 internal/tendency") 
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 external/prescribedState") 
#os.system("h5diff -r "+args.dir1+"Pipeline."+args.pindex+".h5 "+args.dir2+"Pipeline."+args.pindex+".h5 external/externalState") 
