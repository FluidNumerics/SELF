#!/usr/bin/env python

import os
import argparse



parser = argparse.ArgumentParser(description='Compare two HDF5 files from SELF-Fluids model output using the hdiff tools.')

parser.add_argument('file1', metavar='file1', type=str,
                     help='First file for comparison')

parser.add_argument('file2', metavar='file2', type=str,
                     help='Second file for comparison')

args = parser.parse_args()



os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /model_output/x_momentum")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /model_output/y_momentum")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /model_output/z_momentum")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /model_output/density")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /model_output/density_weighted_temperature")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /model_output/density_weighted_tracer")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /model_output/pressure")


os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /static/x_momentum")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /static/y_momentum")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /static/z_momentum")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /static/density")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /static/density_weighted_temperature")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /static/density_weighted_tracer")
os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /static/pressure")

os.system("h5diff --relative=1.0E-6 -r "+args.file1+" "+args.file2+" /model_conditions/drag")
