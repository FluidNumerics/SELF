# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
#
# Maintainers : support@fluidnumerics.com
# Official Repository : https://github.com/FluidNumerics/self/
#
# Copyright © 2024 Fluid Numerics LLC
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
#
# NOTE #
#
# If you encounter the following error on debian/ubuntu:
#
# ```
# libGL error: MESA-LOADER: failed to open swrast: /usr/lib/dri/swrast_dri.so: cannot open shared object file: No such file or directory (search paths /usr/lib/x86_64-linux-gnu/dri:\$${ORIGIN}/dri:/usr/lib/dri, suffix _dri)
# libGL error: failed to load driver: swrast
# 2024-10-18 01:47:34.339 (   3.885s) [    754DF00AD740]vtkXOpenGLRenderWindow.:651    ERR| vtkXOpenGLRenderWindow (0x3b03f80): Cannot create GLX context.  Aborting.
# ERROR:root:Cannot create GLX context.  Aborting.
# ```
# See https://askubuntu.com/questions/1352158/libgl-error-failed-to-load-drivers-iris-and-swrast-in-ubuntu-20-04 for resolution

import numpy as np
import matplotlib
import pyself.model3d as model3d
import pyvista as pv
import os
import sys
import subprocess
import glob


pmin = -7.5e2
pmax = 7.5e2

def get_sorted_files(pattern):
    files = glob.glob(pattern)
    files.sort(key=os.path.getmtime)  # Sort by modification time
    return files

colormap = matplotlib.colormaps['Reds']
# output video name
video_name = "euler3d_blastwave.mp4"

# Specify the directory you want to search in
directory_path = "/scratch/joe/build/examples/" 

# Get a list of all files in the directory
all_files = os.listdir(directory_path)

# Filter the list to include only .dat files
pickup_files = get_sorted_files(f"{directory_path}/solution.*.h5")

# If you want to get the full path of each file
pickup_files = [os.path.join(directory_path, f) for f in pickup_files]

model = model3d.model()
model.load(pickup_files[0])

pv.start_xvfb()

# pl.show(screenshot="pressure.png")
k = 0
# Example usage of reading each file
for pickup_file in pickup_files:
    print(pickup_file)
    if(k==0):
        model.load(pickup_file)
        pl = pv.Plotter(off_screen=True)
       # slices = model.pvdata.slice(normal=[1, 0, 0])
        pl.add_mesh(model.pvdata,
                    scalars="P",
                    cmap=colormap)
        pl.show(auto_close=False)
        pl.camera_position = 'yz'
        pl.open_movie(video_name)
    else:
        model.update_from_file(pickup_file)
       # slices = model.pvdata.slice(normal=[1, 0, 0])


    pl.write_frame()
    k+=1

pl.close()
