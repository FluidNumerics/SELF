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
import self.model2d as model2d
import pyvista as pv
import os
import sys
import subprocess

# Specify the directory you want to search in
directory_path = "/scratch/joe/build/examples/" 

# Get a list of all files in the directory
all_files = os.listdir(directory_path)

# Filter the list to include only .dat files
pickup_files = [f for f in all_files if f.endswith('.h5')]

# If you want to get the full path of each file
pickup_files = [os.path.join(directory_path, f) for f in pickup_files]

model = model2d.model()
# Example usage of reading each file
for pickup_file in pickup_files:
    print(pickup_file)
    model.load(pickup_file)
    colormap = matplotlib.colormaps['coolwarm']

    # Plot the field in 2-d
    pv.start_xvfb()
    pl = pv.Plotter()
    pl.add_mesh(model.pvdata,
                scalars="pressure",
                cmap=colormap,
                clim=[-1.1e-4, 1.1e-4])
    pl.camera_position = 'xy'
    output = pickup_file.split("/")[-1].replace("h5","eps")
    pl.save_graphic(f"./{output}")

    # Convert eps to png using ghostscript
    subprocess.run(f'gs -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -r600 -sDEVICE=pngalpha -sOutputFile={output.replace("eps","png")} {output}',
       shell=True,check=True)


# Create a gif from the png
print("Generating gif from frames")
subprocess.run("ffmpeg -y -framerate 20 -s 1920x1080 -pattern_type glob -i '*.png' -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" linear_euler2d.mp4",shell=True,check=True)
