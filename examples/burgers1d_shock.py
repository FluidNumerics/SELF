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

import numpy as np
import matplotlib.pyplot as plt

def read_teplot_curve(file_path):
    """
    Reads a Tecplot curve file and plots the data.

    Parameters:
        file_path (str): The path to the Tecplot curve file.
    """
    # Initialize lists to store data
    x_data = []
    y_data = []

    time = -1.0
    # Open the Tecplot file
    with open(file_path, 'r') as file:
        for line in file:
            # Skip comments and empty lines
            if line.strip().startswith('#') or not line.strip():
                if line.strip().startswith('#TIME'):
                    time = float(line.split()[1])
                continue
            
            # Split the line into components
            data = line.split()
            x_data.append(float(data[0]))
            y_data.append(float(data[1]))

    # Convert lists to numpy arrays for plotting
    x_data = np.array(x_data)
    y_data = np.array(y_data)

    return x_data, y_data, time


import os

# Specify the directory you want to search in
directory_path = '/scratch/joe/build/examples/'# TO DO - set this to the path of your output 

# Get a list of all files in the directory
all_files = os.listdir(directory_path)

# Filter the list to include only .dat files
tecplot_files = [f for f in all_files if f.endswith('.curve')]

# If you want to get the full path of each file
tecplot_files = [os.path.join(directory_path, f) for f in tecplot_files]

# Print the list of Tecplot files
print('Tecplot files:', tecplot_files)

# Example usage of reading each file
for tecplot_file in tecplot_files:
    x,u,t = read_teplot_curve(tecplot_file)

    # Create a new figure with a dark background
    plt.style.use('dark_background')  # Set the background to dark
    plt.figure(figsize=(10, 6))
    plt.plot(x,u, marker='o', linestyle='-', markersize=3, color='white', label='u vs. x')

    # Set x and y limits if provided
    plt.xlim((-0.01, 1.01))
    plt.ylim((-0.01, 1.01))

    # Customize the axes
    plt.title('Burgers Equation - Traveling Shock (s=0.5)', color='white',fontsize=14)
    plt.xlabel('x', color='white',fontsize=12)
    plt.ylabel('u', color='white',fontsize=12)

    plt.text(0.95, 0.95, f"time = {t}", transform=plt.gca().transAxes, ha='right', va='bottom')

    # Customize ticks and grid
    plt.xticks(color='white')
    plt.yticks(color='white')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig(f"{tecplot_file}.png")
    plt.close()
