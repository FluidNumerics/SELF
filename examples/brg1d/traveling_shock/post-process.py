#!/usr/bin/env python

import self.model1d as model1d
import inspect, os.path
import matplotlib.pyplot as plt

plt.style.use('dark_background')
plt.switch_backend('agg')

# Get full path to examples/
# From https://stackoverflow.com/questions/2632199/how-do-i-get-the-path-of-the-current-executed-file-in-python
filename = inspect.getframeinfo(inspect.currentframe()).filename
path     = os.path.dirname(os.path.abspath(filename))

model = model1d.model()

for io in os.listdir(path):
    if io.endswith(".h5"):

        # Load model data
        print(io)
        model.load(f'{path}/{io}')
        iterate = io.split(".")[1]

        # Create a new figure
        plt.figure()

        # plot this time level
        plt.plot( model.geom.vis_x.flatten(), model.vis_solution.flatten(), '-o', color='white', linewidth=1, markersize=3, label="u vs. x" )
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('x')
        plt.ylabel('u')
        plt.legend()
        plt.grid(color='gray', linestyle='--', linewidth=1)
        
        # to do - get the time level from the hdf5 output

        # Save the figure to a PNG file
        plt.savefig(f'solution_{iterate}.png')

        # Close the figure to free up memory
        plt.close()

#model.load(f'{path}/data/1d/solution.h5')
#attr=dir(model)





## Interpolate from quadrature mesh to plot mesh
#plt.plot( model.geom.vis_x.flatten(), model.vis_solution.flatten() )
#plt.show()
