#!/usr/bin/env python
#

#import dask_array as da


class interp:
    def __init__(self):

        # Control points for interpolation
        # These are the points that model data
        # is bound to
        self.N = 0 # Polynomial degree
        self.controlNodeType = None
        self.controlPoints = None

        # Target points for interpolation
        self.M = 0 # Polynomial degree
        self.targetNodeType = None
        self.targetPoints = None

        # Quadrature points for integration
        self.qWeights = None
        
        # Derivative matrix for mapping
        # function at control points to its derivative
        # at control points
        self.dMatrix = None

        # Derivative matrix for DG operations
        self.dgMatrix = None

        # Matrix to interpolate to boundaries of computational domain [-1,1]
        self.bMatrix = None

    def load(self, hdf5File):
        """Loads in interpolant data from SELF model output"""
        import h5py

        f = h5py.File(hdf5File, 'r')
        if 'interp' in list(f.keys()):
            self.controlPoints = f['interp/controlpoints']
            self.bMatrix = f['interp/bmatrix']
            self.dMatrix = f['interp/dmatrix'] 
            self.dgMatrix = f['interp/dgmatrix']
            self.iMatrix = f['interp/imatrix']
            self.qweights = f['interp/qweights']
            self.M = self.iMatrix.shape[0]-1
            self.N = self.controlPoints.shape[0]-1 # Polynomial degree is number of quadrature points minus 1

        else:
            print(f"Error: /interp group not found in {hdf5File}.")
            return 1

        return 0


