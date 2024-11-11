#!/usr/bin/env python
#

import self.lagrange as lagrange


# class semline:
#     def __init__(self):

#         self.interp = lagrange.interp()
#         self.nElem = 0 # Number of elements
#         self.x = None # physical coordinates at quadrature points
#         # self.dxds = None # Covariant basis vectors at quadrature points
#         # self.dsdx = None # Contravariant basis vectors at quadrature points
#         # self.J = None # Jacobian at quadrature points
#         self.daskChunkSize=1000 # number of elements per dask chunk

#     def load(self, hdf5File):
#         """Loads in interpolant and geometry data from SELF model output"""
#         import h5py
#         import dask.array as da

#         self.interp.load(hdf5File)

#         f = h5py.File(hdf5File, 'r')


#         if 'controlgrid' in list(f.keys()):

#             d = f['controlgrid/geometry/x/interior'] 
#             self.nElem = d.shape[0]
#             nvar = d.shape[1]
#             N = d.shape[2]
#             self.x = da.from_array(d, chunks=(self.daskChunkSize,nvar,N))

#         else:
#             print(f"Error: /controlgrid group not found in {hdf5File}.")
#             return 1
        
#         if 'targetgrid' in list(f.keys()):
#             d = f['targetgrid/geometry/x/interior'] 
#             self.nElem = d.shape[0]
#             nvar = d.shape[1]
#             N = d.shape[2]
#             self.vis_x = da.from_array(d, chunks=(self.daskChunkSize,nvar,N))
#         else:
#             print(f"Error: /targetgrid group not found in {hdf5File}.")
#             return 1

#         return 0


class semquad:
    def __init__(self):

        self.interp = lagrange.interp()
        self.nElem = 0 # Number of elements
        self.x = None # physical x-coordinates at quadrature points
        self.y = None # physical y-coordinates at quadrature points
        # self.dxds = None # Covariant basis vectors at quadrature points
        # self.dsdx = None # Contravariant basis vectors at quadrature points
        # self.J = None # Jacobian at quadrature points
        self.daskChunkSize=1000 # number of elements per dask chunk

    def load(self, hdf5File):
        """Loads in interpolant and geometry data from SELF model output"""
        import h5py
        import dask.array as da

        self.interp.load(hdf5File)

        f = h5py.File(hdf5File, 'r')
        if 'controlgrid' in list(f.keys()):

            d = f['controlgrid/geometry/x_dim1'] 
            self.nElem = d.shape[0]
            N = d.shape[2]
            self.x = da.from_array(d, chunks=(self.daskChunkSize,N,N))
            d = f['controlgrid/geometry/x_dim2'] 
            self.y = da.from_array(d, chunks=(self.daskChunkSize,N,N))
            self.x_name = "x"
            self.x_units = f['controlgrid/geometry/metadata/units/1']
            self.y_name = "y"
            self.y_units = f['controlgrid/geometry/metadata/units/1']


        else:
            print(f"Error: /controlgrid group not found in {hdf5File}.")
            return 1

        return 0


