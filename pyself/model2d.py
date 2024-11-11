#!/usr/bin/env python
#


# Other SELF modules
import self.geometry as geometry

class model:
    def __init__(self):
        self.solution = None
        self.pvdata = None # Pyvista data
        self.varnames = None
        self.varunits = None
        self.geom = geometry.semquad()

    def load(self, hdf5File):
        """Loads in 2-D model from SELF model output"""
        import h5py
        import dask.array as da

        self.geom.load(hdf5File)

        f = h5py.File(hdf5File, 'r')
        self.varnames = []
            
        if 'controlgrid' in list(f.keys()):

            controlgrid = f['controlgrid']
            for group_name in controlgrid.keys():

                if( group_name == 'geometry' ):
                    continue

                group = controlgrid[group_name]
                # Create a list to hold data for this group
                setattr(self, group_name, [])
                group_data = getattr(self, group_name)
                print(f"Loading {group_name} group")
                
                # Load metadata information
                if( 'metadata' in list(group.keys()) ):
                    for v in group[f"metadata/name"].keys():

                        name = group[f"metadata/name/{v}"].asstr()[()][0]
                        try:
                            units = group[f"metadata/units/{v}"].asstr()[()][0]
                        except:
                            units = "error"

                        group_data.append({
                            "name": name,
                            "units": units,
                            'data': None
                        })
                else:
                    print(f"Error: /controlgrid/{group_name}/metadata group not found in {hdf5File}.")
                    return 1

                for k in group.keys():
                    k_decoded = k.encode('utf-8').decode('utf-8')
                    if k == 'metadata':
                        continue
                    else:
                        print(f"Loading {k_decoded} field")
                        # Load the actual data
                        d = group[k]
                        N = d.shape[2]

                        # Find index for this field
                        i = 0
                        for sol in group_data:
                            if sol['name'] == k_decoded:
                                break
                            else:
                                i += 1

                        group_data[i]['data'] = da.from_array(d, chunks=(self.geom.daskChunkSize, N, N))

            self.generate_pyvista()

        else:
            print(f"Error: /controlgrid group not found in {hdf5File}.")
            return 1

        return 0 

    def generate_pyvista(self):
        """Generates pyvista polyData for each solution variable for plotting"""
        import numpy as np
        import pyvista as pv

        (nelem, nx, ny) = self.solution[0]['data'].shape
        n_points = nelem*nx*ny
        n_faces = nelem*(nx-1)*(ny-1)

        # Need to use the plot mesh to create a flat list of (x,y,z=0) points
        # number of points = (M+1)*(M+1)*nelem
        # dimension ordering (i,j,iel)
        # Get the x-y points in flattened array for building unstructured data
        np_points = np.zeros((n_points,3))
        np_points[:,0] = self.geom.x.flatten()
        np_points[:,1] = self.geom.y.flatten()

        # Need to construct the faces from here..
        # Number of faces = M*M*nelem
        faces = np.zeros((n_faces,5),dtype=np.int64)
        fid = 0
        for iel in range(0,nelem):
            for j in range(0,ny-1):
                for i in range(0,nx-1):
                    # lower left corner
                    n0 = i + nx*( j + ny*iel )
                    # lower right corner
                    n1 = i+1 + nx*( j + ny*iel )

                    # upper right corner
                    n2 = i+1 + nx*( j+1 + ny*iel )

                    # upper left corner
                    n3 = i + nx*( j+1 + ny*iel )

                    faces[fid,:] = [4, n0, n1, n2, n3]
                    fid += 1
                    
        self.pvdata = pv.PolyData(np_points, faces)

        # Load fields into pvdata
        k = 0
        for attr in self.__dict__:
            if not attr in ['pvdata','varnames','varunits','geom']:
                controlgroup = getattr(self, attr)
                #print(f"Loading {attr} into pvdata")
                for var in controlgroup:
                 #   print(f"Loading {var['name']} into pvdata")
                    self.pvdata.point_data.set_array(var['data'].flatten(),var['name'])
                    k+=1
            
        print(self.pvdata)
