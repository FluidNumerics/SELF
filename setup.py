from setuptools import setup

setup(
    name='pyself',
    version='0.0.1',
    description='A python interface for the Spectral Element Library in Fortran',
    url='https://github.com/fluidnumerics/self',
    author='Fluid Numerics',
    author_email='support@fluidnumerics.com',
    license='3-Clause BSD with Attribution',
    packages=['pyself'],
    install_requires=['h5py>=3.7.0',
                      'dask',
                      'pyvista',
                      'imageio[ffmpeg]'],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3'
    ],
)
