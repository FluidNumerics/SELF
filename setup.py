from setuptools import setup

setup(
    name='pyself',
    version='0.1.0',
    description='A python interface for the Spectral Element Library in Fortran',
    url='https://github.com/fluidnumerics/self',
    author='Dr. Joe Schoonover',
    author_email='joe@fluidnumerics.com',
    license='Researcher Software License',
    packages=['self'],
    install_requires=['h5py>=3.7.0',
                      'dask',
                      'pyvista'],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3'
    ],
)
