# Spectral Element Library in Fortran (SELF)
Copyright 2017-2020 Fluid Numerics LLC

## What do I need to know to help?
If you are looking to help to with a code contribution, SELF uses the following programming languages and tools
* Fortran (2008 Standard)
* [HIP](link needed) for GPU acceleration
* [MPI](link needed) for Multi-Process parallelization (multi-core, multi-gpu)
* [CMake](https://cmake.org) (Build System)
* [HDF5](https://www.hdfgroup.org/solutions/hdf5/)

If you are interested in making a code contribution and would like to learn more about the technologies that we use, check out the list below.

* [SELF Specifications](./SPECIFICATIONS.md)
* [Understanding the SELF Software Layout](link needed)
* [Commit Guidelines](link needed)

## Branching Model

## Contribute Code
1. Find an issue that you are interested in addressing or a feature that you would like to add.
2. Fork the repository associated with the issue to your local GitHub organization. This means that you will have a copy of the repository under your-GitHub-username/SELF.
3. Clone the repository to your local machine using git clone https://github.com/github-username/SELF.git.
4. Create a new branch for your contribution. To help maintainers easily determine the type of code contribution, name your branch using bugfix/issue-NN or feature/issue-NN prefixes.
5. Make the appropriate changes for the issue you are trying to address or the feature that you want to add.
6. When committing your changes, follow the commit guidelines when writing your commit messages.
7. You are encouraged to run the SELF test harness before opening a pull request. These tests will help you make sure you're changes meet formatting guidelines, build successfully, and produce test results within acceptable tolerance levels. If you don't have a GPU, that's ok. Your contribution will be fully tested when you open a pull request.
8. Open a pull request with the upstream SELF repository. In the title, reference the issue number that you worked on. Include a detailed description of the changes you made and why. If you have recommendations for updates to documentation as a result of your changes, please indicate so. If you've added a new routine, you will need to work with the maintainers to develop tests when integrating your new feature in. 

### Code formatting
Each pull request is checked for formatting before running other tests. The `feq-parse` project uses [`fprettify`](https://pypi.org/project/fprettify/) for formatting fortran source code. We have included a configuration file in the `feq-parse` repository (`fprettify.config`) that can be used for ensuring formatting correctness. 

You can run the following to format code to conform to the expected format for `feq-parse`.

```
fprettify  './src/' --config-file ./fprettify.config --recursive --case 1 1 1 1
fprettify  './test/' --config-file ./fprettify.config --recursive --case 1 1 1 1
fprettify  './example/' --config-file ./fprettify.config --recursive --case 1 1 1 1
```
