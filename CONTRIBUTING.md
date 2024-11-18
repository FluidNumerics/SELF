# Spectral Element Library in Fortran (SELF)
Copyright 2017-2024 Fluid Numerics LLC

## Contribute Code
1. Find an issue that you are interested in addressing or a feature that you would like to add.
2. Fork the repository associated with the issue to your local GitHub organization. This means that you will have a copy of the repository under your-GitHub-username/SELF.
3. Clone the repository to your local machine using git clone https://github.com/github-username/SELF.git.
4. Create a new branch for your contribution. To help maintainers easily determine the type of code contribution, name your branch using bugfix/issue-NN or feature/issue-NN prefixes.
5. Make the appropriate changes for the issue you are trying to address or the feature that you want to add.
6. When committing your changes, follow the commit guidelines when writing your commit messages.
7. You are encouraged to run the SELF tests using `ctest` on your local system before opening a pull request. If you don't have a GPU, that's ok. Your contribution will be fully tested when you open a pull request.
8. Open a pull request with the upstream SELF repository. In the title, reference the issue number that you worked on. Include a detailed description of the changes you made and why. If you have recommendations for updates to documentation as a result of your changes, please indicate so. If you've added a new routine, you will need to work with the maintainers to develop tests when integrating your new feature in. 

### Code formatting
Each pull request is checked for formatting before running other tests. The `self` project uses [`fprettify`](https://pypi.org/project/fprettify/) for formatting fortran source code. We have included a configuration file in the `self` repository (`fprettify.config`) that can be used for ensuring formatting correctness.

#### Using pre-commit
SELF comes with a `.pre-commmit-config.yaml` file that can be used with [`pre-commit`](https://pre-commit.com/). The benefit of using pre-commit is that this automates applying formatting to all Fortran files in SELF with each commit. To use `pre-commit` :

1. Install `pre-commit` 

```
pip install pre-commit fprettify
```

2. Configure your pre-commit hooks. This command must be run from the root directory of the `self` repository

```
pre-commit install
```

#### Manual formatting
You can run the following to format code to conform to the expected format for `self`.

```
fprettify  './src/' --config-file ./fprettify.config --recursive --case 1 1 1 1
fprettify  './test/' --config-file ./fprettify.config --recursive --case 1 1 1 1
fprettify  './examples/' --config-file ./fprettify.config --recursive --case 1 1 1 1
```
