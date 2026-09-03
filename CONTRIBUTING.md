# Spectral Element Library in Fortran (SELF)
Copyright 2017-2025 Fluid Numerics LLC

## Contribute Code
1. Find an issue that you are interested in addressing or a feature that you would like to add.
2. Fork the repository associated with the issue to your local GitHub organization. This means that you will have a copy of the repository under your-GitHub-username/SELF.
3. Clone the repository to your local machine using git clone https://github.com/github-username/SELF.git.
4. Create a new branch for your contribution. To help maintainers easily determine the type of code contribution, name your branch using bugfix/issue-NN or feature/issue-NN prefixes.
5. Make the appropriate changes for the issue you are trying to address or the feature that you want to add.
6. When committing your changes, follow the commit guidelines when writing your commit messages.
7. You are encouraged to run the SELF tests using `ctest` on your local system before opening a pull request. If you don't have a GPU, that's ok. Your contribution will be fully tested when you open a pull request.
8. Open a pull request with the upstream SELF repository. In the title, reference the issue number that you worked on. Include a detailed description of the changes you made and why. If you have recommendations for updates to documentation as a result of your changes, please indicate so. If you've added a new routine, you will need to work with the maintainers to develop tests when integrating your new feature in. 

### Opening a pull request
Pull requests use the template in `.github/PULL_REQUEST_TEMPLATE.md`, and the `pr-template-check` job verifies that you have filled it in. Four sections are required.

1. `## Scope` states what the pull request changes. A reviewer reads this to know what they are being asked to look at, and a change that is not described here does not belong in the pull request.
2. `## Out of scope` lists problems you found while working and are deliberately not fixing. Reference an issue for each, or say that it is not yet filed. Write "None" if you found nothing. Anything listed here is off limits for review comments, so be specific.
3. `## Issues resolved` links the issue this pull request closes, with a closing keyword such as `Fixes #123` so that GitHub closes the issue on merge. Every pull request resolves a filed issue; if none exists, open one first.
4. `## Tests introduced` names the new tests and what each establishes. For a bug fix, name the test that fails against the old code and say why the existing suite did not catch the bug. A bug that reached `main` is evidence of a gap in the tests, and closing that gap is part of the fix.

A maintainer can bypass the check on an urgent fix by applying the `skip-pr-checks` label.

### Coding and documentation style
Formatting is one thing and style is another. `fprettify` handles indentation, alignment, and keyword case; it cannot see how a comment or a documentation page is written. The `style-check` job covers that, scoring each file a pull request changes.

The rules, and the reasoning behind each one, are in [the style guide](docs/Contributing/StyleGuide.md). To run the checks yourself:

```shell
python3 .github/scripts/style_check.py       # src/, test/, examples/
python3 .github/scripts/docs_style_check.py  # docs/ and the root pages
```

Both accept `--files` for an explicit list and `--stats` to report the measured metrics without failing. The `pre-commit` hooks described below run them on the files each commit touches.

Contributors working with a coding agent should also read `.github/copilot-instructions.md` and the per-directory files under `.github/instructions/`, which carry the same rules in the form the agents read. Claude Code users can select the matching output style with `/output-style self`.

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

## Testing code
The `FluidNumerics/SELF` main repository has CI implemented to execute all tests under the `test/` and `examples/` subdirectories. If you are adding new features, make sure that you include new tests and examples that exercise your contribution completely. By "completely" we mean that all new lines of code are exercises **and** you have sufficiently sophisticated tests to cover the range of possible inputs and scenarios you intend to support.

You can test locally on your system after building the code. Suppose you've built the code with the following commands

```shell
mkdir ~/SELF/build
cd ~/SELF/build
cmake {CMAKE_OPTIONS} ~/SELF/
make -j
```

You can then use `ctest` to run all of the tests from your build directory

```shell
ctest ~/SELF/build
```

At the least, we ask that you verify that all tests pass on *your* system, since that likely matters most to you. In our CI environment, we will test on a variety of platforms to ensure SELF passes tests on systems we want to support. 

### Cluster Access
If you contribute code that resolves issues labeled with "free compute reward", you can reach out to support@fluidnumerics.com to request an account on [Fluid Numerics' Galapagos Cluster](https://galapagos.fluidnumerics.com). This will provide you with access to a variety of GPU accelerator platforms for testing in addition to a variety of curated software environments managed by Fluid Numerics. Fluid Numerics provides these accounts to help lower the barrier to entry in being a contributor on the SELF repository and to promote collaboration with the broader community. 

