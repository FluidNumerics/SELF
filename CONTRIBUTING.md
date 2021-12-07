# Spectral Element Library in Fortran (SELF)
Copyright 2017-2020 Fluid Numerics LLC

## What do I need to know to help?
If you are looking to help to with a code contribution, SELF uses the following programming languages and tools
* Fortran (2008 Standard)
* [HIP](link needed) & [HIPfort](link needed) for GPU acceleration
* [MPI](link needed) for Multi-Process parallelization (multi-core, multi-gpu)
* [CMake](https://cmake.org) (Build System)
* [HDF5](https://www.hdfgroup.org/solutions/hdf5/)

If you are interested in making a code contribution and would like to learn more about the technologies that we use, check out the list below.

* [SELF Specifications](./SPECIFICATIONS.md)
* [Understanding the SELF Software Layout](link needed)
* [Commit Guidelines](link needed)

## Branching Model


## How do I make a contribution?

### Contributing Code as a Maintainer
Maintainers have access to the SELF Jira board that is managed by OctopusSkeleton. When contributing as a maintainer, you must include JIRA ticket numbers in commits, branches, and pull requests.


Developers just have to reference Jira issue keys in commits, branches, pull requests, etc, as described in the table below.

In all cases, the issue key must use to the default Jira key format – that is, two or more uppercase letters ([A-Z][A-Z]+), followed by a hyphen and the issue number. For example, ABC-123.

 	Tool	Instructions
Commits	Bitbucket, GitLab, GitHub, GitHub Enterprise, Fisheye
Other service providers or SCM tools	Include the issue key in the commit message.
For example, a commit message like this "TIS-1 Initial commit" will automatically transition the TIS-1 issue from 'To Do' to 'In Progress'.
Branches	Bitbucket, GitLab, GitHub, GitHub Enterprise, Fisheye
Other service providers or SCM tools	Include the issue key in the branch name when you create the branch.
If you create the branch from the Development panel in a Jira issue, the issue key is added automatically.
For example, if you name your branch "TIS-2_feature", the TIS-2 issue in Jira will automatically transition from 'To Do' to 'In Progress'. (Note that Git doesn't allow spaces in branch names.)
Pull requests	Bitbucket, GitLab, GitHub,
GitHub Enterprise
Other service providers or SCM tools	Do at least one of the following:
Include a commit in the pull request that has the issue key in the commit message. Note, the commit cannot be a merge commit.
Include the issue key in the pull request title.
Ensure that the source branch name includes the issue key.
If you create the pull request from the Branches dialog of the Development panel in a Jira issue, the issue key is added automatically.
For example, if you create a pull request that has "TIS-3" in the title, the TIS-3 issue will automatically transition from 'In Progress' to 'In Review'. If you reopen, decline, or merge the pull request, it will also transition the TIS-3 issue accordingly.
Reviews	Crucible	Include the issue key in the review title when you create the review.
For example, if you name your review "TIS-4 New story" and start the review, the TIS-4 issue will automatically transition from 'In Progress' to 'In Review'. If you reject, abandon, or close the review, it will also transition the TIS-4 issue accordingly.
Builds	Bamboo
Bitbucket Pipelines
Other service providers, SCM tools, or CI/CD Pipelines	For Bamboo, a build is automatically linked to an issue if one of the build's commits includes the issue key in its commit message. The issue key must be included in the commit to activate this feature.
For Pipelines, simply include the issue key in the branch name.
Deployment	Bamboo
Bitbucket Pipelines
Other service providers, SCM tools, or CI/CD Pipelines	A deployment to an environment, such as production or testing, is linked if a commit associated with the deploy contains the issue key in its commit message. The issue key must be included in the commit to activate this feature.

### Contribute Code (Non-maintainer)
1. Find an issue that you are interested in addressing or a feature that you would like to add.
2. Fork the repository associated with the issue to your local GitHub organization. This means that you will have a copy of the repository under your-GitHub-username/SELF.
3. Clone the repository to your local machine using git clone https://github.com/github-username/SELF.git.
4. Create a new branch for your contribution. To help maintainers easily determine the type of code contribution, name your branch using bugfix/issue-NN or feature/issue-NN prefixes.
5. Make the appropriate changes for the issue you are trying to address or the feature that you want to add.
6. When committing your changes, follow the commit guidelines when writing your commit messages.
7. You are encouraged to run the SELF test harness before opening a pull request. These tests will help you make sure you're changes meet formatting guidelines, build successfully, and produce test results within acceptable tolerance levels. If you don't have a GPU, that's ok. Your contribution will be fully tested when you open a pull request.
8. Open a pull request with the upstream SELF repository. In the title, reference the issue number that you worked on. Include a detailed description of the changes you made and why. If you have recommendations for updates to documentation as a result of your changes, please indicate so. If you've added a new routine, you will need to work with the maintainers to develop tests when integrating your new feature in. 

### Code Formatting
To help code maintain consistent formatting, use the `fprettify` toolkit (`pip3 install fprettify`). To format your code before submitting changes,
```
fprettify -i 2 -l 120 --whitespace-type true --whitespace-comma false src/*.F90
```
Verify that no warnings are printed to screen during this process.

All pull requests will be checked against this fprettify call to verify that these format standards are met.

## What does the Code of Conduct mean for me?
Our Code of Conduct means that you are responsible for treating everyone on the project with respect and courtesy regardless of their identity. If you are the victim of any inappropriate behavior or comments as described in our Code of Conduct, we are here for you and will do the best to ensure that the abuser is reprimanded appropriately, per our code.
