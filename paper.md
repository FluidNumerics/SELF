---
title: 'Spectral Element Library in Fortran : A portable, spectrally accurate, multi-GPU accelerated conservation law solver written in Fortran'
tags:
  - Fortran
  - Computational Fluid Dynamics
  - Conservation Laws
  - High performance computing
authors:
  - name: Joseph A. Schoonover
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Garrett Byrd
    equal-contrib: true
    affiliation: 1
  - name: Siddhartha Bishnu
  
affiliations:
 - name: Fluid Numerics, United States
   index: 1
date: 13 August 2025
bibliography: paper.bib

---

# Summary
Hyperbolic and parabolic partial differential equations can be cast as a generic conservation
law, defined by an array of conservative fluxes and non-conservative source terms.
A variety of physical phenomena can be modeled using such conservation laws, which allows for
the development of a software framework that can be used to solve a variety of conservation laws.
This allows for the development of a flexible framework that can be used to define partial differential
equation solvers to model a variety of physical phenomena. The Spectral Element Library in Fortran (SELF)
is an object-oriented implementation of the spectral element method applied to generic conservation laws.
SELF consists of core algorithms for computing weak and strong form differential operations on unstructured
isoparametric grids and higher level abstract models. The abstract models provide a suite of explicit 
Runge-Kutta time integrators and a standard forward stepping pipeline. End users of SELF can create their
own solvers by defining a minimal set of pure fortran functions, which can reduce the time-to-science. All 
methods included in SELF are implemented on CPU and GPU platforms. Domain decomposition is built in to 
SELF, which allows for easy scaling, simply by running SELF programs with multiple MPI ranks. GPU acceleration
is supported by HIP or CUDA, which permits portability across venddor platforms.


# Statement of need

It is becoming apparent that there are now multiple viable platforms for high performance scientific computing,
including GPU accelerated systems with either AMD or Nvidia GPUs. Hardware refreshes at major high performance 
computing facilities, which provide compute cycles to academic researchers, require costly software porting
activities which distract researchers from research activities.

* Why a generic software ?
* Why multiple GPU vendor support ?
* 


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge Jorge Galvez-Vallejo and Jonathan Moore for their support of the 
Spectral Element Library in Fortran.

# References