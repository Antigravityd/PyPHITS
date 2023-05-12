---
title: "PyPHITS: A Python porcelain for JAEA's PHITS"
tags:
  - Python
  - physics
  - nuclear physics
  - monte carlo
  - radiation
authors:
  - name: Duncan Wilkie
    affiliation: "1, 2, 3"
affiliations:
  - name: Undergraduate, Louisiana State University Department of Physics and Astronomy, USA
    index: 1
  - name: Undergraduate, Louisiana State University Department of Mathematics, USA
    index: 2
  - name: Research Scientist, Atlantis Industries, USA
    index: 3
date: 12 May 2023
bibliography: paper.bib
---

# Summary

The Japanese Atomic Energy Agency's Particle Heavy-Ion Transport code System (PHITS) [@sato2018]
is a Monte Carlo radiation transport program, simulating the effects of physical processes from electron scattering to nuclear fission.
It's used extensively, for example, in particle accelerator design and radiation therapy analysis.
Applying PHITS to study the biological effects of space radiation, it has become apparent that the program could become
substantially more useful if better-integrated into modern computing environments.
Such integration is this program's _raison d'être._

# Statement of need

PHITS has its own file format---to use it, one must learn its intricacies, many of which appear to be anachronisms from the days
when Fortran was spelled with capital letters (its manual is 342 pages, as a rough severity estimate).
This is a testament to the utility of the software, to the reliability of its implementation, _a la_ Lindy's Law.
Nevertheless, this remains an impediment to its continued adoption.
It impacts existing experts too, as mainframe-style programming impedes many interactive and meta-level workflows
single-user machines have enabled in recent decades.
AI integration is the specific example motivating this work: using PHITS simulations alongside machine learning frameworks
like TensorFlow or Keras proves difficult, as each training step must output to, interpret, and modify several PHITS-specific data formats.

PyPHITS rectifies these issues, by providing an API for interacting with PHITS from Python.
The language has become a _de facto_ standard choice for many physical, data-analytic, and machine-learning workflows;
the plethora of extant libraries in the language mean that many of these workflows are now a function-call away,
rather than requiring a _ad hoc,_ purpose-specific code be written first.
It also presents a user-interface that's reasonably well-conformant to design standards in the Python ecosystem.
Any new user of PHITS familiar with Python will be able to focus on the semantics of the program, rather than the syntax of input files.
Moreover, his strengthened knowledge of Python through use of PHITS will be transferrable to any other program.

Not providing an abstracted interface seems to be an exception to the rule for particle transport codes.
Geant4 [@allison2016; @allison2006; @allison2003] provides geant4py,
and FLUKA [@ahdida2022; @battistoni2015] has the graphical "flair" environment [@vlachoudis2009].
Some people nevertheless use PHITS in preference to these other programs, so correcting what these other projects percieve as a problem
will bring those advantages to those for whom the interface is too great a cost to bear.

There is recent work in space medical physics involving PHITS-assisted machine learning [@taylor2023].
The implementations thus far are single-use and poorly-tested, and artificially restrict the search space in interest of
software engineering expedience.
This is despite substantial commercial and defense interest in the service such models provide.
PyPHITS provides a framework by which such models can be hope to be repeatably and reliably trained,
for wider classes of geometries, with better results, and with less waste of developer and scientist time.

The project is open-source, so that technologies which further advance the objectives oulined above may be quickly
and easily developed.
One such possibility is integration with existing 3D modelling software (Blender is open-source also, and extensible in Python)
for geometry production.
This would enable a pipeline from engineering CAD model to particle transport simulation to physical conclusions
that's nearly hands-off, making the nuclear physics power of PHITS a part of an iterative hardware design process.

# References

# Acknowledgements
