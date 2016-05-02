reasyns - REActive SYnthesis for Nonlinear Systems
==================================================

Introduction
============

under construction...

Installation
============

First, clone the repo: git@github.com:jadecastro/reasyns.git

This package requires MATLAB 2012b or later.  In addition, it requires installation of the following dependencies.  It may be convenient to install these within the provided 'lib' folder in the project directory, to make the path set-up easier.

1) SeDuMi (version 1.3 tested):
     http://sedumi.ie.lehigh.edu/?page_id=58

2) Ellipsoidal Toolbox (version 1.1.3 tested):
     http://systemanalysisdpt-cmc-msu.github.io/ellipsoids/

3) Multi-Parametric Toolbox (version 2.6.3 tested) 
     http://people.ee.ethz.ch/~mpt/2/downloads/

4) Drake (either Binary or Source installation):
     https://github.com/RobotLocomotion/drake/wiki

5) Mosek (optional, but preferred):
     https://www.mosek.com/resources/downloads
     A license for Mosek can requested free of charge for academic users or purchased for non-academic users.  
     If Mosek is not installed, then SeDuMi is used instead as the optimization engine.

Running an Example
===================

To run the example, first set up the path using the script "reasynsPath.m".  The paths assume that the dependencies reside in the 'lib' folder, but these can be adjusted as needed.

under construction...
