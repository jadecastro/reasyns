#reasyns - REActive SYnthesis for Nonlinear Systems

Introduction
============

Reasyns generates a library of atomic controllers for nonlinear systems satisfying a high-level controller (a finite state machine) synthesized from a reactive mission specification.  It takes in a controller high-level controller, system model, and workspace topology and generates a set of controllers that satisfy the specification with respect to the nonlinear system. 

The toolbox takes in controller/workspace information in the same format as required by the [LTLMoP toolkit](http://ltlmop.github.io/), and generates controllers that may be executed using LTLMoP.  Currently, the execution is implemented on a [separate branch](https://github.com/jadecastro/LTLMoP/tree/reasyns_fast) (which requires no Matlab dependency to run).  The system model is specified in a similar format as required by the [Drake toolbox](http://drake.mit.edu/). 

Installation
============

This package requires MATLAB 2012b or later and installation of the following dependencies.  It may be convenient to install these within the provided `\lib` folder in the project directory, to make the path set-up easier.

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
     
Notes: 
- The source version of Drake is included here as a submodule; note that you may need to run `git submodule init && git submodule update` to fetch the source, then follow [these instructions](http://drake.mit.edu/installation.html) to compile it on your machine.  You may also opt to download a [precompiled release](https://github.com/RobotLocomotion/drake/releases).  If you do, be sure to reflect these changes in `reasynsPath.m`.
- MPT and Ellipsoids come bundled with their own versions of SeDuMi, which do not work with Drake or reasyns. The paths are removed within `reasynsPath.m` (if you've installed these toolboxes in any other location, be sure to change their paths within `reasynsPath.m`).

Running an Example
===================

To run an example, first set up the path, `reasynsPath.m`, then load an example into the workspace. Next, load the parameters and system model (for the included example, this may be done by running `box_pushing.m`), and then run the synthesis engine `synthesizeAtomicControllers(sys, filePath, configAndProblemDomainName, options)`. 

Note that each example in the `/examples` directory requires a region file (`.regions`), a configuration file (`.config`), and an automaton file (`.aut`).  These are of the same format as those generated using LTLMoP (refer to [this tutorial](https://github.com/VerifiableRobotics/LTLMoP/wiki/Tutorial) for information on getting started). 

System models are encoded in the format of "Drake System" objects, of which there are several examples included within Drake.  The `DubinsCar.m` class located in the `/models` folder contains one such example.

The controllers are saved as `.mat` files contained within a folder called `<example path>/reasyns_controllers`.  To use a controller, simply input the name of the desired `.mat` file as an argument in the MotionControllerHandler field.  Refer to the `box_pushing.config` example for specifics.

Finally, execution of the controller requires [this branch](https://github.com/jadecastro/LTLMoP/tree/reasyns_fast) of LTLMoP.  From the `LTLMoP/src` folder, LTLMoP can be called by typing `python specEditor.py <path/to/your/example>`.

