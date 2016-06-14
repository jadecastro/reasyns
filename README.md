#reasyns - REActive SYnthesis for Nonlinear Systems

Introduction
============

Reasyns generates a library of atomic controllers for nonlinear systems satisfying a high-level controller (a finite state machine) synthesized from a reactive mission specification.  It takes in a controller high-level controller, system model, and workspace topology and generates a set of controllers that satisfy the specification with respect to the nonlinear system. 

The toolbox takes in controller/workspace information in the same format as required by the [LTLMoP toolkit](http://ltlmop.github.io/), and generates controllers that may be executed using LTLMoP.  Currently, the execution is implemented on a [separate branch](https://github.com/jadecastro/LTLMoP/tree/reasyns_fast) (which requires no Matlab dependency to run).  The system model is specified in a similar format as required by the [Drake toolbox](http://drake.mit.edu/). 

Installation
============

This package requires MATLAB 2012b or later (setup requires python 2.7 or later). 

Installation is done by running the setup script using the command `python setup.py` inside the project directory.

Upon running `setup.py`, the following dependencies will be automatically installed within the `\lib` folder in the project directory and an initialization script, `reasyns_init.m`, will be created.  If any of the dependencies exist, then the installation can be skipped, and the script will then simply create `reasyns_init.m`.

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
     
Note: 
- MPT and Ellipsoids come bundled with their own versions of SeDuMi, which do not work with Drake or reasyns. The paths are removed within `reasynsPath.m` (if you've installed these toolboxes in any other location, be sure to change their paths within `reasynsPath.m`).

Running an Example
==================

We will refer to the included `box_pushing` example (in the `/examples` directory) as a use case for illustrating how to synthesize and execute controllers. 

Synthesis
---------
1) Set up the path by running `reasyns_init.m`.
2) Run `box_pushing.m` to load the synthesis settings and system model into the workspace.
3) Run the synthesis engine,
```
synthesizeAtomicControllers(sys, filePath, configAndProblemDomainName, options)
```
4) After a set of funnels has been found satisfying the task, the synthesizer will prompt whether or not to continue generating funnels. Choosing this option will increase funnel coverage over the state space, but could take considerable additional time.

Note that the controllers are saved within `.mat` files inside the `reasyns_controllers` directory within the `/box_pushing` directory.

Execution
---------
Execution of the controller requires LTLMoP; in particular, [this branch](https://github.com/jadecastro/LTLMoP/tree/reasyns_fast) of LTLMoP.  
1) Set up the configuration file.  For the `box_pushing` example, this only requires that the name of the desired controller be inserted as an argument in the MotionControllerHandler field within the `configs/box_pushing.config` file.  
2) From within the `LTLMoP/src` folder, LTLMoP can be called by typing `python specEditor.py <path/to/your/example>`.

Creating an Example
===================

High-Level Controller
---------------------
Reasyns requires a high-level controller from which a palette of atomic controllers will be synthesized.  Examples (`/examples` directory) require a region file (`.regions`), a configuration file (`.config`), and an explicit-state automaton file (`.aut`) in JTLV format.  These files may be generated using LTLMoP (refer to [this tutorial](https://github.com/VerifiableRobotics/LTLMoP/wiki/Tutorial) for information on getting started). 

System Model
------------
System models are encoded in the format of `DrakeSystem` objects, of which there are several examples included within Drake.  The `DubinsCar.m` class located in the `/models` folder contains one such example.
1) Construct a class derived from `DrakeSystem` or one of its subclasses.  The class should include the following parameters:
- sysparams.n : number of states
- sysparams.m : number of inputs

And the following methods:
- dynamics
- getInitialState
- state2SEconfig

Refer to the examples (e.g. `UnicyclePlant.m`) for details.
2) Create symbolic gradients by running the function `generateGradients` using the new model. Type `help generateGradients` for details.
