
filePath = [examplesPath,'/plane'];
configAndProblemDomainName = 'plane';

% Model parameters
sysparams(1).n = 4; % number of states
sysparams(1).m = 1; % number of inputs
sysparams(1).v = .75;
sysparams(1).max_wdot = 5;  % maximum angular acceleration of the robot
sysparams(1).limsNonRegState = [-pi; pi];
sysparams(1).isCyclic = [0 0 1 0]';

% Waypoint steering controller parameters
sysparams(1).e = 0.001;
sysparams(1).closeEnough = 0.2;
sysparams(1).distAccept = 0.1;

% Feedback controller params
sysparams(1).Q = 100*diag([0.01 0.01 0.01 0.01]); %0.01*eye(sysparams.n);
% sysparams(1).ctrloptions.R = 0.007;
sysparams(1).R = 0.5;
sysparams(1).Qf = 20*eye(sysparams(1).n);  % Final ellipsoid in the transition funnels. Require it to be a ball.
sysparams(1).Qf_join = 50*eye(sysparams(1).n);  % Final ellipsoid in the join funnels. Require it to be a ball.

% Sampling parameters
sysparams(1).Qrand = 10.0*eye(sysparams(1).n);
options.Ncover = 10000;
options.Nterm = 10;  % number of consecutive failures before coverage terminates
options.coverPct = 0.8;
options.deflationAmount = 0.05; % deflate each region by this amount (in meters) for initial/goal point sampling.

% Define the system 
sys(1) = PlanePlant(sysparams);

% Number of funnels/controllers for each state
options.maxFunnelsTrans = 1;
options.maxFunnelsInward = 0;
% options.maxFunnelsInward(4) = 1;
options.maxFunnelsReactJoin = 10;

% Trajectory parameters
options.TstepTraj = 0.02;
options.Tfin = 1000;  % Absolute cutoff time
options.maxTrajLength = 200/options.TstepTraj; 

% RRT parameters
options.pathLengthRRT = 0.3;
options.maxNodes = 10;
options.sampleSkipColl = 2;
options.gaussWeight = 0.95;  % Weight [0-1] on Gaussian sampling biasing at qGoal
options.M = diag([0.01 0.01 100 100]);    % Covariance matrixc

% Final rho
options.rhof = 0.02;

% Downsampling
options.sampSkipColl = 5;  % skipped samples in collision check -- higher value speed up collision checks but may miss parts of the trajectory
options.sampSkipFun = 200;  % skipped samples in funnel computation
options.sampSkipValid = 5;  % skipped samples in misbehavior check

options.maxFunTrials = 40;  % number of tries before aborting the current funnel
options.maxTrials1 = 30; % set to a high value
options.maxTrials2 = 1;
options.maxTrials3 = 2;  % set to a low value because funFail takes care of final point interations.
options.maxNonRegTrials = 30;

% Options for computing barrier function-based funnels
options.flagBuildConformingFunnel = false;
options.d = 6;  % order of the barrier function
