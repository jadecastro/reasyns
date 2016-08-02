
filePath = [examplesPath,'/dubins_car'];
configAndProblemDomainName = 'dubins';

% Model parameters
sysparams(1).n = 3; % number of states
sysparams(1).m = 1; % number of inputs
sysparams(1).v = 0.08;
sysparams(1).max_w = 3000;  % maximum angular velocity of the robot
sysparams(1).limsNonRegState = [-pi; pi];
sysparams(1).isCyclic = [0 0 1]';

% Waypoint steering controller parameters
sysparams(1).e = 0.001;
sysparams(1).closeEnough = 0.02;
sysparams(1).distAccept = 0.01;

% Feedback controller params
sysparams(1).Q = 100*diag([0.01 0.01 0.01]); %0.01*eye(sysparams.n);
% sysparams(1).ctrloptions.R = 0.007;
sysparams(1).R = 0.05*1e0;
sysparams(1).Qf = 10000*eye(sysparams(1).n);  % To enable 'fast' rejection of bad trajectories (prior to the funnel-creation step), we require the final ellipsoid to be a ball

% Sampling parameters
sysparams(1).Qrand = 0.05*eye(sysparams(1).n);
options.Ncover = 10000;
options.Nterm = 10;  % number of consecutive failures before coverage terminates
options.coverPct = 0.8;
options.deflationAmount = 0.0005; % deflate each region by this amount (in meters) for initial/goal point sampling.

% Define the system 
sys(1) = UnicyclePlant(sysparams);

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
options.pathLengthRRT = 0.2;
options.maxNodes = 200;
options.sampleSkipColl = 2;
options.gaussWeight = 0.9;  % Weight [0-1] on Gaussian sampling biasing at qGoal
options.M = diag([0.1, 0.1, 1]);    % Covariance matrix

% Final rho
options.rhof = .5;

% Downsampling
options.sampSkipColl = 5;  % skipped samples in collision check -- higher value speed up collision checks but may miss parts of the trajectory
options.sampSkipFun = 200;  % skipped samples in funnel computation
options.sampSkipValid = 5;  % skipped samples in misbehavior check

options.maxFunTrials = 40;  % number of tries before aborting the current funnel
options.maxTrials1 = 30; % set to a high value
options.maxTrials2 = 5;
options.maxTrials3 = 2;  % set to a low value because funFail takes care of final point interations.
options.maxNonRegTrials = 30;

% Options for computing barrier function-based funnels
options.flagBuildConformingFunnel = false;
options.d = 6;  % order of the barrier function
