
filePath = [examplesPath,'/plane'];
configAndProblemDomainName = 'plane';

% load the trajectory
load('plane_trajectories.mat');

% Model parameters
sysparams(1).n = 4; % number of states
sysparams(1).m = 1; % number of inputs
sysparams(1).v = .75;
sysparams(1).max_wdot = 5;  % maximum angular acceleration of the robot
sysparams(1).limsNonRegState = [-pi; pi];
sysparams(1).isCyclic = [0 0 1 0]';

% Feedback controller params
sysparams(1).Q = 1*diag([0.01 0.01 0.01 0.01]); %0.01*eye(sysparams.n);
% sysparams(1).ctrloptions.R = 0.007;
sysparams(1).R = 0.5;
sysparams(1).Qf = 2*eye(sysparams(1).n);  % Final ellipsoid in the transition funnels. Require it to be a ball.

% Region parameters
options.deflationAmount = 0.05; % deflate each region by this amount (in meters) for initial/goal point sampling.

% Define the system 
sys(1) = PlanePlant(sysparams);

% Final rho
options.rhof = 0.1;

% Downsampling
options.sampSkipColl = 5;  % skipped samples in collision check -- higher value speed up collision checks but may miss parts of the trajectory
options.sampSkipFun = 200;  % skipped samples in funnel computation
options.sampSkipValid = 5;  % skipped samples in collision check

% Create a log file
options.doLog = false;
