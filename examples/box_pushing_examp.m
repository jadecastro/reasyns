
%
% A box-transportation example 
% Refer to: DeCastro and Kress-Gazit, HSCC 2016 for details.
%

filePath = '/home/jon/Dropbox/Repos/LTLMoP/src/examples/box_pushing';
configAndProblemDomainName = 'box_pushing';

% load automaton
% TODO: input dynamics propositions also
%autFile_complexMap
%autFile_inspectionMap
% autFile_twoRegions
% [aut, transTmp, transTmpNew, transTmpNewWithoutSelfLoops, transWithoutSelfLoops] = ...
%     processAutFileFastSlow('/home/jon/Dropbox/Repos/LTLMoP/src/examples/box_pushing/box_pushing.aut');

load aut_boxPushing
aut = renumberStates(aut);
aut.f = {1, 1, 1, 1, 1, 2, 2};  % enforce certain dynamics depending on the current transition (region activation
%aut.f = {3, 3, 3, 3, 3, 3, 3};  % enforce certain dynamics depending on the current transition (region activation

Ntrans = length(aut.trans);
Nmodes = length(aut.q);

% load map
%regFile_complexMap
regFile_boxPushing
% regFile_twoRegions
% aut.f = {1, 1};

% reg(1) = Region(vReg{1},calibMatrix);
% reg(2) = Region(vReg{2},calibMatrix);
% reg(3) = Region(vReg{3},calibMatrix);
% regBnd = Region(vBnd{1},calibMatrix);

% Model parameters - Unicycle
polyMdlFun{1} = @(t,x,u) CreateKinematicsPoly(t,x,u);
mdlFun{1} = @(t,x,options,p1,p2) CreateKinematicsNLWayptCtrl2(t,x,options,p1,p2);
ctrlFun{1} = @(x,options,p1,p2) CreateKinematicsNLWayptCtrl2_ctrl(x,options,p1,p2);
drakeplant{1} = DubinsPlant();

% Model parameters
sysparams(1).n = 3; % number of states
sysparams(1).m = 2; % number of inputs
sysparams(1).limsNonRegState = [-pi; pi];
sysparams(1).isCyclic = [0 0 1]';
sysparams(1).H = eye(2);
sysparams(1).l = 1.0;  % for car model
sysparams(1).x_look = 1;
sysparams(1).e = 0.1;
sysparams(1).closeEnough = 0.5;
sysparams(1).distAccept = 0.1;

% Feedback controller params
sysparams(1).ctrloptions.Q = 1*diag([0.01 0.01 0.01]); %0.01*eye(sysparams.n);
% sysparams(1).ctrloptions.R = 0.007;
sysparams(1).ctrloptions.R = 0.05;
sysparams(1).ctrloptions.Qf = 1*eye(max([sysparams.n]));

% Model parameters - Holonomic Robot
polyMdlFun{2} = @(t,x,u) HolonomicPoly(t,x,u);
mdlFun{2} = @(t,x,options,p1,p2) HolonomicWayptCtrl(t,x,options,p1,p2);
ctrlFun{2} = @(x,options,p1,p2) HolonomicWayptCtrl_ctrl(x,options,p1,p2);
drakeplant{2} = HolonomicPlant();
sysparams(2) = sysparams(1);


% Model parameters - Dubins Car with fwd velocity as a parameter
polyMdlFun{3} = @(t,x,u) CreateKinematicsVelAsParamPoly(t,x,u);
mdlFun{3} = @(t,x,options,p1,p2) CreateKinematicsVelAsParamNLWayptCtrl2(t,x,options,p1,p2);
ctrlFun{3} = @(x,options,p1,p2) CreateKinematicsNLWayptCtrl2_ctrl(x,options,p1,p2);
drakeplant{3} = DubinsPlantVelocityAsParam();

sysparams(3) = sysparams(1);
sysparams(3).n = 4; % number of states
sysparams(3).isCyclic = [0 0 1 0]';
sysparams(3).limsNonRegState = [-pi 4; pi 4];
sysparams(3).ctrloptions.Q = 1*diag([0.01 0.01 0.01 1e-7]); %0.01*eye(sysparams.n);

% Model parameters
polyMdlFun{4} = @(t,x,u) CreateKinematicsPoly(t,x,u);
mdlFun{4} = @(t,x,options,p1,p2) CreateKinematicsNLWayptCtrl2_bigOmega(t,x,options,p1,p2);
ctrlFun{4} = @(x,options,p1,p2) CreateKinematicsNLWayptCtrl2_bigOmega_ctrl(x,options,p1,p2);
drakeplant{4} = DubinsPlant();
sysparams(4) = sysparams(1);


sys(1) = SystemDynamics(polyMdlFun{1},mdlFun{1},ctrlFun{1},drakeplant{1},sysparams(1).H,sysparams(1));
sys(2) = SystemDynamics(polyMdlFun{2},mdlFun{2},ctrlFun{2},drakeplant{1},sysparams(2).H,sysparams(2));
sys(3) = SystemDynamics(polyMdlFun{3},mdlFun{3},ctrlFun{3},drakeplant{1},sysparams(3).H,sysparams(3));
sys(4) = SystemDynamics(polyMdlFun{4},mdlFun{4},ctrlFun{4},drakeplant{4},sysparams(4).H,sysparams(4));

% Sampling parameters
options.Qrand = 0.1;
options.Ncover = 10000;
options.Nterm = 10;  % number of consecutive failures before coverage terminates
options.coverPct = 0.8;

% Number of funnels/controllers for each state
options.maxFunnelsTrans(1:Nmodes) = 1;
options.maxFunnelsInward(1:Nmodes) = zeros(1,Nmodes);
% options.maxFunnelsInward(4) = 1;
options.maxFunnelsReactJoin(1:Nmodes) = 10;

% TODO: need?
TstepTraj = 0.02;
options.maxTrajLength = 200/TstepTraj; % 10 seconds; in attempt to avoid 

% Tree depth
% depthTrans = 1;
% depthInward = 1;

% Downsampling
options.sampSkipColl = 5;  % skipped samples in collision check -- higher value speed up collision checks but may miss parts of the trajectory
options.sampSkipFun = 50;  % skipped samples in funnel computation
options.sampSkipValid = 5;  % skipped samples in misbehavior check

options.maxFunTrials = 40;  % number of tries before giving up
options.maxTrials1 = 30; % set to a high value
options.maxTrials2 = 20;
options.maxTrials3 = 2;  % set to a low value because funFail takes care of final point interations.
options.maxNonRegTrials = 30;
