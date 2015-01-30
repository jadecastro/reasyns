%
% ReaSyNS - REActive SYnthesis for Nonlinear Systems
%

reasynsPath

warning off

format compact
clear all
close all

clk = fix(clock);
fid = fopen(['reasyns_log_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6)),'.txt'],'w');

% load automaton
autFile_threeRegionsExamp

Ntrans = length(aut.trans);
Nmodes = length(aut.q);

% load map
regFile_threeRegionsExamp

reg(1) = region(vReg{1});
reg(2) = region(vReg{2});
reg(3) = region(vReg{3});
regBnd = region(vBnd{1});

% Model parameters
polyMdlFun = @(t,x,u) CreateKinematicsPoly(t,x,u);
mdlFun = @(t,x,options,p1,p2) CreateKinematicsNLWayptCtrl2(t,x,options,p1,p2);
ctrlFun = @(x,options,p1,p2) CreateKinematicsNLWayptCtrl2_ctrl(x,options,p1,p2);
sysparams.n = 3; % number of states
sysparams.m = 2; % number of inputs
sysparams.limsNonRegState = [-pi; pi];
sysparams.isCyclic = [0 0 1]';
sysparams.H = eye(2);
sysparams.l = 1.0;  % for car model
sysparams.x_look = 1;
sysparams.e = 0.4;
sysparams.closeEnough = 0.5;
sysparams.distAccept = 0.1;
sys = systemdynamics(polyMdlFun,mdlFun,ctrlFun,sysparams.H,sysparams);

% Feedback controller params
options.ctrloptions_trans.Q = eye(3);
options.ctrloptions_trans.R = 1;
options.ctrloptions_trans.Qf = 20*eye(sysparams.n);

% Sampling parameters
options.Qrand = 2;
options.Ncover = 10000;
options.Nterm = 10;  % number of consecutive failures before coverage terminates
options.coverPct = 0.8;

% Number of funnels/controllers for each state
options.maxFunnelsTrans(1:Nmodes) = 1;
options.maxFunnelsInward(1:Nmodes) = [0 0 0 0];
options.maxFunnelsReactJoin(1:Nmodes) = 10;

% TODO: need?
TstepTraj = 0.02;
options.maxTrajLength = 10/TstepTraj; % 10 seconds; otherwise we're probably spiraling

% Tree depth
% depthTrans = 1;
% depthInward = 1;

% Downsampling
options.sampSkipColl = 5;  % skipped samples in collision check -- higher value speed up collision checks but may miss parts of the trajectory
options.sampSkipFun = 3;  % skipped samples in funnel computation
options.sampSkipValid = 5;  % skipped samples in misbehavior check

options.maxFunTrials = 40;  % number of tries before giving up
options.maxTrials1 = 30; % set to a high value
options.maxTrials2 = 20;
options.maxTrials3 = 2;  % set to a low value because funFail takes care of final point interations.
options.maxNonRegTrials = 30;

trans = vertcat(aut.trans{:});

% for imode = 1:Nmodes
%     qReg = aut.q{imode};
%     qCover{imode} = getCoverPts(vReg,{qReg},1,options.Ncover,sysparams.H,n,sysparams.limsNonRegState);
% end

figure(5)
plot(reg)
% axis([min(regBnd{1}.v(:,1)) max(regBnd{1}.v(:,1)) min(regBnd{1}.v(:,2)) max(regBnd{1}.v(:,2))])


%%  Main

tic


%%
NmodesReach = Nmodes;
% NmodesReach = 1;

tic
toc
fprintf(fid,'%f : Starting ....\n',toc);
j = 0;
for iModeToPatch = 1:NmodesReach
    for itrans = find(trans(:,1)==iModeToPatch)'
        j = j+1;
        ac{j} = [];
    end
end

%%
joinedModes = false*ones(Nmodes,1);  % keep a list of modes which need joining
patchedModes = false*ones(Nmodes,1);  % keep a list of modes which need joining

%%
for iModeToGo = 1:NmodesReach
    
    patchedAndUnvistedModesToPatch = ~(patchedModes(1:iModeToGo));
    
    for iPatchTry = 1:10  % max number of tries
        for iModeToPatch =  find(patchedAndUnvistedModesToPatch)'
            
            %% Reach Operation
            fprintf(fid,'%f : Attempting to patch and join mode %d.\n',toc,iModeToPatch);
            reachIncomplete = true;
            % while reachIncomplete
            [ac_trans,lastTrans] = reachOp_new(sys,reg,regBnd,aut,ac,iModeToPatch,options);
            if ~isempty(lastTrans)
                disp('Reach was unsuccessful. Repeating the Reach.')
                reachIncomplete = true;
                joinedModes(iModeToPatch) = false;
                patchedModes(iModeToPatch) = false;
            else
                disp('Reach was successful.')
                reachIncomplete = false;
                joinedModes(iModeToPatch) = false;
                patchedModes(iModeToPatch) = true;
                tmp = find(trans(:,2)==iModeToPatch)';
                if length(tmp) == 1
                    if ~isempty(ac{tmp})
                        joinedModes(iModeToPatch) = true;
                        break
                    end
                end
            end
            % end
            toc
            
            if ~patchedModes(iModeToPatch)
                iPreModeToPatch = [];
                if ~isempty([ac{lastTrans}])
                    iPreModeToPatch = trans(lastTrans,1);
                end
                disp(['... Failed to perform Reach for mode ',num2str(iModeToPatch),'. Need to re-patch.'])
                toc
                fprintf(fid,'%f : Failed to perform Reach for mode %d. Re-performing the Reach operation\n',toc,iModeToPatch);
                break
            else
                toc
                disp('Finished Reach operation.')
                fprintf(fid,'%f : Finished Reach operation for mode %d.\n',toc,iModeToPatch);
                j = 0;
                for i = find(trans(:,1)==iModeToPatch)'
                    j = j+1;
                    ac{i} = ac_trans(j);
                end
            end
        end
        
        %% Sequence Operation
        ac_inward = sequenceOp_new(sys,reg,regBnd,aut,iModeToPatch,options);
        fprintf(fid,'%f : Finished Sequence operation for mode %d.\n',toc,iModeToPatch);
        
        if ~patchedModes(iModeToGo)
            % update the list of modes which have already been joined but are affected by the patch and concatenate with the list of unvisited modes
            patchedAndUnvistedModesToPatch = find(~(patchedModes(iModeToPatch:iModeToGo))) + (iModeToPatch-1);  % everything left on the todo list
            patchedAndUnvistedModesToPatch = [setdiff(trans(trans(:,1)==iModeToPatch,2),iModeToGo+1:Nmodes); patchedAndUnvistedModesToPatch];  % every successor to the current mode which had been joined
            if ~isempty(iPreModeToPatch)
                patchedAndUnvistedModesToPatch = [setdiff(trans(trans(:,1)==iPreModeToPatch,2),iModeToGo+1:Nmodes); patchedAndUnvistedModesToPatch];
                patchedAndUnvistedModesToPatch = [setdiff(iPreModeToPatch,iModeToGo+1:Nmodes); patchedAndUnvistedModesToPatch];
            end
            
            patchedAndUnvistedModesToPatch = unique(patchedAndUnvistedModesToPatch);
            patchedAndUnvistedModesToPatch
            patchedModes = [~ismember(1:iModeToGo,patchedAndUnvistedModesToPatch)'; patchedModes(iModeToGo+1:Nmodes)];
            
            %%
            if ~isempty(iPreModeToPatch)  % if join failed, also patch the pre mode
                disp(['... Now patching the pre-mode (mode ',num2str(iPreModeToPatch),') ...'])
                toc
                fprintf(fid,'%f : Now Patching the pre-mode %d to the failed mode (%d).\n',toc,iPreModeToPatch,iModeToPatch);
                
                [ac_trans,lastTrans] = reachOp_new(sys,reg,regBnd,aut,ac,iPreModeToPatch,options);
                if ~isempty(lastTrans)
                    disp('Reach was unsuccessful. Repeating the Reach.')
                    reachIncomplete = true;
                    joinedModes(iPreModeToPatch) = false;
                    patchedModes(iPreModeToPatch) = false;
                else
                    disp('Reach was successful.')
                    reachIncomplete = false;
                    joinedModes(iPreModeToPatch) = false;
                    patchedModes(iPreModeToPatch) = true;
                    tmp = find(trans(:,2)==iPreModeToPatch)';
                    if length(tmp) == 1
                        if ~isempty(ac{tmp})
                            joinedModes(iPreModeToPatch) = true;
                            break
                        end
                    end
                end
                toc
                
                if ~patchedModes(iPreModeToPatch)
                    disp(['... Failed to join mode ',num2str(iModeToPatch),'.'])
                    toc
                    fprintf(fid,'%f : Failed to perform Reach for mode %d. Quitting.\n',toc,iModeToPatch);
                    error('Failed to perform Reach');
                else
                    toc
                    fprintf(fid,'%f : Finished Reach operation for mode %d.\n',toc,iPreModeToPatch);
                    j = 0;
                    for i = find(trans(:,1)==iPreModeToPatch)'
                        j = j+1;
                        ac{i} = ac_trans(j);
                    end
                end
                fprintf(fid,'%f : ... Finished Reach operation for pre-mode %d.\n',toc,iPreModeToPatch);
                ac_inward = sequenceOp_new(sys,reg,regBnd,aut,iPreModeToPatch,options);
                fprintf(fid,'%f : ... Finished Sequence operation for pre-mode %d.\n',toc,iPreModeToPatch);
            end
        end
        if all(patchedModes(1:iModeToGo))
            break
        end
    end
end
toc
fprintf(fid,'%f : Finished Reach operation for all modes. Now performing Join operation.\n',toc)

save test


%%
for imode = 1:NmodesReach
    patchedAndUnvistedModesToJoin = ~(joinedModes(1:imode));
    
    for iJoinTry = 1:10  % max number of tries
        for iModeToJoin = find(patchedAndUnvistedModesToJoin)'  % recursively patch successor modes, if necessary
            disp(['... Attempting to join mode ',num2str(iModeToJoin)])
            toc
            fprintf(fid,'%f : Attempting to join mode %d.\n',toc,iModeToJoin);
             % Join Operation
            joinOp
            
            
            if ~joinedModes(iModeToJoin)
                disp(['... Failed to join mode ',num2str(iModeToJoin),'. Need to re-patch.'])
                toc
                fprintf(fid,'%f : Failed to join mode %d. Re-performing the Reach operation\n',toc,iModeToJoin);
                break
            else
                toc
                fprintf(fid,'%f : Successfully joined mode %d.\n',toc,iModeToJoin);
            end
        end
        if ~joinedModes(imode) 
            % update the list of modes which have already been joined but are affected by the patch and concatenate with the list of unvisited modes
            patchedAndUnvistedModesToJoin = find(~(joinedModes(iModeToJoin:imode))) + (iModeToJoin-1);  % everything left on the todo list
            patchedAndUnvistedModesToJoin = [setdiff(trans(trans(:,1)==iModeToJoin,2),imode+1:Nmodes); patchedAndUnvistedModesToJoin];  % every successor to the current mode which had been joined
            if ~joinedModes(iModeToJoin) 
                patchedAndUnvistedModesToJoin = [setdiff(trans(trans(:,1)==iPreModeToJoin,2),imode+1:Nmodes); patchedAndUnvistedModesToJoin];
                patchedAndUnvistedModesToJoin = [setdiff(iPreModeToJoin,imode+1:Nmodes); patchedAndUnvistedModesToJoin];
            end
            
            patchedAndUnvistedModesToJoin = unique(patchedAndUnvistedModesToJoin);
            patchedAndUnvistedModesToJoin
            joinedModes = [~ismember(1:imode,patchedAndUnvistedModesToJoin)'; joinedModes(imode+1:Nmodes)];
            
            % Patch Operation
            disp(['... Patching the failed current mode (mode ',num2str(iModeToJoin),') ...'])
            toc
            fprintf(fid,'%f : Patching the failed current mode %d.\n',toc,iModeToJoin);
            iModeToPatch = iModeToJoin;
            reachOp
            toc
            fprintf(fid,'%f : ... Finished Reach operation for mode %d.\n',toc,iModeToJoin);
            sequenceOp
            fprintf(fid,'%f : ... Finished Sequence operation for mode %d.\n',toc,iModeToJoin);
            
            if ~joinedModes(iModeToJoin)  % if join failed, also patch the pre mode
                disp(['... Now patching the pre-mode (mode ',num2str(iPreModeToJoin),') ...'])
                toc
                fprintf(fid,'%f : Now Patching the pre-mode %d to the failed mode (%d).\n',toc,iPreModeToJoin,iModeToJoin);
                iModeToPatch = iPreModeToJoin;
                reachOp
                toc
                fprintf(fid,'%f : ... Finished Reach operation for pre-mode %d.\n',toc,iPreModeToJoin);
                sequenceOp
                fprintf(fid,'%f : ... Finished Sequence operation for pre-mode %d.\n',toc,iPreModeToJoin);
            end
        end
        if all(joinedModes(1:imode))
            break
        end
    end
end
toc
fprintf(fid,'%f : Finished Join operation for all modes. Now performing Reactive Join operation.\n',toc);

save test_joined

%%
reactJoinedModes = false*ones(Nmodes,1);

for iModeToReactJoin = unique(trans(isReactive,1))'
    disp(['... Attempting to reactively join the mode ',num2str(iModeToReactJoin),'.'])
    % Reactive Join Operation
    for iReactTry = 1:10
        reactivejoinOp
        if reactJoinedModes(iModeToReactJoin), break, end
        toc
        fprintf(fid,'%f : ... Reactive Join operation failed for mode %d.  Retrying... (%d)\n',toc,iModeToReactJoin,iReactTry);
    end
    toc
    fprintf(fid,'%f : Finished Reactive Join operation for mode %d.\n',toc,iModeToReactJoin);
end
toc
fprintf(fid,'%f : Finished Reactive Join operation for all modes. Success!\n',toc);


%%
clk = fix(clock);
if strcmp(modelType,'unicycle')
    eval(['save antsy5_uni_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6))])
elseif strcmp(modelType,'car')
    eval(['save antsy5_car_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6))])
end

toc
fclose(fid);
