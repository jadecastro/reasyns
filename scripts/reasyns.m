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

% make copies of files that we may be modifying later on
resetAllInputFiles

% load automaton
% TODO: input dynamics propositions also
%autFile_complexMap
%autFile_inspectionMap
% autFile_twoRegions
% [aut, transTmp, transTmpNew, transTmpNewWithoutSelfLoops, transWithoutSelfLoops] = ...
%     processAutFileFastSlow('/home/jon/Dropbox/Repos/LTLMoP/src/examples/box_pushing/box_pushing.aut');

load aut_boxPushing
aut = renumberStates(aut);
% aut.f = {2, 1, 1, 1, 1, 2, 2};  % enforce certain dynamics depending on the current transition (region activation
aut.f = {3, 3, 3, 3, 3, 3, 3};  % enforce certain dynamics depending on the current transition (region activation

Ntrans = length(aut.trans);
Nmodes = length(aut.q);

% load map
%regFile_complexMap
regFile_boxPushing
% regFile_twoRegions
% aut.f = {1, 1};

% reg(1) = region(vReg{1});
% reg(2) = region(vReg{2});
% reg(3) = region(vReg{3});
% regBnd = region(vBnd{1});

% Model parameters - Unicycle
polyMdlFun{1} = @(t,x,u) CreateKinematicsPoly(t,x,u);
mdlFun{1} = @(t,x,options,p1,p2) CreateKinematicsNLWayptCtrl2(t,x,options,p1,p2);
ctrlFun{1} = @(x,options,p1,p2) CreateKinematicsNLWayptCtrl2_ctrl(x,options,p1,p2);
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

% Model parameters - Holonomic Robot
polyMdlFun{2} = @(t,x,u) HolonomicPoly(t,x,u);
mdlFun{2} = @(t,x,options,p1,p2) HolonomicWayptCtrl(t,x,options,p1,p2);
ctrlFun{2} = @(x,options,p1,p2) HolonomicWayptCtrl_ctrl(x,options,p1,p2);
sysparams(2) = sysparams(1);

% Model parameters - Dubins Car with fwd velocity as a parameter
polyMdlFun{3} = @(t,x,u) CreateKinematicsVelAsParamPoly(t,x,u);
mdlFun{3} = @(t,x,options,p1,p2) CreateKinematicsVelAsParamNLWayptCtrl2(t,x,options,p1,p2);
ctrlFun{3} = @(x,options,p1,p2) CreateKinematicsNLWayptCtrl2_ctrl(x,options,p1,p2);
sysparams(3) = sysparams(1);
sysparams(3).n = 4; % number of states
sysparams(3).isCyclic = [0 0 1 0]';
sysparams(3).limsNonRegState = [-pi 4; pi 4];

sys(1) = systemdynamics(polyMdlFun{1},mdlFun{1},ctrlFun{1},sysparams(1).H,sysparams(1));
sys(2) = systemdynamics(polyMdlFun{2},mdlFun{2},ctrlFun{2},sysparams(2).H,sysparams(2));
sys(3) = systemdynamics(polyMdlFun{3},mdlFun{3},ctrlFun{3},sysparams(3).H,sysparams(3));

% Feedback controller params
options.ctrloptions_trans.Q = 1*diag([0.01 0.01 0.01 1e-8]); %0.01*eye(sysparams.n);
% options.ctrloptions_trans.R = 0.007;
options.ctrloptions_trans.R = 0.05;
options.ctrloptions_trans.Qf = 1*eye(max([sysparams.n]));

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

trans = vertcat(aut.trans{:});

% for imode = 1:Nmodes
%     qReg = aut.q{imode};
%     qCover{imode} = getCoverPts(vReg,{qReg},1,options.Ncover,sysparams.H,n,sysparams.limsNonRegState);
% end

figure(5)
plot(reg)
axis equal
figure(500)
plot(reg)
hold on
plot(regDefl)
axis equal
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
    ac_inward{iModeToPatch} = [];
    for itrans = find(trans(:,1)==iModeToPatch)'
        j = j+1;
        ac_trans{j} = [];
    end
end

%%
joinedModes = false*ones(Nmodes,1);  % keep a list of modes which need joining
patchedModes = false*ones(Nmodes,1);  % keep a list of modes which need joining

%%
for iModeToGo = 1:NmodesReach
    
    patchedAndUnvistedModesToPatch = find(~(patchedModes(1:iModeToGo)));
    counter = 0;
    
    for iPatchTry = 1:20  % max number of tries
        for iModeToPatch =  patchedAndUnvistedModesToPatch'
            
            %% Reach Operation
            fprintf(fid,'%f : Attempting to patch and join mode %d.\n',toc,iModeToPatch);
            reachIncomplete = true;
            % while reachIncomplete
            [ac_tmp,lastTrans] = reachOp_new(sys,reg,regDefl,regBnd,aut,ac_trans,iModeToPatch,options);
            if ~isempty(lastTrans)
                disp('Reach was unsuccessful. Repeating the Reach.')
                reachIncomplete = true;
                joinedModes(iModeToPatch) = false;
                patchedModes(iModeToPatch) = false;
                counter = counter+1;
            else
                disp('Reach was successful.')
                counter = 0;
                reachIncomplete = false;
                joinedModes(iModeToPatch) = false;
                patchedModes(iModeToPatch) = true;
                tmp = find(trans(:,2)==iModeToPatch)';
                if length(tmp) == 1
                    if ~isempty(ac_trans{tmp})  % we have just joined the funnel from the previous transition
                        joinedModes(iModeToPatch) = true;
                    end
                end
                
            end
            
            iModeToGo
            patchedModes
            joinedModes
            iModeToPatch
            % end
            toc
            
            if ~patchedModes(iModeToPatch)
                iPreModeToPatch = [];
                if ~isempty([ac_trans{lastTrans}]) && counter > 10
                    counter = 0;
                    iPreModeToPatch = trans(lastTrans,1);
                    iPreModeToPatch
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
                for i = find(trans(:,1)==iModeToPatch)'  % TODO: ordering is preserved, as this agrees with what is in reachOp_new, but need to make more robust
                    j = j+1;
                    
                    % If there are two or more outgoing transitions for this
                    % mode, find a bounding polytope for the non-reactively
                    % composable part
                    % TODO: function...
%                     if sum(trans(:,1) == iModeToPatch) > 1
%                         
%                         % Find a new polytope for the region
%                         buildNewRegion
%                         
%                         % Update the Region file
%                         addNewRegionToFile
%                         
%                         % Update the Structuredslugs file
%                         modifySpecForNewRegion
%                     end
                    
                    ac_trans{i} = ac_tmp(j);
                    i
                    ac_trans
                end
            end
        end
        
        %% Sequence Operation
        ac_tmp = sequenceOp_new(sys,reg,regDefl,regBnd,aut,ac_trans,iModeToPatch,options);
        ac_inward{iModeToPatch} = ac_tmp;
        fprintf(fid,'%f : Finished Sequence operation for mode %d.\n',toc,iModeToPatch);
        
        if ~patchedModes(iModeToGo)
            % update the list of modes which have already been joined but are affected by the patch and concatenate with the list of unvisited modes
            patchedAndUnvistedModesToPatch = find(~(patchedModes(iModeToPatch:iModeToGo))) + (iModeToPatch-1);  % everything left in the todo list
            tmp = setdiff(trans(trans(:,1)==iModeToPatch,2),find(~joinedModes));   % include everything that that *could be* affected by re-patching the current mode; in particular, those that have been joined.
            patchedAndUnvistedModesToPatch = [setdiff(tmp,iModeToGo+1:Nmodes); patchedAndUnvistedModesToPatch];  % every successor to the current mode which had been joined
            if ~isempty(iPreModeToPatch)
                tmp = setdiff(trans(trans(:,1)==iPreModeToPatch,2),find(~joinedModes));  % include everything that that *could be* affected by re-patching the pre mode; in particular, those that have been joined.  
                %TODO: we can relax by assuming no effect, then checking for containment
                patchedAndUnvistedModesToPatch = [setdiff(tmp,iModeToGo+1:Nmodes); patchedAndUnvistedModesToPatch];  % do not flag modes which we know haven't yet been visited
                patchedAndUnvistedModesToPatch = [setdiff(iPreModeToPatch,iModeToGo+1:Nmodes); patchedAndUnvistedModesToPatch];  % don't forget the pre-mode!
            end
            
            patchedAndUnvistedModesToPatch = unique(patchedAndUnvistedModesToPatch);
            patchedAndUnvistedModesToPatch
            patchedModes = [~ismember(1:iModeToGo,patchedAndUnvistedModesToPatch)'; patchedModes(iModeToGo+1:Nmodes)];
            patchedModes
            
            %%
            if ~isempty(iPreModeToPatch)  % if join failed, also patch the pre mode
                disp(['... Now patching the pre-mode (mode ',num2str(iPreModeToPatch),') ...'])
                toc
                fprintf(fid,'%f : Now Patching the pre-mode %d to the failed mode (%d).\n',toc,iPreModeToPatch,iModeToPatch);
                
                [ac_tmp,lastTrans] = reachOp_new(sys,reg,regDefl,regBnd,aut,ac_trans,iPreModeToPatch,options);
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
                        if ~isempty(ac_trans{tmp})  % we have just joined the funnel from the previous transition
                            joinedModes(iPreModeToPatch) = true;
                        end
                    end
                end
                toc
                
                iModeToGo
                patchedModes
                joinedModes
                iModeToPatch
                
                if ~patchedModes(iPreModeToPatch)
                    disp(['... Failed to join mode ',num2str(iPreModeToPatch),'.'])
                    toc
                    fprintf(fid,'%f : Failed to perform Reach for mode %d. Quitting.\n',toc,iPreModeToPatch);
                    error('Failed to perform Reach');
                else
                    toc
                    patchedAndUnvistedModesToPatch = setdiff(patchedAndUnvistedModesToPatch, iPreModeToPatch);
                    fprintf(fid,'%f : Finished Reach operation for mode %d.\n',toc,iPreModeToPatch);
                    j = 0;
                    for i = find(trans(:,1)==iPreModeToPatch)'
                        j = j+1;
                        ac_trans{i} = ac_tmp(j);
                        i
                        ac_trans
                    end
                end
                fprintf(fid,'%f : ... Finished Reach operation for pre-mode %d.\n',toc,iPreModeToPatch);
                ac_tmp = sequenceOp_new(sys,reg,regDefl,regBnd,aut,ac_trans,iPreModeToPatch,options);
                ac_inward{iModeToPatch} = ac_tmp;
                fprintf(fid,'%f : ... Finished Sequence operation for pre-mode %d.\n',toc,iPreModeToPatch);
            end
        end
        if all(patchedModes(1:iModeToGo))
            disp(['finished ',num2str(iModeToGo)]);
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
            [ac_tmp,lastTrans] = joinOp_new(sys,reg,regDefl,regBnd,aut,ac_trans,iModeToJoin,options);
            if ~isempty(lastTrans)
                disp('Join was unsuccessful. Repeating.')
                joinIncomplete = true;
                joinedModes(iModeToJoin) = false;
            else
                disp('Join was successful.')
                joinIncomplete = false;
                joinedModes(iModeToJoin) = true;
                tmp = find(trans(:,2)==iModeToJoin)';
                j=0;
                i = find(trans(:,1)==iModeToJoin)';
                ac_inward{i} = [ac_inward{i}; ac_tmp];
                
%                 for i = find(trans(:,1)==iModeToJoin)'  % TODO: ordering is preserved, as this agrees with what is in reachOp_new, but need to make more robust
%                     j = j+1;
%                     ac_inward{i} = ac_tmp(j);
%                 end
%                 if length(tmp) == 1
%                     if ~isempty(ac_inward{tmp})  % we have just joined the funnel from the previous transition
%                         joinedModes(iModeToJoin) = true;
%                     end
%                 end
            end
            
            
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
%         if ~joinedModes(imode)
%             % update the list of modes which have already been joined but are affected by the patch and concatenate with the list of unvisited modes
%             patchedAndUnvistedModesToJoin = find(~(joinedModes(iModeToJoin:imode))) + (iModeToJoin-1);  % everything left on the todo list
%             patchedAndUnvistedModesToJoin = [setdiff(trans(trans(:,1)==iModeToJoin,2),imode+1:Nmodes); patchedAndUnvistedModesToJoin];  % every successor to the current mode which had been joined
%             if ~joinedModes(iModeToJoin)
%                 patchedAndUnvistedModesToJoin = [setdiff(trans(trans(:,1)==iPreModeToJoin,2),imode+1:Nmodes); patchedAndUnvistedModesToJoin];
%                 patchedAndUnvistedModesToJoin = [setdiff(iPreModeToJoin,imode+1:Nmodes); patchedAndUnvistedModesToJoin];
%             end
%             
%             patchedAndUnvistedModesToJoin = unique(patchedAndUnvistedModesToJoin);
%             patchedAndUnvistedModesToJoin
%             joinedModes = [~ismember(1:imode,patchedAndUnvistedModesToJoin)'; joinedModes(imode+1:Nmodes)];
%             
%             % Patch Operation
%             disp(['... Patching the failed current mode (mode ',num2str(iModeToJoin),') ...'])
%             toc
%             fprintf(fid,'%f : Patching the failed current mode %d.\n',toc,iModeToJoin);
%             iModeToPatch = iModeToJoin;
%             reachOp
%             toc
%             fprintf(fid,'%f : ... Finished Reach operation for mode %d.\n',toc,iModeToJoin);
%             ac_inward = sequenceOp_new(sys,reg,regBnd,aut,iPreModeToPatch,options);
%             fprintf(fid,'%f : ... Finished Sequence operation for mode %d.\n',toc,iModeToJoin);
%             
%             if ~joinedModes(iModeToJoin)  % if join failed, also patch the pre mode
%                 disp(['... Now patching the pre-mode (mode ',num2str(iPreModeToJoin),') ...'])
%                 toc
%                 fprintf(fid,'%f : Now Patching the pre-mode %d to the failed mode (%d).\n',toc,iPreModeToJoin,iModeToJoin);
%                 iModeToPatch = iPreModeToJoin;
%                 reachOp
%                 toc
%                 fprintf(fid,'%f : ... Finished Reach operation for pre-mode %d.\n',toc,iPreModeToJoin);
%                 sequenceOp
%                 fprintf(fid,'%f : ... Finished Sequence operation for pre-mode %d.\n',toc,iPreModeToJoin);
%             end
%         end
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

for imode = 1:NmodesReach
    if sum(trans(:,1) == imode) > 1 
        
        [ac_tmp, bc_tmp, lastTrans, existingReg, newRegArray, reg] = reactiveJoinOp_new1(sys,reg,regDefl,regBnd,aut,ac_trans,[],imode,calibMatrix,options);
        
        if ~isempty(lastTrans)
            disp('Reactive Join was unsuccessful. ')
            joinIncomplete = true;
            joinedModes(imode) = false;
        else
            disp('Reactive Join was successful.')
            joinIncomplete = false;
            joinedModes(imode) = false;
            tmp = find(trans(:,2)==imode)';
            j=0;
            % i = find(trans(:,1)==imode)';
            for j = 1:length(ac_tmp)
                ac_react{imode}(j) = ac_tmp(j);
                bc_react{imode}(j) = bc_tmp(j);
            end
            
            ac_inward{imode} = ac_react{imode};
            
%             for i = find(trans(:,1)==iModeToJoin)'  % TODO: ordering is preserved, as this agrees with what is in reachOp_new, but need to make more robust
%                 j = j+1;
%                 ac_react{i} = ac_tmp(j);
%             end
%             if length(tmp) == 1
%                 if ~isempty(ac_react{tmp})  % we have just joined the funnel from the previous transition
%                     joinedModes(iModeToJoin) = true;
%                 end
%             end
        end
        
        iModeToPatch = imode;
        for i = find(trans(:,1)==imode)'
            % If there are two or more outgoing transitions for this mode, find a bounding polytope for the non-reactively composable part
            % TODO: function...
            % Find a new polytope for the region
            %buildNewRegion
            
%             if ~isempty(newReg)
%                 % Update the Region file
%                 addNewRegionToFile
%                 
%                 % Update the Structuredslugs file
%                 modifySpecForNewRegion
%             end
        end
    end
end


% for iModeToReactJoin = unique(trans(isReactive,1))'
%     disp(['... Attempting to reactively join the mode ',num2str(iModeToReactJoin),'.'])
%     % Reactive Join Operation
%     for iReactTry = 1:10
%         reactivejoinOp
%         if reactJoinedModes(iModeToReactJoin), break, end
%         toc
%         fprintf(fid,'%f : ... Reactive Join operation failed for mode %d.  Retrying... (%d)\n',toc,iModeToReactJoin,iReactTry);
%     end
%     toc
%     fprintf(fid,'%f : Finished Reactive Join operation for mode %d.\n',toc,iModeToReactJoin);
% end
% toc
% fprintf(fid,'%f : Finished Reactive Join operation for all modes. Success!\n',toc);


%%
% clk = fix(clock);
% if strcmp(modelType,'unicycle')
%     eval(['save antsy5_uni_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6))])
% elseif strcmp(modelType,'car')
%     eval(['save antsy5_car_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6))])
% end

toc
fclose(fid);
