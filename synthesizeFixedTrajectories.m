function [ac_trans] = synthesizeFixedTrajectories(sys, filePath, configAndProblemDomainName, stateTraject, inputTraject, options)
%
% synthesizeFixedTrajectories:
% Given an automaton, collection of regions and trajectories, this function
% synthesizes atomic controllers collecting time-varying LQR controllers
% and funnels for each controller.
%
% Inputs:
%   - sys: The DrakeSystem model containing the equations of motion of the
%   system.
%   - filePath: A string naming the file path where all the files
%   pertaining to the problem domain and configuration exist.  e.g.
%   '/path/to/reasyns/examples/box_pushing/'.
%   - configAndProblemDomainName: A string with the uniform name given to
%   all the files in the 'filePath' directory.  The required files are:
%      <configAndProblemDomainName>.m
%      <configAndProblemDomainName>.aut
%      <configAndProblemDomainName>.regions
%      <configAndProblemDomainName>.spec
%      configs/<configAndProblemDomainName>.config
%   - stateTraject: An p-by-q structure of Trajects containing a trajectory
%   of system state variables for predecessor mode p to the successor mode
%   q in terms of the state designations in the automaton structure `aut`.
%   - stateTraject: An p-by-q structure of Trajects containing a trajectory
%   of system inputs for predecessor mode p to the successor mode q in
%   terms of the state designations in the automaton structure `aut`.
%   - options: An options structure containing the options required to
%   build the atomic controllers and funnels.  See
%   "plane_with_trajectory.m" for details of this structure's contents.
%
% Output:
%   - ac_trans: A p-by-q structure of QuadraticACs containing the
%   trajectories, controllers and funnels for each controller.
%

warning off
format compact
close all

fileName = [filePath,'/',configAndProblemDomainName];
if ~isdir([filePath,'/reasyns_log']), mkdir([filePath,'/reasyns_log']); end
if ~isdir([filePath,'/reasyns_controllers']), mkdir([filePath,'/reasyns_controllers']); end

% Load the ordered set of regions.
[reg, regDefl, ~] = processRegFile(filePath, configAndProblemDomainName, options);
%processRegFile(filePath, configAndProblemDomainName, options, 1/40)  
% without the calibration matrix -- instead provide a scalar defining the relationship between pixels to lab units

% Load the automaton.
aut = processAutFileFastSlow([filePath, '/', configAndProblemDomainName, '.aut'],reg);
trans = vertcat(aut.trans{:});

Ntrans = length(aut.trans);
Nmodes = max([aut.state{:}]);

aut.f = cell(1,Ntrans);  % enforce certain dynamics depending on the current transition (region activation) 
[aut.f{:}] = deal(1);

options.maxFunnelsTrans(1:Nmodes) = 1;

% Set the clock (used for file naming).
clk = fix(clock);

% Prepare a log file.
if options.doLog
    fid = fopen([filePath,'/reasyns_log/reasyns_log_', date, '_', num2str(clk(4)), num2str(clk(5)), num2str(clk(6)), '.txt'], 'w');
end

% Make copies of files that we may be modifying later on.
resetAllInputFiles

% Choose SDP solver.
if checkDependency('mosek')
    options.solver = @spot_mosek;
    options.solver_name = 'mosek';
elseif checkDependency('sedumi')
    options.solver = @spot_sedumi;
    options.solver_name = 'sedumi';
else
    error('synthesizeFixedTrajectories: Please install either MOSEK or SeDuMi and try again.');
end

figure(1)
clf
plot(reg)
axis equal
hold on

figure(3)
clf
plot(reg)
hold on
plot(regDefl)
axis equal
% axis([min(regBnd{1}.v(:,1)) max(regBnd{1}.v(:,1)) min(regBnd{1}.v(:,2)) max(regBnd{1}.v(:,2))])


%%  Main

tic


%%
% NmodesReach = Nmodes;  % Commented out for now since we are only this
                         % function with a single trajectory.
NmodesReach = 4;

tic
toc
if options.doLog, fprintf(fid, '%f : Starting ....\n', toc); end
j = 0;
for statePre = 1:NmodesReach
    ac_inward{statePre} = [];  % N.B. ac_inward remains empty in this synthesis function.
    for itrans = find(trans(:,1) == statePre)'
        j = j+1;
        ac_trans{j} = [];
    end
end

%%
joinedStates = false*ones(Nmodes,1);  % keep a list of modes which need joining
patchedStates = false*ones(Nmodes,1);  % keep a list of modes which need joining

%%
for indexToGo = 1:NmodesReach
    
    patchedAndUnvistedStatesToPatch = find(~(patchedStates(1:indexToGo)));
    counter = 0;
    
    for statePre =  patchedAndUnvistedStatesToPatch'
        
        if options.doLog, fprintf(fid,'%f : Attempting to compute funnels for State %d.\n', toc, statePre); end
        % while reachIncomplete
        
        [indexTransVect, indexPostVect] = findTransitionsWithNonRepeatingRegions(aut, statePre);
        indexTransVect = indexTransVect';
        
        % remove any transitions, for whose region pairs, a funnel has already been created.
        newIndexTransVect = indexTransVect;
        for i = 1:length(ac_trans)
            if isempty(ac_trans{i}), continue, end
            preState = ac_trans{i}.pre;
            ac_trans{i}.post
            for postState = ac_trans{i}.post
                for indexTrans = indexTransVect'
                    if (aut.label{preState} ~= aut.label{aut.trans{indexTrans}(1)}), continue, end
                    if (aut.label{postState} ~= aut.label{aut.trans{indexTrans}(2)}), continue, end
                    newIndexTransVect = setdiff(newIndexTransVect, indexTrans);
                end
            end
        end
        
        errorMsg = [];
        if ~isempty(newIndexTransVect)
            for indexTrans = newIndexTransVect'
                statePost = trans(indexTrans,2);
                
                x0 = stateTraject{statePre, statePost};
                u0 = inputTraject{statePre, statePost};
                [ac_tmp, errorMsg] = computeFunnel(sys, x0, u0, reg, aut, ac_trans, statePre, statePost, options);
                
                indexSet = vertcat(aut.label{indexPostVect}) == vertcat(aut.label{vertcat(aut.state{:}) == statePre});
                statePostVect = vertcat(aut.state{indexSet});
                for i = 1:length(ac_tmp)
                    ac_tmp(i).setTransition(statePre, statePostVect);
                end
                ac_trans{indexTrans} = ac_tmp;
            end
        end
        
        if ~isempty(errorMsg)
            disp(['synthesizeFixedTrajectories: Transition funnel computation was unsuccessful for State ',num2str(statePre),'. Repeating.'])
            joinedStates(statePre) = false;
            patchedStates(statePre) = false;
            counter = counter+1;
        else
            disp(['synthesizeFixedTrajectories: Transition funnel computation was successful for State ',num2str(statePre),'.'])
            counter = 0;
            joinedStates(statePre) = false;
            patchedStates(statePre) = true;
            
            joinedStates(statePre) = false;  % If 'true' then we have just joined the funnel from the previous transition
            if false  % 'false' here means we always just assume the worst case and that no funnels are composed.
                % TODO: fix the following to flag the state as 'joined' only if it passes a check for composition
                for jpre = 1:length(ac_trans)
                    if isempty(ac_trans{jpre}), continue, end
                    if ~any(ac_trans{jpre}.post == statePre), continue, end
                    tmp = ac_trans{jpre}.ellipsoid;
                    ellToCompose = tmp(end);
                    for jpost = 1:length(ac_trans)
                        if isempty(ac_trans{jpost}), continue, end
                        if ~any(ac_trans{jpost}.pre == statePost), continue, end
                        if ~ac_trans{jpost}.funnelContainsEllipsoid(sys, ellToCompose), continue, end
                        joinedStates(statePre) = true;
                    end
                end
            end
            
        end
        
        disp(toc);
        
        if patchedStates(statePre) , continue, end
        
        iPreStateToPatch = [];
        if ~isempty([ac_trans{lastTrans}]) && counter > 10
            counter = 0;
            iPreStateToPatch = trans(lastTrans,1);
        end
        disp(['synthesizeFixedTrajectories: Failed to create a transition funnel for State ',num2str(statePre),'. '])
        toc
        if options.doLog, fprintf(fid,'%f : Failed to create a transition funnel for State  %d. Re-computing the transition funnel.\n', toc, statePre); end
        break
    end
end
toc


%%
toc
if options.doLog, fclose(fid); end

% Save the controller data to file.
log_string = [filePath,'/reasyns_controllers/reasyns_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6))];
saveReasynsProblemInstance(log_string, sys, aut, ac_inward, ac_trans)
%eval(['save(''',filePath,'/reasyns_controllers/reasyns_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6)),''');'])

