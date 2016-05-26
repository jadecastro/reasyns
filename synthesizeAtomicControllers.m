function synthesizeAtomicControllers(sys, filePath, configAndProblemDomainName, options)
%
% ReaSyNS - REActive SYnthesis for Nonlinear Systems
%

warning off
format compact
close all

fileName = [filePath,'/',configAndProblemDomainName];
if ~isdir([filePath,'/reasyns_log']), mkdir([filePath,'/reasyns_log']); end
if ~isdir([filePath,'/reasyns_controllers']), mkdir([filePath,'/reasyns_controllers']); end

% Load the ordered set of regions
[reg,regDefl,regBnd] = processRegFile(filePath, configAndProblemDomainName, options);
%processRegFile(filePath, configAndProblemDomainName, options, 1/40)  % without the calibration matrix -- instead provide a scalar defining the relationship between pixels to lab units

% Load the automaton
aut = processAutFileFastSlow([filePath,'/',configAndProblemDomainName,'.aut'],reg);
trans = vertcat(aut.trans{:});

Ntrans = length(aut.trans);
Nmodes = max([aut.state{:}]);

aut.f = cell(1,Ntrans);  % enforce certain dynamics depending on the current transition (region activation) 
[aut.f{:}] = deal(1);

options.maxFunnelsTrans(1:Nmodes) = options.maxFunnelsTrans;
options.maxFunnelsInward(1:Nmodes) = options.maxFunnelsInward;
options.maxFunnelsReactJoin(1:Nmodes) = options.maxFunnelsReactJoin;

% Prepare a log file
clk = fix(clock);
fid = fopen([filePath,'/reasyns_log/reasyns_log_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6)),'.txt'],'w');

% make copies of files that we may be modifying later on
resetAllInputFiles


% Choose SDP solver
if checkDependency('mosek')
    options.solver = @spot_mosek;
    options.solver_name = 'mosek';
elseif checkDependency('sedumi')
    options.solver = @spot_sedumi;
    options.solver_name = 'sedumi';
else
    error('Please install either MOSEK or SeDuMi and try again.');
end

figure(1)
clf
plot(reg)
axis equal
hold on

figure(5)
clf
plot(reg)
axis equal
hold on

figure(500)
clf
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
for statePre = 1:NmodesReach
    ac_inward{statePre} = [];
    for itrans = find(trans(:,1)==statePre)'
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
    
    for iPatchTry = 1:20  % max number of tries
        for statePre =  patchedAndUnvistedStatesToPatch'
    

            %% Reach Operation
            fprintf(fid,'%f : Attempting to compute transition funnels for State %d.\n',toc,statePre);
            reachIncomplete = true;
            % while reachIncomplete
            
            [indexTransVect, indexPostVect] = findTransitionsWithNonRepeatingRegions(aut,statePre);
            indexTransVect = indexTransVect';
            
            % remove any transitions, for whose region pairs, a funnel has already been created. 
            newIndexTransVect = indexTransVect;
            for i = 1:length(ac_trans)
                if ~isempty(ac_trans{i})
                    preState = ac_trans{i}.pre;
                    for postState = ac_trans{i}.post
                        for indexTrans = indexTransVect'
                            if (aut.label{preState} == aut.label{aut.trans{indexTrans}(1)} && aut.label{postState} == aut.label{aut.trans{indexTrans}(2)})
                                newIndexTransVect = setdiff(newIndexTransVect,indexTrans);
                            end
                        end
                    end
                end
            end
            
            lastTrans = [];
            for indexTrans = newIndexTransVect'
                statePost = trans(indexTrans,2);
                
                [ac_tmp,lastTrans] = computeTransitionFunnel(sys,reg,regDefl,regBnd,aut,ac_trans,statePre,statePost,options);
                
                statePostVect = vertcat(aut.state{vertcat(aut.label{indexPostVect}) == vertcat(aut.label{vertcat(aut.state{:}) == statePre})});
                for i = 1:length(ac_tmp)
                    ac_tmp(i).setTransition(statePre,statePostVect);
                end
                ac_trans{indexTrans} = ac_tmp;
            end
            
            if ~isempty(lastTrans)
                disp(['Transition funnel computation was unsuccessful for State ',num2str(statePre),'. Repeating.'])
                reachIncomplete = true;
                joinedStates(statePre) = false;
                patchedStates(statePre) = false;
                counter = counter+1;
            else
                disp(['Transition funnel computation was successful for State ',num2str(statePre),'.'])
                counter = 0;
                reachIncomplete = false;
                joinedStates(statePre) = false;
                patchedStates(statePre) = true;
   
                joinedStates(statePre) = false;  % 'true' if we have just joined the funnel from the previous transition
                if false  % 'false' means we will assume nothing has joined
                    % TODO: fix the following to flag the state as 'joined' only if it passes a check for composition
                    % TODO: cleanup
                    for jpre = 1:length(ac_trans)
                        if ~isempty(ac_trans{jpre})
                            if any(ac_trans{jpre}.post == statePre)
                                tmp = ac_trans{jpre}.ellipsoid;
                                ellToCompose = tmp(end);
                                for jpost = 1:length(ac_trans)
                                    if ~isempty(ac_trans{jpost})
                                        if any(ac_trans{jpost}.pre == statePost) && ac_trans{jpost}.funnelContainsEllipsoid(sys,ellToCompose)
                                            joinedStates(statePre) = true;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
            end

%             patchedStates'
%             joinedStates'
%             statePre
            % end
            toc
            
            if ~patchedStates(statePre)
                iPreStateToPatch = [];
                if ~isempty([ac_trans{lastTrans}]) && counter > 10
                    counter = 0;
                    iPreStateToPatch = trans(lastTrans,1);
%                     iPreStateToPatch
                end
                disp(['... Failed to create a transition funnel for State ',num2str(statePre),'. '])
                toc
                fprintf(fid,'%f : Failed to create a transition funnel for State  %d. Re-computing the transition funnel.\n',toc,statePre);
                break
            end
        end
        
        if all(patchedStates(1:indexToGo))
%             disp(['finished ',num2str(indexToGo)]);
            break
        end
    end
end
toc


%%
indexState = 0;
for statePost = horzcat(aut.state{:})
    
    foundAc = false;
    for indexTrans = 1:length(aut.trans)
        if ~isempty(ac_trans{indexTrans})
            if any(ac_trans{indexTrans}.pre == statePost), foundAc = true; break; end
        end
    end
    
    if ~foundAc
        indexState = indexState + 1;
        joinedStates(indexState) = 1;
    else
        for iJoinTry = 1:10  % max number of tries
            disp(['... Attempting to create inward joining funnels for State ',num2str(statePost)])
            toc
            fprintf(fid,'%f : Attempting to create inward joining funnels for State %d.\n',toc,statePost);
            
            % Join Operation
            
            [ac_tmp,lastTrans] = computeInwardJoiningFunnel(sys,reg,regDefl,regBnd,aut,ac_trans,statePost,options);
            
            for i = 1:length(ac_tmp)
                ac_tmp(i).setTransition(statePost);
            end
            
            if ~isempty(lastTrans)
                disp(['Inward joining funnel computation was unsuccessful for State ',num2str(statePost),'. Repeating.'])
                joinIncomplete = true;
                joinedStates(indexState) = false;
            else
                disp(['Inward joining funnel computation was successful for State ',num2str(statePost),'.'])
                joinIncomplete = false;
                joinedStates(indexState) = true;
                ac_inward{indexState} = [ac_inward{indexState}; ac_tmp];
            end
            
            
            if ~joinedStates(indexState)
                disp(['... Failed to create inward joining funnels for State ',num2str(statePost),'.'])
                toc
                fprintf(fid,'%f : Failed to create inward joining funnels for State %d.\n',toc,statePost);
                break
            else
                toc
                fprintf(fid,'%f : Successfully computed inward joining funnels for State %d.\n',toc,statePost);
            end
            
            if all(joinedStates(1:indexState))
                break
            end
        end
    end
end


%% Prompt the user to continue generating additional inward funnels

doInwardReactive = false;
userResponse = [];
while ~isstr(userResponse) || ~all(ismember(lower(userResponse),['y','n']))
    userResponse = input('Transition and inward joining funnel generation completed.  Create reactive inward funnels? [y/N]: ','s');
    if isempty(userResponse), userResponse = 'n'; end
end
if lower(userResponse) == 'y', doInwardReactive = true; end


%%
if doInwardReactive
    reactJoinedStates = false*ones(Nmodes,1);
    
    for state = horzcat(aut.state{:})
        if sum(trans(:,1) == state) > 1
            
            [ac_tmp, bc_tmp, lastTrans, existingReg, newRegArray, reg] = computeInwardReactiveFunnel(sys,reg,regDefl,regBnd,aut,ac_trans,[],state,options,fileName);
            
            for i = 1:length(ac_tmp)
                ac_tmp(i).setTransition(state);
            end
            for i = 1:length(bc_tmp)
                bc_tmp(i).setTransition(state);
            end
            
            if ~isempty(lastTrans)
                disp(['Computation of inward reactive funnels was unsuccessful for State ',num2str(state),'. '])
                joinIncomplete = true;
                joinedStates(state) = false;
            else
                disp(['Computation of inward reactive funnels was successful for State ',num2str(state),'.'])
                joinIncomplete = false;
                joinedStates(state) = false;
                tmp = find(trans(:,2)==state)';
                j=0;
                % i = find(trans(:,1)==state)';
                for j = 1:length(ac_tmp)
                    ac_react{state}(j) = ac_tmp(j);
                    bc_react{state}(j) = bc_tmp(j);
                end
                
                ac_inward{state} = [ac_inward{state}; ac_react{state}];
                
            end
        end
    end
end

%%
toc
fclose(fid);

saveReasynsProblemInstance([filePath,'/reasyns_controllers/reasyns_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6))],sys,aut,ac_inward,ac_trans)
%eval(['save(''',filePath,'/reasyns_controllers/reasyns_',date,'_',num2str(clk(4)),num2str(clk(5)),num2str(clk(6)),''');'])

