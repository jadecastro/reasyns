function verification(sys, filePath, configAndProblemDomainName, options)
%
% ReaSyNS - REActive SYnthesis for Nonlinear Systems
%

warning off
format compact
close all

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


%%
NmodesReach = Nmodes;
% NmodesReach = 1;

j = 0;
for statePre = 1:NmodesReach
    ac_inward{statePre} = [];
    for itrans = find(trans(:,1)==statePre)'
        j = j+1;
        ac_trans{j} = [];
    end
end

%%
patchedStates = false*ones(Nmodes,1);  % keep a list of modes which need joining

%%
for indexToGo = 1:NmodesReach
    
    patchedAndUnvistedStatesToPatch = find(~(patchedStates(1:indexToGo)));
    counter = 0;
    
    for iPatchTry = 1:20  % max number of tries
        for statePre =  patchedAndUnvistedStatesToPatch'
    

            %% Reach Operation
            % while reachIncomplete
            
            [indexTransVect, indexPostVect] = findTransitionsWithNonRepeatingRegions(aut,statePre);
            indexTransVect = indexTransVect';
            
            % remove any transitions, for whose region pairs, a funnel has already been created. 
            if ~isempty(indexTransVect)
                newIndexTransVect = indexTransVect;
                for i = 1:length(ac_trans)
                    if ~isempty(ac_trans{i})
                        preState = ac_trans{i}.pre;
                        ac_trans{i}.post
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
                if ~isempty(newIndexTransVect)
                    for indexTrans = newIndexTransVect'
                        failed = true;
                        while failed
                            statePost = trans(indexTrans,2);
                            try
                                [ac_tmp,lastTrans] = computeTransitionFunnel(sys,reg,regDefl,regBnd,aut,ac_trans,statePre,statePost,options);
                                failed = false;
                            catch ME
                                disp(ME.message)
                                ac_tmp = [];
                            end
                            statePostVect = vertcat(aut.state{vertcat(aut.label{indexPostVect}) == vertcat(aut.label{vertcat(aut.state{:}) == statePre})});
                            for i = 1:length(ac_tmp)
                                ac_tmp(i).setTransition(statePre,statePostVect);
                            end
                            ac_trans{indexTrans} = ac_tmp;
                        end
                    end
                end
                
                if ~isempty(lastTrans)
                    disp(['Transition funnel computation was unsuccessful for State ',num2str(statePre),'. Repeating.'])
                    patchedStates(statePre) = false;
                    counter = counter+1;
                else
                    disp(['Transition funnel computation was successful for State ',num2str(statePre),'.'])
                    counter = 0;
                    patchedStates(statePre) = true;
                end
                
            end

        end
        
        if all(patchedStates(1:indexToGo))
            break
        end
    end
end


