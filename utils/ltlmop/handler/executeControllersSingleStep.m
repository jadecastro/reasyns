function [vx,vy,w,errorMsg,acLastData] = executeControllersSingleStep(aut,sysArray,ac_trans,ac_inward,x,t,currReg,nextReg,acLastData)

persistent t_offset t_base tend ell t_trials

global u_sav x_sav x0_sav

delayTime = -0.1;

timeBaseFlag = false;  % if 'true', use time as a basis for choosing the TVLQR states; otherwise use a weighted Euclidean distance.
useLastAC = true;

%%%%%%%%%%%%%%%%
% Account for our manual shifting of the map to the left because of the "dead zone" in the Vicon field 
x(1) = x(1) + 0.0043;
x(2) = x(2) + 0.1629;
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% Account for our manual shifting of the map manually
x(1) = x(1) - 0.9774;
x(2) = x(2) + 0.0258;
%%%%%%%%%%%%%%%%

x(3) = x(3) + 0.2;  % theta bias needed to account for misalignment of the youBot's coordinate frame wrt. its true orientation

prevCtrl = false;

if isempty(t_base)
    t_base = clock;
end

vx = 0;  vy = 0;  w = 0.001;
errorMsg = 'nothing to see here...';

teval = 0.02;

try
    if size(x,2) > size(x,1), x = x'; end
    
    trans = vertcat(aut.trans{:});
    currStateVec = find([aut.q{:}]==currReg);
    nextStateVec = find([aut.q{:}]==nextReg);
    currTrans = find(ismember(trans(:,1),currStateVec) & ismember(trans(:,2),nextStateVec), 1);
    if isempty(currTrans)
        error('you have specified an invalid transition!')
    end
    currState = trans(currTrans,1);  nextState = trans(currTrans,2);
    
    % deal with the possibility of sys being a cell array of systems.
    if length(sysArray) > 1
        sys = sysArray(aut.f{currTrans});
    else
        sys = sysArray;
    end
    
    % identify the funnel we're in and get the ellipse index
    iTrans = false;
    iIn = false;
    if isinternal(ac_trans{currTrans}, x,'u',sys)
%     if isinternal(ac_trans{currTrans}, x,'u')
        iTrans = currTrans;
        ac = ac_trans{currTrans};
    end
    if iTrans == false && ~isempty(ac_inward)
        for funIdx = 1:length(ac_inward{currState})
            if isinternal(ac_inward{currState}(funIdx), x,'u',sys)
                %         if isinternal(ac_inward{currState}, x,'u')
                iIn = currState;
                ac = ac_inward{currState}(funIdx);
            end
        end
    end
    
    % if no funnels in the new region/transition OR have temporarily left a funnel, then keep activating the previous one.
    currStateOld = [];  nextStateOld = [];
    if ~iTrans && ~iIn
        if useLastAC
            % propagate the appropriate old states on to the next iteration
            if ~isempty(acLastData{5}) && ~isempty(acLastData{6})
                currStateOld = acLastData{5};  nextStateOld = acLastData{6};
            else
                currStateOld = acLastData{3};  nextStateOld = acLastData{4};
            end
            currTransOld = trans(:,1) == currStateOld & trans(:,2) == nextStateOld;
            if ~isempty(acLastData{1})
                ac = ac_trans{currTransOld};
            else
                ac = ac_inward{currStateOld};
            end
            iTrans = acLastData{1};  iIn = acLastData{2};
            prevCtrl = true;
            disp('WARNING: no funnels found! Using the previous controller.')
        else
            iTrans = currTrans;
            ac = ac_trans{currTrans};
            disp('WARNING: no funnels found! Forcing an (unverified) controller.')
        end
    end
    
    
    % determine the new time offset if something has changed
    if timeBaseFlag
        acData = {iTrans, iIn, currState, nextState, currStateOld, nextStateOld};
        for i = 1:length(acData)-2,
            cmp(i) = acData{i}==acLastData{i};
        end
        cmp
        if (~all(cmp) || ~(all(ismember([acLastData{1:4}],[acData{1:4}])) && all([acLastData{3:4}] == [acData{3:4}]))) && ~prevCtrl
            t_trials = getTimeVec(ac.x0);
            tend = t_trials(end);
            t_base = clock
            ell = ellipsoid(ac);
            
            for i = 1:length(ell)
                if isinternal(ell(i),x,'u')
%                 if isinternal(ell(i),x,'u')
                    t_offset = t_trials(i);
                    break
                end
            end
        end
        
        acLastData = acData;
        
        if isempty(t)  % if time is not given, we compute it here
            t = etime(clock,t_base) + delayTime
        end
        t_offset
        teval = min(tend,t + t_offset);
    else
        acData = {iTrans, iIn, currState, nextState, currStateOld, nextStateOld};
        for i = 1:length(acData)-2,
            cmp(i) = acData{i}==acLastData{i};
        end
        cmp
        if ~all(cmp) || ~(all(ismember([acLastData{1:4}],[acData{1:4}])) && all([acLastData{3:4}] == [acData{3:4}]))
            t_base = clock
            ell = ellipsoid(ac);
            t_trials = getTimeVec(ac.x0);
        end
        acLastData = acData;
        
        minDelta = inf;
        weights = [1;1;0.2];
        for i = 2:3:length(ell)
            xtmp = double(ac.x0, t_trials(i));
            testMinDist = min([
                norm(weights.*(xtmp(1:3) - (x(1:3)+[0 0 2*pi]'))) 
                norm(weights.*(xtmp(1:3) - x(1:3))) 
                norm(weights.*(xtmp(1:3) - (x(1:3)-[0 0 2*pi]')))]);
            if testMinDist < minDelta
                teval = t_trials(i);
                minDelta = testMinDist;
            end
        end
        if isempty(t)  % if time is not given, we compute it here
            t = etime(clock,t_base)
        end
    end
    
    % compute a command
    teval
    K = double(ac.K,teval)
    x0 = double(ac.x0,teval)
    u0 = double(ac.u0,teval);
    
    K = K(end-length(u0)+1:end,:);
    
    u_test = [(K*(x+[0 0 2*pi]' - x0)) (K*(x - x0)) (K*(x-[0 0 2*pi]' - x0))];
    u_idx = abs(u_test(end,:)) == min(abs(u_test(end,:)));
    u_ctrl = u_test(:,u_idx(1))
    u = u0 + 1*u_ctrl
%     u(2) = max(min(u(2),16),-16);
    
    if aut.f{currTrans} == 1  % we're using diff-drive dynamics
        vx = u(1);  w = max(min(u(2),0.2),-0.2);
    elseif aut.f{currTrans} == 2  % we're using holonomic-drive dynamics
        u_local = [cos(-x(3)) -sin(-x(3));sin(-x(3)) cos(-x(3))]*u(1:2)
        vx = 1*u_local(1);  vy = 1*u_local(2);  w = u(3);
    else
        error('unhandled system identifier!')
    end
    
    %u_sav = [u_sav; t u'];
    x_sav = [x_sav; t x'];
    x0_sav = [x0_sav; t x0']
    
catch ME
    rethrow(ME)
    errorMsg = ME.message;
end

