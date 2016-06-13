function [randPt] = getCenterRand(sys,reg,ac_trans,varargin)
% Draw a sample within the configuration space of a system while remaining
% within a given region. If a configuration vector is given, a sample is
% drawn from a Gaussian centered about that point. Otherwise, the Cartesian
% directions are sampled from a Gaussian centered about the centroid of the
% region and any other (bounded) directions are sampled uniformly.

n = sys.sysparams.n;
limsNonRegState = sys.sysparams.limsNonRegState;
stateLimits = sys.sysparams.stateLimits;
isCyclic = sys.sysparams.isCyclic;
Qrand = sys.sysparams.Qrand;

minNonRegStates = limsNonRegState(1,:);
maxNonRegStates = limsNonRegState(2,:);
minStates = stateLimits(1,:);
maxStates = stateLimits(2,:);

[idxRegStates,idxNonRegStates] = getRegNonRegStates(sys);

figure(90), plot(reg)
hold on

if isempty(varargin)
    tmpNonRegStates = minNonRegStates + (maxNonRegStates - minNonRegStates).*rand(1,length(idxNonRegStates));
else
    qCenter = varargin{1}';
    tmpNonRegStates = qCenter(idxNonRegStates);
end
    
if isempty(varargin)
    vReg = reg.v;
    maxReg = max(vReg);
    minReg = min(vReg);
    tmpRegStates = 1/2*(maxReg + minReg)' + Qrand(idxRegStates,idxRegStates)*randn(length(idxRegStates),1);
else
    tmpRegStates = qCenter(idxRegStates)' + Qrand(idxRegStates,idxRegStates)*randn(length(idxRegStates),1);
end

% If in obstacle or outside invariant, recompute its location in a deterministic fashion
for indx = 2:100
    testState = zeros(n,1);
    testState(idxRegStates) = tmpRegStates;
    testState(idxNonRegStates) = tmpNonRegStates;
    
%     figure(90), plot(tmpRegStates(1),tmpRegStates(2),'ro')
    
    isect = isinside(reg,sys,testState);
    if isect
        if ~isempty(ac_trans)
            res = true;
            for i = 1:length(ac_trans)
                res = res & isinternal(ac_trans(i),[tmpRegStates; tmpNonRegStates],'u');
            end
            if res
                break
            else
                tmpNonRegStates = minNonRegStates' + (maxNonRegStates - minNonRegStates)'.*rand(length(idxNonRegStates),1);
            end
        else
            break
        end
    end
    if isempty(varargin)
        tmpRegStates = 1/2*(maxReg + minReg)' + Qrand(idxRegStates,idxRegStates)*randn(length(idxRegStates),1);
    else
        tmpRegStates = qCenter(idxRegStates)' + Qrand(idxRegStates,idxRegStates)*randn(length(idxRegStates),1);
    end
end

if indx == 100
    
    tmpStates = minStates' + (maxStates - minStates)'.*rand(n,1);
    
    for indx = 2:100
        isect = isinside(reg,sys,tmpStates);
        if isect
            if ~isempty(ac_trans)
                res = true;
                for i = 1:length(ac_trans)
                    res = res & isinternal(ac_trans(i),tmpStates,'u');
                end
                if res
                    break
                end
            else
                break
            end
        end
        
        tmpStates = minStates' + (maxStates - minStates)'.*rand(n,1);
        tmpRegStates = sys.state2SEconfig([],tmpStates,[]);
        figure(90), plot(tmpRegStates(1),tmpRegStates(2),'ro')
    end
end

if indx == 100
    error('Can''t find a non-colliding point.')
end

randPt = zeros(1,n);
randPt(idxRegStates) = tmpRegStates; 
randPt(idxNonRegStates) = tmpNonRegStates;

if isempty(randPt)
    error('Cannot obtain a point.')
end
