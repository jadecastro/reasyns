function [randPt] = getCenterRand_new(sys,reg,ac_trans,options,varargin)

H = sys.params.H;
n = sys.params.n;
limsNonRegState = sys.params.limsNonRegState;
isCyclic = sys.params.isCyclic;
Qrand = options.Qrand;

if length(H) ~= 2, error('Workspaces of dimension other than two not yet supported. Sorry.'); end
if n == length(H), error('Empty non-region state vector not yet supported. Sorry.'); end
if size(limsNonRegState,2) ~= n - length(H), error('Dimension mismatch in limsNonRegState.'); end

minNonRegStates = limsNonRegState(1,:);
maxNonRegStates = limsNonRegState(2,:);

if nargin < 5
    tmpNonRegStates = minNonRegStates + (maxNonRegStates - minNonRegStates).*rand(1,n-length(H));
else
    qCenter = varargin{1}';
    tmpNonRegStates = qCenter(length(H)+1:n);
end
    
if nargin < 5
    vReg = reg.v;
    maxReg = max(vReg);
    minReg = min(vReg);
    tmpRegStates = (H\((1/2*(maxReg + minReg) + Qrand*randn(1,length(H)))*H)')';
else
    tmpRegStates = (H\((qCenter(1:length(H))  + Qrand*randn(1,length(H)))*H)')';
end

% If in obstacle or outside invariant, recompute its location in a deterministic fashion
for indx = 2:100
    isect = isinside(reg,sys,tmpRegStates);
    if isect
        if ~isempty(ac_trans)
            res = true;
            for i = 1:length(ac_trans)
                res = res & isinternal(ac_trans(i),[tmpRegStates tmpNonRegStates]','u');
            end
            if res
                break
            else
                tmpNonRegStates = minNonRegStates + (maxNonRegStates - minNonRegStates).*rand(1,n-length(H));
            end
        else
            break
        end
    end
    if nargin < 5
        tmpRegStates = (H\((1/2*(maxReg + minReg) + Qrand*randn(1,length(H)))*H)')';
    else
        tmpRegStates = (H\((qCenter(1:length(H))  + Qrand*randn(1,length(H)))*H)')';
    end
end

if indx == 100
    error('Can''t find a non-colliding point.')
end

randPt = [tmpRegStates tmpNonRegStates];

if isempty(randPt)
    error('Cannot obtain a point.')
end
