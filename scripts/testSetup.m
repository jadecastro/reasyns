global t x

% Controller parameters
x_look = 1;
e = 0.2;
gotopt = 1;         % Initialize waypoint index
closeEnough = 0.05;

wlim = 3;

n = 4;      % Dimension of state vector

t = msspoly('t');
x = msspoly('x',n); 

dX = 1.0;

if n > 2
    triangleRegion
elseif n == 2
    triangleRegion1
end

regAutom = [
    1;
    2];
edgeAutom = [
    1 2;
    1 1];  % <-- need this for "stay in r1"
%    2 2];  % <-- need this for "stay in r2"

iEdgeAut = 1;

if edgeAutom(iEdgeAut,1) == edgeAutom(iEdgeAut,2)
    ctrlType = 'stay';  % We're staying in the current state
else
    ctrlType = 'leave';  % We're leaving the current state
end

% Specify regions in this transition
Xinit = Xreg{edgeAutom(iEdgeAut,1)};
Xgoal = Xreg{edgeAutom(iEdgeAut,2)};

% remove intersecting boundaries
Xtmp = Xgoal;
for i = 1:length(Xinit.int)
    for j = 1:length(Xtmp.ext)
        if double(Xinit.int(i) + Xtmp.ext(j)) == 0;
            Xtmp.ext(j) = [];
            Xtmp.extVert(:,j) = [];
            break
        end
    end
end

if strmatch(ctrlType,'stay')
    X = Xinit;
elseif strmatch(ctrlType,'leave')
    % Define the composite region (init+goal)
    X.ext = [Xinit.ext];
    X.extVert = Xinit.extVert;
    X.int = [];
    X.intVert = Xinit.intVert;
    for i = 1:length(Xtmp.ext)
        X.int{i} = [Xinit.int Xtmp.ext(i)];
        X.intVert{i} = [Xinit.intVert Xtmp.extVert(:,i)];
    end
end

reg = X;
clear X
