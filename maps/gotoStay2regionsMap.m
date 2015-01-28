

% Set up map
pix2m = 10/(100*2.65);

% Use a test map for now
bndPos = [145, 73]*pix2m;
bndPoints = [0, 0; 261, 0; 261, 250; 0, 250]*pix2m;
r1pos = [0, 0]*pix2m;
r1points = [0, 0; 90, 0; 90, 134; 0, 134]*pix2m;
r2pos = [90, 0]*pix2m;
r2points = [0, 0; 90, 0; 90, 134; 0, 134]*pix2m;

vReg{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:); 
    r1pos+r1points(3,:); 
    r1pos+r1points(4,:)];
vReg{2} = [
    r2pos+r2points(1,:); 
    r2pos+r2points(2,:); 
    r2pos+r2points(3,:); 
    r2pos+r2points(4,:)];
vBndOuter{1} = [
    vReg{1}(1,:); 
    vReg{2}(2,:); 
    vReg{2}(3,:); 
    vReg{1}(4,:)];
vBnd = vBndOuter;

% Bloated regions
dBloat = 1.2;
% positive bloat: smaller avoid regions
vRegB{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)-dBloat*[1 0]; 
    r1pos+r1points(3,:)-dBloat*[1 0]; 
    r1pos+r1points(4,:)];
vRegB{2} = [
    r2pos+r2points(1,:)-dBloat*[-1 0];
    r2pos+r2points(2,:); 
    r2pos+r2points(3,:); 
    r2pos+r2points(4,:)-dBloat*[-1 1]];
% negative bloat: larger interior regions
vRegBN{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)+dBloat*[1 0]; 
    r1pos+r1points(3,:)+dBloat*[1 0]; 
    r1pos+r1points(4,:)];
vRegBN{2} = [
    r2pos+r2points(1,:)+dBloat*[-1 0];
    r2pos+r2points(2,:); 
    r2pos+r2points(3,:); 
    r2pos+r2points(4,:)+dBloat*[-1 1]];

% These regions are defined for checking if funnels are misbehaving... need
% some conservatism
vBndOuterB{1} = [
    vBndOuter{1}(1,:)+0.1*[-1 -1]; 
    vBndOuter{1}(2,:)+0.1*[1 -1]; 
    vBndOuter{1}(3,:)+0.1*[1 1]; 
    vBndOuter{1}(4,:)+0.1*[-1 1]];
vRegB2{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)-(dBloat+0.1)*[1 0]; 
    r1pos+r1points(3,:)-(dBloat+0.1)*[1 0]; 
    r1pos+r1points(4,:)];
vRegB2{2} = [
    r2pos+r2points(1,:)-(dBloat+0.1)*[-1 0]; 
    r2pos+r2points(2,:); 
    r2pos+r2points(3,:); 
    r2pos+r2points(4,:)-(dBloat+0.1)*[-1 1]];

dBloat = 0;
vRegBN2{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)+(dBloat+0.1)*[1 0]; 
    r1pos+r1points(3,:)+(dBloat+0.1)*[1 0]; 
    r1pos+r1points(4,:)];
vRegBN2{2} = [
    r2pos+r2points(1,:)+(dBloat+0.1)*[-1 0]; 
    r2pos+r2points(2,:); 
    r2pos+r2points(3,:); 
    r2pos+r2points(4,:)+(dBloat+0.1)*[-1 1]];


% Transitions
clear trans
trans{1} = [1 2];
trans{2} = [2 2];

% Automaton definition
clear aut
aut.q{1} = 1;
aut.q{2} = 2;
aut.trans{1} = [1 2];
aut.trans{2} = [2 2];

isReactive(1) = true;
isReactive(2) = false;

% Construct polytopes
pReg{1} = polytope(vReg{1});
pReg{2} = polytope(vReg{2});

pRegB{1} = polytope(vRegB{1});
pRegB2{1} = polytope(vRegB2{1});
pRegBN{1} = polytope(vRegBN{1});
pRegB{2} = polytope(vRegB{2});
pRegB2{2} = polytope(vRegB2{2});
pRegBN{2} = polytope(vRegBN{2});
pRegBN2{1} = polytope(vRegBN2{1});
[v,c] = double(pRegBN2{1});
hRegBN2{1} = hyperplane(v',c');
pRegBN2{2} = polytope(vRegBN2{2});
[v,c] = double(pRegBN2{2});
hRegBN2{2} = hyperplane(v',c');

pBnd{1} = polytope(vBndOuter{1});
[v,c] = double(pBnd{1});
hBnd{1} = hyperplane(v',c');
pBndB{1} = polytope(vBndOuterB{1});
[v,c] = double(pBndB{1});
hBndB{1} = hyperplane(v',c');

x = msspoly('x',2);
for i = 1:length(pReg)
    [H,K] = double(pReg{i});
    mssReg{i} = (H*x(1:2)-K)' + eps*sum(x);
    [H,K] = double(pRegB{i});
    mssRegB{i} = (H*x(1:2)-K)' + eps*sum(x);
    if ~isempty(double(pRegBN{i}))
        [H,K] = double(pRegBN{i});
        mssRegBN{i} = (H*x(1:2)-K)' + eps*sum(x);
    else
        mssRegBN{i} = mssReg{i};
    end
    
    reg{i}.p = pReg{i};
    reg{i}.pB = pRegB{i};
    reg{i}.pB2 = pRegB2{i};
    reg{i}.hBN2 = hRegBN2{i};
    reg{i}.mssExt = -mssReg{i};
    reg{i}.mssExtB = -mssRegBN{i};
    reg{i}.mssInt{1} = [];
    reg{i}.mssIntB{1} = [];    
    reg{i}.v = vReg{i};
    reg{i}.vB = vRegB{i};
    reg{i}.vBN = vRegBN{i};
end

for i = 1:length(pBnd)
    [H,K] = double(pBnd{i});
    mssBnd{i} = (H*x(1:2)-K)' + eps*sum(x);
    
    regBnd{i}.p = pBnd{i};
    regBnd{i}.hB = hBndB{i};
    regBnd{i}.mss = -mssBnd{i};
    regBnd{i}.v = vBnd{i};
end

% Construct composite region for each transition
regX{1}.mssExt = regBnd{1}.mss;
regX{1}.mssExtB = regBnd{1}.mss;
regX{1}.mssInt{1} = [];
regX{1}.mssIntB{1} = [];
regX{1}.pExt = regBnd{1}.p;
regX{1}.hExtB = regBnd{1}.hB;
regX{1}.pInt = [];
regX{1}.pIntB = [];
regX{1}.pIntB2 = [];

regX{2}.mssExt = regBnd{1}.mss;
regX{2}.mssExtB = regBnd{1}.mss;
regX{2}.mssInt{1} = mssReg{1};
regX{2}.mssIntB{1} = mssRegB{1};
regX{2}.pExt = pReg{2};
regX{2}.hExtB = reg{2}.hBN2;
regX{2}.pInt(1) = pReg{1};
regX{2}.pIntB(1) = pRegB{1};
regX{2}.pIntB2(1) = pRegB2{1};
