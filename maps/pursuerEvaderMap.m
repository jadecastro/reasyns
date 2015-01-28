

% Set up map
pix2m = 10/(100*2.65);

% Use a test map for now
bndPos = [145, 73]*pix2m;
bndPoints = [0, 0; 261, 0; 261, 250; 0, 250]*pix2m;
r1pos = [0, 0]*pix2m;
r1points = [0, 0; 80, 0; 80, 240; 0, 240]*pix2m;  % (counterclockwise)
r2pos = [80, 120]*pix2m;
r2points = [0, 0; 60, 0; 180, 120; 0, 120]*pix2m;
r3pos = [80, 0]*pix2m;
r3points = [0, 0; 180, 0; 60, 120; 0, 120]*pix2m;
r4pos = [225, 60]*pix2m;
r4points = [0, 0; 80, 0; 80, 120; 0, 120]*pix2m;
r5pos = [334 73]*pix2m;
r5size = [73 249]*pix2m;

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
vReg{3} = [
    r3pos+r3points(1,:); 
    r3pos+r3points(2,:); 
    r3pos+r3points(3,:); 
    r3pos+r3points(4,:)];
vReg{5} = [
    r4pos+r4points(1,:); 
    r4pos+r4points(2,:); 
    r4pos+r4points(3,:); 
    r4pos+r4points(4,:)];
vReg{4} = [
    vReg{3}(2,:); 
    [vReg{5}(2,1) vReg{3}(2,2)];
    vReg{5}(2,:); 
    vReg{5}(1,:); 
    vReg{5}(4,:); 
    vReg{5}(3,:); 
    [vReg{5}(3,1) vReg{2}(3,2)]; 
    vReg{2}(3,:); 
    vReg{2}(2,:);
    ];
vBndOuter{1} = [
    vReg{1}(1,:); 
    vReg{4}(2,:); 
    vReg{4}(7,:); 
    vReg{1}(4,:)];
vBnd = vBndOuter;

% Bloated regions
dBloat = 0.2;
% positive bloat: smaller avoid regions
vRegB{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)+dBloat*[-1 0]; 
    r1pos+r1points(3,:)+dBloat*[-1 0]; 
    r1pos+r1points(4,:)];
vRegB{2} = [
    r2pos+r2points(1,:)+dBloat*[1 0]; 
    r2pos+r2points(2,:)+dBloat*[-1.414 0]; 
    r2pos+r2points(3,:)+dBloat*[-1.414 0]; 
    r2pos+r2points(4,:)+dBloat*[1 0]];
vRegB{3} = [
    r3pos+r3points(1,:)+dBloat*[1 0]; 
    r3pos+r3points(2,:)+dBloat*[-1.414 0]; 
    r3pos+r3points(3,:)+dBloat*[-1.414 0]; 
    r3pos+r3points(4,:)+dBloat*[1 0]];
vRegB{5} = [
    r4pos+r4points(1,:)+dBloat*[1 1]; 
    r4pos+r4points(2,:)+dBloat*[0 1]; 
    r4pos+r4points(3,:)+dBloat*[0 -1]; 
    r4pos+r4points(4,:)+dBloat*[1 -1]];
vRegB{4} = [
    vReg{4}(1,:)+dBloat*[1.414 0]; 
    vReg{4}(2,:)+dBloat*[0 0]; 
    vReg{4}(3,:)+dBloat*[0 -1]; 
    vReg{4}(4,:)+dBloat*[-1 -1]; 
    vReg{4}(5,:)+dBloat*[-1 1]; 
    vReg{4}(6,:)+dBloat*[0 1]; 
    vReg{4}(7,:)+dBloat*[0 0]; 
    vReg{4}(8,:)+dBloat*[1.414 0]; 
    vReg{4}(9,:)+dBloat*[1.414 0]; 
    ];
% negative bloat: larger interior regions
vRegBN{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)-dBloat*[-1 0]; 
    r1pos+r1points(3,:)-dBloat*[-1 0]; 
    r1pos+r1points(4,:)];
vRegBN{2} = [
    r2pos+r2points(1,:)-dBloat*[1 0]; 
    r2pos+r2points(2,:)-dBloat*[-1.414 0]; 
    r2pos+r2points(3,:)-dBloat*[-1.414 0]; 
    r2pos+r2points(4,:)-dBloat*[1 0]];
vRegBN{3} = [
    r3pos+r3points(1,:)-dBloat*[1 0]; 
    r3pos+r3points(2,:)-dBloat*[-1.414 0]; 
    r3pos+r3points(3,:)-dBloat*[-1.414 0]; 
    r3pos+r3points(4,:)-dBloat*[1 0]];
vRegBN{5} = [
    r4pos+r4points(1,:)-dBloat*[1 1]; 
    r4pos+r4points(2,:)-dBloat*[0 1]; 
    r4pos+r4points(3,:)-dBloat*[0 -1]; 
    r4pos+r4points(4,:)-dBloat*[1 -1]];
vRegBN{4} = [
    vReg{4}(1,:)-dBloat*[1.414 0]; 
    vReg{4}(2,:)-dBloat*[0 0]; 
    vReg{4}(3,:)-dBloat*[0 -1]; 
    vReg{4}(4,:)-dBloat*[-1 -1]; 
    vReg{4}(5,:)-dBloat*[-1 1]; 
    vReg{4}(6,:)-dBloat*[0 1]; 
    vReg{4}(7,:)-dBloat*[0 0]; 
    vReg{4}(8,:)-dBloat*[1.414 0]; 
    vReg{4}(9,:)-dBloat*[1.414 0]; 
    ];

% These regions are defined for checking if funnels are misbehaving... need
% some conservatism
vBndOuterB{1} = [
    vBndOuter{1}(1,:)+0.1*[-1 -1]; 
    vBndOuter{1}(2,:)+0.1*[1 -1]; 
    vBndOuter{1}(3,:)+0.1*[1 1]; 
    vBndOuter{1}(4,:)+0.1*[-1 1]];
vRegB2{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)+(dBloat+0.1)*[-1 0]; 
    r1pos+r1points(3,:)+(dBloat+0.1)*[-1 0]; 
    r1pos+r1points(4,:)];
vRegB2{2} = [
    r2pos+r2points(1,:)+(dBloat+0.1)*[1 0]+0.1*[0 1]; 
    r2pos+r2points(2,:)+(dBloat+0.1)*[-1 0]+0.1*[-(2*0.707-1) 1]; 
    r2pos+r2points(3,:)+(dBloat+0.1)*[-1.414 0]; 
    r2pos+r2points(4,:)+(dBloat+0.1)*[1 0]];
vRegB2{3} = [
    r3pos+r3points(1,:)+(dBloat+0.1)*[1 0]; 
    r3pos+r3points(2,:)+(dBloat+0.1)*[-1.414 0]; 
    r3pos+r3points(3,:)+(dBloat+0.1)*[-1 0]+0.1*[-(2*0.707-1) -1]; 
    r3pos+r3points(4,:)+(dBloat+0.1)*[1 0]+0.1*[0 -1]];
vRegB2{5} = [
    r4pos+r4points(1,:)+(dBloat+0.1)*[1 1]; 
    r4pos+r4points(2,:)+(dBloat+0.1)*[0 1]; 
    r4pos+r4points(3,:)+(dBloat+0.1)*[0 -1]; 
    r4pos+r4points(4,:)+(dBloat+0.1)*[1 -1]];
vRegB2{4} = [
    vReg{4}(1,:)+(dBloat+0.1)*[1.414 0]; 
    vReg{4}(2,:)+(dBloat+0.1)*[0 0]; 
    vReg{4}(3,:)+(dBloat+0.1)*[0 -1]; 
    vReg{4}(4,:)+(dBloat+0.1)*[-1 -1]; 
    vReg{4}(5,:)+(dBloat+0.1)*[-1 1]; 
    vReg{4}(6,:)+(dBloat+0.1)*[0 1]; 
    vReg{4}(7,:)+(dBloat+0.1)*[0 0]; 
    vReg{4}(8,:)+(dBloat+0.1)*[1.414 0]; 
    vReg{4}(9,:)+(dBloat+0.1)*[1.414 0]; 
    ];
vRegBN2{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)-(dBloat+0.1)*[-1 0]; 
    r1pos+r1points(3,:)-(dBloat+0.1)*[-1 0]; 
    r1pos+r1points(4,:)];
vRegBN2{2} = [
    r2pos+r2points(1,:)-((dBloat+0.1)*[1 0]+0.1*[0 1]); 
    r2pos+r2points(2,:)-((dBloat+0.1)*[-1 0]+0.1*[-(2*0.707-1) 1]); 
    r2pos+r2points(3,:)-(dBloat+0.1)*[-1.414 0]; 
    r2pos+r2points(4,:)-(dBloat+0.1)*[1 0]];
vRegBN2{3} = [
    r3pos+r3points(1,:)-(dBloat+0.1)*[1 0]; 
    r3pos+r3points(2,:)-(dBloat+0.1)*[-1.414 0]; 
    r3pos+r3points(3,:)-((dBloat+0.1)*[-1 0]+0.1*[-(2*0.707-1) -1]); 
    r3pos+r3points(4,:)-((dBloat+0.1)*[1 0]+0.1*[0 -1])];
vRegBN2{5} = [
    r4pos+r4points(1,:)-(dBloat+0.1)*[1 1]; 
    r4pos+r4points(2,:)-(dBloat+0.1)*[0 1]; 
    r4pos+r4points(3,:)-(dBloat+0.1)*[0 -1]; 
    r4pos+r4points(4,:)-(dBloat+0.1)*[1 -1]];
vRegBN2{4} = [
    vReg{4}(1,:)-(dBloat+0.1)*[1.414 0]; 
    vReg{4}(2,:)-(dBloat+0.1)*[0 0]; 
    vReg{4}(3,:)-(dBloat+0.1)*[0 -1]; 
    vReg{4}(4,:)-(dBloat+0.1)*[-1 -1]; 
    vReg{4}(5,:)-(dBloat+0.1)*[-1 1]; 
    vReg{4}(6,:)-(dBloat+0.1)*[0 1]; 
    vReg{4}(7,:)-(dBloat+0.1)*[0 0]; 
    vReg{4}(8,:)-(dBloat+0.1)*[1.414 0]; 
    vReg{4}(9,:)-(dBloat+0.1)*[1.414 0]; 
    ];

% Transitions
trans{1} = [1 2];
trans{2} = [1 3];
trans{3} = [2 4];
trans{4} = [3 4];
trans{5} = [4 5];
trans{6} = [5 6];
trans{7} = [6 7];
trans{8} = [6 8];
trans{9} = [7 1];
trans{10} = [8 1];

% Automaton definition
clear aut
aut.q{1} = 1;
aut.q{2} = 2;
aut.q{3} = 3;
aut.q{4} = 4;
aut.q{5} = 5;
aut.q{6} = 4;
aut.q{7} = 2;
aut.q{8} = 3;
aut.trans{1} = [1 2];
aut.trans{2} = [1 3];
aut.trans{3} = [2 4];
aut.trans{4} = [3 4];
aut.trans{5} = [4 5];
aut.trans{6} = [5 6];
aut.trans{7} = [6 7];
aut.trans{8} = [6 8];
aut.trans{9} = [7 1];
aut.trans{10} = [8 1];

% Construct polytopes
pReg{1} = polytope(vReg{1});
pReg{2} = polytope(vReg{2});
pReg{3} = polytope(vReg{3});
pReg{4} = polytope(vReg{4});
pReg{5} = polytope(vReg{5});

pRegB{1} = polytope(vRegB{1});
pRegB2{1} = polytope(vRegB2{1});
pRegBN{1} = polytope(vRegBN{1});
pRegB{2} = polytope(vRegB{2});
pRegB2{2} = polytope(vRegB2{2});
pRegBN{2} = polytope(vRegBN{2});
pRegB{3} = polytope(vRegB{3});
pRegB2{3} = polytope(vRegB2{3});
pRegBN{3} = polytope(vRegBN{3});
pRegB{4} = polytope(vRegB{4});
pRegB2{4} = polytope(vRegB2{4});
pRegBN{4} = polytope(vRegBN{4});
pRegB{5} = polytope(vRegB{5});
pRegB2{5} = polytope(vRegB2{5});
pRegBN{5} = polytope(vRegBN{5});
pRegBN2{1} = polytope(vRegBN2{1});
[v,c] = double(pRegBN2{1});
hRegBN2{1} = hyperplane(v',c');
pRegBN2{2} = polytope(vRegBN2{2});
[v,c] = double(pRegBN2{2});
hRegBN2{2} = hyperplane(v',c');
pRegBN2{3} = polytope(vRegBN2{3});
[v,c] = double(pRegBN2{3});
hRegBN2{3} = hyperplane(v',c');
pRegBN2{4} = polytope(vRegBN2{4});
[v,c] = double(pRegBN2{4});
hRegBN2{4} = hyperplane(v',c');
pRegBN2{5} = polytope(vRegBN2{5});
[v,c] = double(pRegBN2{5});
hRegBN2{5} = hyperplane(v',c');

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
end
% TODO: automate the setting of internal regions??
reg{4}.mssInt{1} = mssReg{5};
reg{4}.mssIntB{1} = mssRegB{5};

for i = 1:length(pBnd)
    [H,K] = double(pBnd{i});
    mssBnd{i} = (H*x(1:2)-K)' + eps*sum(x);
    
    regBnd{i}.p = pBnd{i};
    regBnd{i}.hB = hBndB{i};
    regBnd{i}.mss = -mssBnd{i};
    regBnd{i}.v = vBndOuter{i};
end

%%
% Construct composite region for each transition
clear regX
regX{1}.mssExt = regBnd{1}.mss;
regX{1}.mssExtB = regBnd{1}.mss;
regX{1}.pExt = regBnd{1}.p;
regX{1}.hExtB = regBnd{1}.hB;
regX{2} = regX{1};
regX{3} = regX{1};
regX{4} = regX{1};
regX{5} = regX{1};
regX{6} = regX{1};
regX{7} = regX{1};
regX{8} = regX{1};
regX{9} = regX{1};
regX{10} = regX{1};

for indx = 1:length(aut.trans)
    count = 0;
    for i = 1:length(mssReg)
        if all(i ~= [aut.q{aut.trans{indx}}])
            count = count + 1;
            regX{indx}.mssInt{count} = mssReg{i};
            regX{indx}.mssIntB{count} = mssRegB{i};
            regX{indx}.pInt(count) = pReg{i};
            regX{indx}.pIntB(count) = pRegB{i};
            regX{indx}.pIntB2(count) = pRegB2{i};
        end
    end
end
