

% Set up map
pix2m = 10/(100*2.65);

% Use a test map for now
% bndPos = [145, 73]*pix2m;
bndPoints = [0, 0; 400, 0; 400, 300; 0, 300]*pix2m;
%obstacles
r10pos = [50, 0]*pix2m;
r10points = [0, 0; 170, 0; 60, 180]*pix2m;  % (counterclockwise)
r11pos = [110, 300]*pix2m;
r11points = [0, 0; 130, -140; 250, 0]*pix2m;
r12pos = [330, 0]*pix2m;
r12points = [0, 0; 70, 0; 70, 270]*pix2m;

r1pos = [0, 0]*pix2m;
r1points = [0, 0; 50, 0; 110, 180; 0, 260]*pix2m;
r2pos = [0, 0]*pix2m;
r2points = [110, 180; 40, 300; 0, 300; 0, 260]*pix2m;
r3pos = [0, 0]*pix2m;
r3points = [110, 180; 240, 300; 40, 300]*pix2m;
r4pos = [0, 0]*pix2m;
r4points = [110, 180; 240, 250; 240, 300]*pix2m;
r5pos = [0, 0]*pix2m;
r5points = [110, 180; 160, 0; 240, 160; 240, 250]*pix2m;
r6pos = [0, 0]*pix2m;
r6points = [160, 0; 340, 0; 240, 160]*pix2m;
r7pos = [0, 0]*pix2m;
r7points = [340, 0; 400, 200; 240, 250; 240, 160]*pix2m;
r8pos = [0, 0]*pix2m;
r8points = [400, 200; 400, 260; 240, 280; 240, 250]*pix2m;
r9pos = [0, 0]*pix2m;
r9points = [400, 260; 400, 300; 240, 300; 240, 280]*pix2m;


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
    r3pos+r3points(3,:)];
vReg{4} = [
    r4pos+r4points(1,:); 
    r4pos+r4points(2,:); 
    r4pos+r4points(3,:)];
vReg{5} = [
    r5pos+r5points(1,:); 
    r5pos+r5points(2,:); 
    r5pos+r5points(3,:);
    r5pos+r5points(4,:)];
vReg{6} = [
    r6pos+r6points(1,:); 
    r6pos+r6points(2,:); 
    r6pos+r6points(3,:)];
vReg{7} = [
    r7pos+r7points(1,:); 
    r7pos+r7points(2,:); 
    r7pos+r7points(3,:);
    r7pos+r7points(4,:)];
vReg{8} = [
    r8pos+r8points(1,:); 
    r8pos+r8points(2,:); 
    r8pos+r8points(3,:);
    r8pos+r8points(4,:)];
vReg{9} = [
    r9pos+r9points(1,:); 
    r9pos+r9points(2,:); 
    r9pos+r9points(3,:);
    r9pos+r9points(4,:)];
vReg{10} = [
    r10pos+r10points(1,:); 
    r10pos+r10points(2,:); 
    r10pos+r10points(3,:)];
vReg{11} = [
    r11pos+r11points(1,:); 
    r11pos+r11points(2,:); 
    r11pos+r11points(3,:)];
vReg{12} = [
    r12pos+r12points(1,:); 
    r12pos+r12points(2,:); 
    r12pos+r12points(3,:)];

vBndOuter{1} = [
    bndPoints(1,:); 
    bndPoints(2,:); 
    bndPoints(3,:); 
    bndPoints(4,:)];
vBnd = vBndOuter;

reg1{1} = region(vReg{1});
reg1{2} = region(vReg{2});
reg1{3} = region(vReg{3});
reg1{4} = region(vReg{4});
reg1{5} = region(vReg{5});
reg1{6} = region(vReg{6});
reg1{7} = region(vReg{7});
reg1{8} = region(vReg{8});
reg1{9} = region(vReg{9});
reg1{10} = region(vReg{10});
reg1{11} = region(vReg{11});
reg1{12} = region(vReg{12});
regBndOuter{1} = region(vBndOuter{1});
regBnd1{1} = region(vBnd{1});

% Bloated regions
dBloat = 0.2;

% positive bloat: smaller avoid regions
regBndOuterB{1} = region(vBndOuter{1},0.1);

regB{1} = region(vReg{1},dBloat);
regB{2} = region(vReg{2},dBloat);
regB{3} = region(vReg{3},dBloat);
regB{4} = region(vReg{4},dBloat);
regB{5} = region(vReg{5},dBloat);
regB{6} = region(vReg{6},dBloat);
regB{7} = region(vReg{7},dBloat);
regB{8} = region(vReg{8},dBloat);
regB{9} = region(vReg{9},dBloat);
regB{10} = region(vReg{10});
regB{11} = region(vReg{11});
regB{12} = region(vReg{12});

% negative bloat: larger interior regions
regBN{1} = region(vReg{1},-dBloat);
regBN{2} = region(vReg{2},-dBloat);
regBN{3} = region(vReg{3},-dBloat);
regBN{4} = region(vReg{4},-dBloat);
regBN{5} = region(vReg{5},-dBloat);
regBN{6} = region(vReg{6},-dBloat);
regBN{7} = region(vReg{7},-dBloat);
regBN{8} = region(vReg{8},-dBloat);
regBN{9} = region(vReg{9},-dBloat);
regBN{10} = region(vReg{10});
regBN{11} = region(vReg{11});
regBN{12} = region(vReg{12});

% These regions are defined for checking if funnels are misbehaving... need
% some conservatism
regB2{1} = region(vReg{1},dBloat+0.1);
regB2{2} = region(vReg{2},dBloat+0.1);
regB2{3} = region(vReg{3},dBloat+0.1);
regB2{4} = region(vReg{4},dBloat+0.1);
regB2{5} = region(vReg{5},dBloat+0.1);
regB2{6} = region(vReg{6},dBloat+0.1);
regB2{7} = region(vReg{7},dBloat+0.1);
regB2{8} = region(vReg{8},dBloat+0.1);
regB2{9} = region(vReg{9},dBloat+0.1);
regB2{10} = region(vReg{10},0.1);
regB2{11} = region(vReg{11},0.1);
regB2{12} = region(vReg{12},0.1);

regBN2{1} = region(vReg{1},-dBloat-0.1);
regBN2{2} = region(vReg{2},-dBloat-0.1);
regBN2{3} = region(vReg{3},-dBloat-0.1);
regBN2{4} = region(vReg{4},-dBloat-0.1);
regBN2{5} = region(vReg{5},-dBloat-0.1);
regBN2{6} = region(vReg{6},-dBloat-0.1);
regBN2{7} = region(vReg{7},-dBloat-0.1);
regBN2{8} = region(vReg{8},-dBloat-0.1);
regBN2{9} = region(vReg{9},-dBloat-0.1);
regBN2{10} = region(vReg{10},-0.1);
regBN2{11} = region(vReg{11},-0.1);
regBN2{12} = region(vReg{12},-0.1);



% Transitions (in terms of regions)
clear trans
trans{1} = [1 2];
trans{2} = [2 3];
trans{3} = [3 4];
trans{4} = [4 5];
trans{5} = [5 6];
trans{6} = [6 7];
trans{7} = [7 8];
trans{8} = [8 7];
trans{9} = [7 6];
trans{10} = [6 5];
trans{11} = [5 4];
trans{12} = [4 3];
trans{13} = [3 2];

% Automaton definition
clear aut
aut.q{1} = 1;
aut.q{2} = 2;
aut.q{3} = 3;
aut.q{4} = 4;
aut.q{5} = 5;
aut.q{6} = 6;
aut.q{7} = 7;
aut.q{8} = 8;
aut.q{9} = 7;
aut.q{10} = 6;
aut.q{11} = 5;
aut.q{12} = 4;
aut.q{13} = 3;
aut.trans{1} = [1 2];
aut.trans{2} = [2 3];
aut.trans{3} = [3 4];
aut.trans{4} = [4 5];
aut.trans{5} = [5 6];
aut.trans{6} = [6 7];
aut.trans{7} = [7 8];
aut.trans{8} = [8 9];
aut.trans{9} = [9 10];
aut.trans{10} = [10 11];
aut.trans{11} = [11 12];
aut.trans{12} = [12 13];
aut.trans{13} = [13 2];

isReactive(1) = false;
isReactive(2) = false;
isReactive(3) = false;
isReactive(4) = false;
isReactive(5) = false;
isReactive(6) = false;
isReactive(7) = false;
isReactive(8) = false;
isReactive(9) = false;
isReactive(10) = false;
isReactive(11) = false;
isReactive(12) = false;
isReactive(13) = false;

% Construct polytopes
for i = 1:length(vReg)
    pReg{i} = reg1{i}.p_bnd;
    pRegB{i} = regB{i}.p_bnd;
    pRegB2{i} = regB2{i}.p_bnd;
    pRegBN{i} = regBN{i}.p_bnd;
    pRegBN2{i} = regBN2{i}.p_bnd;
    [v,c] = double(pRegBN2{i});
    hRegBN2{i} = hyperplane(v',c');
end

pBnd{1} = polytope(vBndOuter{1});
[v,c] = double(pBnd{1});
hBnd{1} = hyperplane(v',c');
pBndB{1} = regBndOuterB{1}.p_bnd;
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
    reg{i}.vB = vReg{i};
    reg{i}.vBN = vReg{i};
end
% TODO: automate the setting of internal regions??
reg{3}.mssInt{1} = mssReg{2};
reg{3}.mssIntB{1} = mssRegB{2};

for i = 1:length(pBnd)
    [H,K] = double(pBnd{i});
    mssBnd{i} = (H*x(1:2)-K)' + eps*sum(x);
    
    regBnd{i}.p = pBnd{i};
    regBnd{i}.hB = hBndB{i};
    regBnd{i}.mss = -mssBnd{i};
    regBnd{i}.v = vBnd{i};
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
regX{11} = regX{1};
regX{12} = regX{1};
regX{13} = regX{1};

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

% testing...
% indx = length(aut.trans)+1;
% count = 0;
% for i = 4
%     count = count + 1;
%     regX{indx}.mssInt{count} = mssReg{i};
%     regX{indx}.mssIntB{count} = mssRegB{i};
%     regX{indx}.pInt(count) = pReg{i};
%     regX{indx}.pIntB(count) = pRegB{i};
%     regX{indx}.pIntB2(count) = pRegB2{i};
% end
% indx = length(aut.trans)+2;
% count = 0;
% count = count + 1;
% regX{indx}.mssInt{count} = [];
% regX{indx}.mssIntB{count} = [];
% regX{indx}.pInt = [];
% regX{indx}.pIntB = [];
% regX{indx}.pIntB2 = [];

