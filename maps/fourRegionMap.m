
% Set up map
pix2m = 10/(100*2.65);

% Use a test map for now
bndPos = [145, 73]*pix2m;
bndPoints = [0, 0; 261, 0; 261, 250; 0, 250]*pix2m;
r1pos = [0, 0]*pix2m;
r1points = [0, 0; 80, 0; 80, 80; 0, 80]*pix2m;
r2pos = [45, 135]*pix2m;
r2points = [0, 0; 80, 0; 80, 80; 0, 80]*pix2m;
r3pos = [135, 80]*pix2m;
r3points = [0, 0; 80, 0; 80, 80; 0, 80]*pix2m;
r4pos = [334 73]*pix2m;
r4size = [73 249]*pix2m;

vReg{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:); 
    r1pos+r1points(3,:); 
    r1pos+r1points(4,:)];
vReg{3} = [
    r2pos+r2points(1,:); 
    r2pos+r2points(2,:); 
    r2pos+r2points(3,:); 
    r2pos+r2points(4,:)];
vReg{4} = [
    r3pos+r3points(1,:); 
    r3pos+r3points(2,:); 
    r3pos+r3points(3,:); 
    r3pos+r3points(4,:)];
vReg{2} = [
    [vReg{4}(2,1) vReg{1}(2,2)]; 
    vReg{4}(2,:); 
    vReg{4}(1,:); 
    vReg{4}(4,:); 
    vReg{4}(3,:); 
    [vReg{4}(3,1) vReg{3}(3,2)]; 
    vReg{3}(3,:); 
    vReg{3}(2,:); 
    vReg{3}(1,:); 
    vReg{3}(4,:); 
    [vReg{1}(1,1) vReg{3}(3,2)]; 
    vReg{1}(4,:); 
    vReg{1}(3,:); 
    vReg{1}(2,:)];
% vBnd{1} = [
%     vReg{1}(1,:); 
%     vReg{2}(3,:); 
%     vReg{3}(1,:); 
%     vReg{3}(4,:); 
%     vReg{3}(3,:); 
%     vReg{2}(1,:); 
%     vReg{1}(3,:); 
%     vReg{1}(2,:)];  % Bound is based on all the regions
vBndOuter{1} = [
    vReg{1}(1,:); 
    vReg{2}(1,:); 
    vReg{2}(6,:); 
    vReg{2}(11,:)];
vBnd = vBndOuter;

% Bloated regions
dBloat = 0.2;
% positive bloat: smaller avoid regions
vRegB{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)+dBloat*[-1 0]; 
    r1pos+r1points(3,:)+dBloat*[-1 -1]; 
    r1pos+r1points(4,:)+dBloat*[0 -1]];
vRegB{3} = [
    r2pos+r2points(1,:)+dBloat*[1 1]; 
    r2pos+r2points(2,:)+dBloat*[-1 1]; 
    r2pos+r2points(3,:)+dBloat*[-1 0]; 
    r2pos+r2points(4,:)+dBloat*[1 0]];
vRegB{4} = [
    r3pos+r3points(1,:)+dBloat*[1 1]; 
    r3pos+r3points(2,:)+dBloat*[0 1]; 
    r3pos+r3points(3,:)+dBloat*[0 -1]; 
    r3pos+r3points(4,:)+dBloat*[1 -1]];
vRegB{2} = [
    [vReg{4}(2,1) vReg{1}(2,2)]; 
    vReg{4}(2,:)-dBloat*[0 1]; 
    vReg{4}(1,:)-dBloat*[1 1]; 
    vReg{4}(4,:)-dBloat*[1 -1]; 
    vReg{4}(3,:)-dBloat*[0 -1]; 
    [vReg{4}(3,1) vReg{3}(3,2)]; 
    vReg{3}(3,:)-dBloat*[-1 0]; 
    vReg{3}(2,:)-dBloat*[-1 1]; 
    vReg{3}(1,:)-dBloat*[1 1]; 
    vReg{3}(4,:)-dBloat*[1 0]; 
    [vReg{1}(1,1) vReg{3}(3,2)]; 
    vReg{1}(4,:)-dBloat*[0 -1]; 
    vReg{1}(3,:)-dBloat*[-1 -1]; 
    vReg{1}(2,:)-dBloat*[-1 0]];
% negative bloat: larger interior regions
vRegBN{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)-dBloat*[-1 0]; 
    r1pos+r1points(3,:)-dBloat*[-1 -1]; 
    r1pos+r1points(4,:)-dBloat*[0 -1]];
vRegBN{3} = [
    r2pos+r2points(1,:)-dBloat*[1 1]; 
    r2pos+r2points(2,:)-dBloat*[-1 1]; 
    r2pos+r2points(3,:)-dBloat*[-1 0]; 
    r2pos+r2points(4,:)-dBloat*[1 0]];
vRegBN{4} = [
    r3pos+r3points(1,:)-dBloat*[1 1]; 
    r3pos+r3points(2,:)-dBloat*[0 1]; 
    r3pos+r3points(3,:)-dBloat*[0 -1]; 
    r3pos+r3points(4,:)-dBloat*[1 -1]];

% These regions are defined for checking if funnels are misbehaving... need
% some conservatism
vBndOuterB{1} = [
    vReg{1}(1,:)+0.1*[-1 -1]; 
    vReg{2}(1,:)+0.1*[1 -1]; 
    vReg{2}(6,:)+0.1*[1 1]; 
    vReg{2}(11,:)+0.1*[-1 1]];
vRegB2{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)+(dBloat+0.1)*[-1 0]; 
    r1pos+r1points(3,:)+(dBloat+0.1)*[-1 -1]; 
    r1pos+r1points(4,:)+(dBloat+0.1)*[0 -1]];
vRegB2{3} = [
    r2pos+r2points(1,:)+(dBloat+0.1)*[1 1]; 
    r2pos+r2points(2,:)+(dBloat+0.1)*[-1 1]; 
    r2pos+r2points(3,:)+(dBloat+0.1)*[-1 0]; 
    r2pos+r2points(4,:)+(dBloat+0.1)*[1 0]];
vRegB2{4} = [
    r3pos+r3points(1,:)+(dBloat+0.1)*[1 1]; 
    r3pos+r3points(2,:)+(dBloat+0.1)*[0 1]; 
    r3pos+r3points(3,:)+(dBloat+0.1)*[0 -1]; 
    r3pos+r3points(4,:)+(dBloat+0.1)*[1 -1]];
vRegB2{2} = [
    [vReg{3}(2,1) vReg{1}(2,2)]; 
    vReg{3}(2,:)-(dBloat+0.1)*[0 1]; 
    vReg{3}(1,:)-(dBloat+0.1)*[1 1]; 
    vReg{3}(4,:)-(dBloat+0.1)*[1 0];
    [vReg{1}(1,1) vReg{3}(3,2)]; 
    vReg{1}(4,:)+(dBloat+0.1)*[0 1]; 
    vReg{1}(3,:)+(dBloat+0.1)*[1 1]; 
    vReg{1}(2,:)+(dBloat+0.1)*[1 0]];
vRegBN2{1} = [
    r1pos+r1points(1,:); 
    r1pos+r1points(2,:)-(dBloat+0.1)*[-1 0]; 
    r1pos+r1points(3,:)-(dBloat+0.1)*[-1 -1]; 
    r1pos+r1points(4,:)-(dBloat+0.1)*[0 -1]];
vRegBN2{3} = [
    r2pos+r2points(1,:)-(dBloat+0.1)*[1 1]; 
    r2pos+r2points(2,:)-(dBloat+0.1)*[-1 1]; 
    r2pos+r2points(3,:)-(dBloat+0.1)*[-1 0]; 
    r2pos+r2points(4,:)-(dBloat+0.1)*[1 0]];
vRegBN2{4} = [
    r3pos+r3points(1,:)-(dBloat+0.1)*[1 1]; 
    r3pos+r3points(2,:)-(dBloat+0.1)*[0 1]; 
    r3pos+r3points(3,:)-(dBloat+0.1)*[0 -1]; 
    r3pos+r3points(4,:)-(dBloat+0.1)*[1 -1]];

% Transitions
trans{1} = [1 2];
trans{2} = [2 3];
trans{3} = [3 2];
trans{4} = [2 4];
trans{5} = [4 2];
trans{6} = [2 1];
% trans{7} = [1 1];

% Automaton definition
clear aut
aut.q{1} = 1;
aut.q{2} = 2;
aut.q{3} = 3;
aut.q{4} = 2;
aut.q{5} = 4;
aut.q{6} = 2;
aut.trans{1} = [1 2];
aut.trans{2} = [2 3];
aut.trans{3} = [2 1];
aut.trans{4} = [3 4];
aut.trans{5} = [4 1];
aut.trans{6} = [4 5];
aut.trans{7} = [5 6];
aut.trans{8} = [6 1];
%aut.trans{9} = [1 1];

% Construct polytopes
pReg{1} = polytope(vReg{1});
pReg{2} = polytope(vReg{2});
pReg{3} = polytope(vReg{3});
pReg{4} = polytope(vReg{4});
%pReg{5} = polytope(vReg{5});

pRegB{1} = polytope(vRegB{1});
pRegB2{1} = polytope(vRegB2{1});
pRegBN{1} = polytope(vRegBN{1});
pRegB{2} = polytope(vRegB{2});
pRegB2{2} = polytope(vRegB2{2});
pRegB{3} = polytope(vRegB{3});
pRegB2{3} = polytope(vRegB2{3});
pRegBN{3} = polytope(vRegBN{3});
pRegB{4} = polytope(vRegB{4});
pRegB2{4} = polytope(vRegB2{4});
pRegBN{4} = polytope(vRegBN{4});
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
reg{2}.mssInt{1} = mssReg{1};
reg{2}.mssInt{2} = mssReg{3};
reg{2}.mssInt{3} = mssReg{4};
reg{2}.mssIntB{1} = mssRegB{1};
reg{2}.mssIntB{2} = mssRegB{3};
reg{2}.mssIntB{3} = mssRegB{4};

for i = 1:length(pBnd)
    [H,K] = double(pBnd{i});
    mssBnd{i} = (H*x(1:2)-K)' + eps*sum(x);
    
    regBnd{i}.p = pBnd{i};
    regBnd{i}.hB = hBndB{i};
    regBnd{i}.mss = -mssBnd{i};
    regBnd{i}.v = vBndOuter{i};
end

% Construct composite region for each transition
regX{1}.mssExt = regBnd{1}.mss;
regX{1}.mssExtB = regBnd{1}.mss;
regX{1}.mssInt{1} = mssReg{3};
regX{1}.mssIntB{1} = mssRegB{3};
regX{1}.mssInt{2} = mssReg{4};
regX{1}.mssIntB{2} = mssRegB{4};
regX{1}.pExt = regBnd{1}.p;
regX{1}.hExtB = regBnd{1}.hB;
regX{1}.pInt(1) = pReg{3};
regX{1}.pIntB(1) = pRegB{3};
regX{1}.pIntB2(1) = pRegB2{3};
regX{1}.pInt(2) = pReg{4};
regX{1}.pIntB(2) = pRegB{4};
regX{1}.pIntB2(2) = pRegB2{4};

regX{2} = regX{1};
regX{2}.mssExt = regBnd{1}.mss;
regX{2}.mssExtB = regBnd{1}.mss;
regX{2}.mssInt{1} = mssReg{1};
regX{2}.mssIntB{1} = mssRegB{1};
regX{2}.mssInt{2} = mssReg{4};
regX{2}.mssIntB{2} = mssRegB{4};
regX{2}.pInt(1) = pReg{1};
regX{2}.pIntB(1) = pRegB{1};
regX{2}.pIntB2(1) = pRegB2{1};
regX{2}.pInt(2) = pReg{4};
regX{2}.pIntB(2) = pRegB{4};
regX{2}.pIntB2(2) = pRegB2{4};

regX{3} = regX{1};

regX{4} = regX{2};

regX{5} = regX{1};

regX{6} = regX{1};
regX{6}.mssExt = regBnd{1}.mss;
regX{6}.mssExtB = regBnd{1}.mss;
regX{6}.mssInt{1} = mssReg{1};
regX{6}.mssIntB{1} = mssRegB{1};
regX{6}.mssInt{2} = mssReg{3};
regX{6}.mssIntB{2} = mssRegB{3};
regX{6}.pInt(1) = pReg{1};
regX{6}.pIntB(1) = pRegB{1};
regX{6}.pIntB2(1) = pRegB2{1};
regX{6}.pInt(2) = pReg{3};
regX{6}.pIntB(2) = pRegB{3};
regX{6}.pIntB2(2) = pRegB2{3};

regX{7} = regX{6};

regX{8} = regX{1};

regX{9} = regX{1};
regX{9}.mssExt = -mssReg{1};
regX{9}.mssExtB = -mssRegBN{1};
regX{9}.mssInt{1} = [];
regX{9}.mssIntB{1} = [];
regX{9}.pInt(1) = [];
regX{9}.pIntB(1) = [];
regX{9}.pIntB2(1) = [];

% TODO: Automated method for constructing composite regions
% for itrans = 1:Ntrans
%     regInit = reg{trans{itrans}(1)};
%     regGoal = reg{trans{itrans}(2)};
%     % remove intersecting boundaries
%     flg = 0;
%     for k = 1:length(regInit.mssInt)
%         for i = 1:length(regInit.mssInt{k})
%             for j = 1:length(regGoal.mssExt)
%                 if double(regInit.mssInt{k}(i) + regGoal.mssExt(j)) == 0;
%                     regGoal.mssExt(j) = [];
%                     regGoal.mssExtB(j) = [];
%                     flg = 1;
%                     break
%                 end
%             end
%         end
%     end
%     for i = 1:length(regInit.mssExt)
%         for k = 1:length(regGoal.mssInt)
%             for j = 1:length(regGoal.mssInt{k})
%                 if double(regInit.mssExt(i) + regGoal.mssInt{k}(j)) == 0;
%                     regInit.mssExt(j) = [];
%                     regInit.mssExtB(j) = [];
%                     flg = 2;
%                     break
%                 end
%             end
%         end
%     end
%     % Define the composite region (init+goal)
%     if flg == 1
%         regX.mssExt = [regInit.mssExt];
%         regX.mssExtB = [regInit.mssExtB];
%         regX.mssInt = [];
%         regX.mssIntB = [];
%         for i = 1:length(regGoal.mssExt)
%             regX.mssInt{i} = [regInit.mssInt regGoal.mssExt(i)];
%             regX.mssIntB{i} = [regInit.mssIntB regGoal.mssExtB(i)];
%         end
%     elseif flg == 2
%         regX.mssExt = [regInit.mssExt];
%         regX.mssExtB = [regInit.mssExtB];
%         regX.mssInt = [];
%         regX.mssIntB = [];
%         for i = 1:length(regGoal.mssExt)
%             regX.mssInt{i} = [regInit.mssInt regGoal.mssExt(i)];
%             regX.mssIntB{i} = [regInit.mssIntB regGoal.mssExtB(i)];
%         end
%     end
% end


%%
figure(3)
clf
hold on
axis equal
plot(pReg{2},'g')
plot(pReg{1},'c')
plot(pReg{3},'b')
plot(pReg{4},'y')
%plot(pReg{5},'m')
figure(4)
clf
hold on
axis equal
plot(pReg{2},'g')
plot(pReg{1},'c')
plot(pReg{3},'b')
plot(pReg{4},'y')
%plot(pReg{5},'m')
