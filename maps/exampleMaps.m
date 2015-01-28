
% Set up map
pix2m = 10/(100*2.65);

% Use a test map for now
bndPos = [145, 73]*pix2m;
bndPoints = [0, 0; 261, 0; 261, 250; 0, 250]*pix2m;
r1pos = [0, 0]*pix2m;%[146, 73]*pix2m;
r1points = [0, 0; 90, 0; 90, 134; 0, 134]*pix2m;
r2pos = [175, 100]*pix2m;
r2points = [0, 0; 90, 0; 90, 134; 0, 134]*pix2m;
r3pos = [334 73]*pix2m;
r3size = [73 249]*pix2m;

%vBnd = [bndPos+bndPoints(1,:); bndPos+bndPoints(2,:); bndPos+bndPoints(3,:); bndPos+bndPoints(4,:)];
vReg{1} = [r1pos+r1points(1,:); r1pos+r1points(2,:); r1pos+r1points(3,:); r1pos+r1points(4,:)];
vReg{3} = [r2pos+r2points(1,:); r2pos+r2points(2,:); r2pos+r2points(3,:); r2pos+r2points(4,:)];
% vReg{2} = [[vReg{3}(2,1) vReg{1}(3,2)]; vReg{3}(2,:); [vReg{1}(1,1) vReg{3}(2,2)]; vReg{1}(4,:)];
% vReg{2} = [[vReg{3}(2,1) vReg{1}(2,2)]; vReg{3}(3,:); [vReg{1}(1,1) vReg{3}(3,2)]; vReg{1}(1,:)];
vReg{2} = [[vReg{3}(2,1) vReg{1}(2,2)]; vReg{3}(2,:); vReg{3}(1,:); vReg{3}(4,:); 
    [vReg{1}(1,1) vReg{3}(3,2)]; vReg{1}(4,:); vReg{1}(3,:); vReg{1}(2,:)];
% vReg{4} = [r3pos; r3pos+[0 r3size(2)]; r3pos+r3size; r3pos+[r3size(1) 0]];
vBnd{1} = [vReg{1}(1,:); vReg{2}(3,:); vReg{3}(1,:); vReg{3}(4,:); vReg{3}(3,:); vReg{2}(1,:); vReg{1}(3,:); vReg{1}(2,:)];  % Bound is based on all the regions
vBndOuter{1} = [vReg{1}(1,:); vReg{2}(1,:); vReg{3}(3,:); vReg{2}(5,:)];

% Bloated regions
dBloat = 0.05;
vRegB{1} = [r1pos+r1points(1,:); r1pos+r1points(2,:)-dBloat*[1 0]; r1pos+r1points(3,:)-dBloat*[1 1]; r1pos+r1points(4,:)-dBloat*[0 1]];
vRegB{3} = [r2pos+r2points(1,:)+dBloat*[1 1]; r2pos+r2points(2,:)+dBloat*[0 1]; r2pos+r2points(3,:); r2pos+r2points(4,:)+dBloat*[1 0]];
vRegB{2} = [[vReg{3}(2,1) vReg{1}(2,2)]; vReg{3}(2,:)-dBloat*[0 1]; vReg{3}(1,:)-dBloat*[1 1]; vReg{3}(4,:)-dBloat*[1 0];
    [vReg{1}(1,1) vReg{3}(3,2)]; vReg{1}(4,:)+dBloat*[0 1]; vReg{1}(3,:)+dBloat*[1 1]; vReg{1}(2,:)+dBloat*[1 0]];

% Transitions
trans{1} = [1 2];
trans{2} = [2 3];
trans{3} = [3 2];
trans{4} = [2 1];

% trans{1} = [1 2];
% trans{2} = [2 1];

% Construct polytopes
pReg{1} = polytope(vReg{1});
pReg{2} = polytope(vReg{2});
pReg{3} = polytope(vReg{3});
% pReg{4} = polytope(vReg{4});
%pReg{5} = polytope(vReg{5});

pRegB{1} = polytope(vRegB{1});
pRegB{2} = polytope(vRegB{2});
pRegB{3} = polytope(vRegB{3});

pBnd{1} = polytope(vBndOuter{1});

x = msspoly('x',n);
for i = 1:length(pReg)
    [H,K] = double(pReg{i});
    mssReg{i} = (H*x(1:2)-K)' + eps*sum(x);
    [H,K] = double(pRegB{i});
    mssRegB{i} = (H*x(1:2)-K)' + eps*sum(x);
    
    reg{i}.p = pReg{i};
    reg{i}.pB = pRegB{i};
    reg{i}.mss = mssReg{i};
    reg{i}.mssB = mssRegB{i};
    reg{i}.v = vReg{i};
    reg{i}.vB = vRegB{i};
end
for i = 1:length(pBnd)
    [H,K] = double(pBnd{i});
    mssBnd{i} = (H*x(1:2)-K)' + eps*sum(x);
    
    regBnd{i}.p = pBnd{i};
    regBnd{i}.mss = -mssBnd{i};
    regBnd{i}.v = vBnd{i};
end

figure(3)
clf
hold on
axis equal
plot(pReg{2},'g')
plot(pReg{1},'c')
plot(pReg{3},'b')
% plot(pReg{4},'y')
%plot(pReg{5},'m')
figure(4)
clf
hold on
axis equal
plot(pReg{2},'g')
plot(pReg{1},'c')
plot(pReg{3},'b')
% plot(pReg{4},'y')
%plot(pReg{5},'m')

