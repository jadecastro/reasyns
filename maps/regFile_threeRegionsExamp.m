

% Set up map
pix2m = 10/(100*2.65);

% Use a test map for now
bndPos = [145, 73]*pix2m;
bndPoints = [0, 0; 261, 0; 261, 250; 0, 250]*pix2m;
r1pos = [0, 0]*pix2m;
r1points = [0, 0; 90, 0; 90, 134; 0, 134]*pix2m;
r2pos = [175, 100]*pix2m;
r2points = [0, 0; 90, 0; 90, 134; 0, 134]*pix2m;
% r3pos = [334 73]*pix2m;
% r3size = [73 249]*pix2m;

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
vReg{2} = [
    [vReg{3}(2,1) vReg{1}(2,2)]; 
    vReg{3}(2,:); vReg{3}(1,:); 
    vReg{3}(4,:); 
    [vReg{1}(1,1) vReg{3}(3,2)]; 
    vReg{1}(4,:); vReg{1}(3,:); 
    vReg{1}(2,:)];
vBndOuter{1} = [
    vReg{1}(1,:); 
    vReg{2}(1,:); 
    vReg{3}(3,:); 
    vReg{2}(5,:)];
vBnd = vBndOuter;

