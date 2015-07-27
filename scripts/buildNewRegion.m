
overApproxEll = ellunion_ea(projection(ac_tmp(j),sys));

[offset,Q] = double(overApproxEll);
[eigvec,eigval] = eig(inv(Q));

const = diag(eigval).^-0.5;
angles = atan2(eigvec(2,:),eigvec(1,:));

newVert(:,1) = offset + [cos(angles(2)) -sin(angles(2));sin(angles(2)) cos(angles(2))]*[const(2)/sqrt(2); const(1)/sqrt(2)];
newVert(:,2) = offset + [cos(angles(2)) -sin(angles(2));sin(angles(2)) cos(angles(2))]*[const(2)/sqrt(2); -const(1)/sqrt(2)];
newVert(:,3) = offset + [cos(angles(2)) -sin(angles(2));sin(angles(2)) cos(angles(2))]*[-const(2)/sqrt(2); -const(1)/sqrt(2)];
newVert(:,4) = offset + [cos(angles(2)) -sin(angles(2));sin(angles(2)) cos(angles(2))]*[-const(2)/sqrt(2); const(1)/sqrt(2)];

newVertT = newVert';
newRegTmp = region(newVertT);
epsilon = 0;
while ~isinside(ac_tmp(j),newRegTmp,sys)
    epsilon = epsilon+0.02;
    newRegTmp = region(newVertT,epsilon);
end

newReg = intersect(newRegTmp,reg(aut.q{iModeToPatch}));

vertsToWrite = newReg.v;

% transform to map coordinates
for idx = 1:size(vertsToWrite,1)
    tmpVert(idx,:) = 1/pix2m*(inv(calibMatrix)*[vertsToWrite(idx,1:2)'; 1])';
end
vertsToWrite = tmpVert(:,1:2);