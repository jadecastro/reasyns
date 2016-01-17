
function newVertT = buildNewRegion(ell, isOverApprox)

%overApproxEll = ellunion_ea(projection(ac_tmp(j),sys));

[offset,Q] = double(ell);
[eigvec,eigval] = eig(inv(Q));

const = diag(eigval).^-0.5;
if isOverApprox
    const = const/sqrt(2);
end
angles = atan2(eigvec(2,:),eigvec(1,:));

newVert(:,1) = offset + sqrt(2)*[cos(angles(2)) -sin(angles(2));sin(angles(2)) cos(angles(2))]*[const(2); const(1)];
newVert(:,2) = offset + sqrt(2)*[cos(angles(2)) -sin(angles(2));sin(angles(2)) cos(angles(2))]*[const(2); -const(1)];
newVert(:,3) = offset + sqrt(2)*[cos(angles(2)) -sin(angles(2));sin(angles(2)) cos(angles(2))]*[-const(2); -const(1)];
newVert(:,4) = offset + sqrt(2)*[cos(angles(2)) -sin(angles(2));sin(angles(2)) cos(angles(2))]*[-const(2); const(1)];

newVertT = newVert';
newRegTmp = region(newVertT);
% newRegTmp = region(newVertT,0.06);
newVertT = newRegTmp.v;

% newReg = intersect(newRegTmp,reg(aut.q{iModeToPatch}));
% 
% vertsToWrite = newReg.v;
% 
% % transform to map coordinates
% for idx = 1:size(vertsToWrite,1)
%     tmpVert(idx,:) = 1/pix2m*(inv(calibMatrix)*[vertsToWrite(idx,1:2)'; 1])';
% end
% vertsToWrite = tmpVert(:,1:2);
