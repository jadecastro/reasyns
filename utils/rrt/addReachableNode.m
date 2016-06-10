function qNew = addReachableNode(sys,qNear,U0,stepSize,options)

x0 = computeOpenLoopTrajectory(sys,qNear,U0,stepSize,options);

t = x0.getTimeVec;
qNew = double(x0,t(end))';

%[Xk,t] = downsampleUniformly(x0,options.sampSkipColl);
% t = t';
%Xk = Xk';
% Uk = double(u0,t)';
%qNew = Xk(end,:);

% if length(t) > 1
% 
%     clear xknrm
%     xknrm(1) = 0;
%     for k = 1:length(t)-1, xknrm(k+1) = norm(Xk(k+1,1:2) - Xk(k,1:2)); end
%     indxRem = (cumsum(xknrm) > stepSize);
%     t(indxRem) = [];
%     Xk(indxRem,:) = [];
%     Uk(indxRem,:) = [];
% 
% 
% %     plot(Xk(:,1),Xk(:,2),'g:')
% %     drawnow
% 
%     qNew = Xk(end,:);
% 
% end