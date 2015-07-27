function [ac] = computeAtomicControllerSegmentDubins(u0,x0,sys,ctrloptions,sampSkip,xMssExt)

try
    ac = computeAtomicControllerDubins(u0,x0,sys,ctrloptions,sampSkip,xMssExt);

catch ME
    if strcmp(ME.identifier, 'Drake:PolynomialTrajectorySystem:InfeasibleRho')
        if length(x0) > 10
            % bisect the (sub) trajectory
            [x01, x02] = bisect(x0);
            [u01, u02] = bisect(u0);
            
            %TODO: ensure containment of sequenced funnels - reverse ordering 
            ac2 = computeAtomicControllerSegmentDubins(u02,x02,sys,ctrloptions,sampSkip,xMssExt);
            t0 = ac2.P.getTimeVec();
            ctrloptions.Qf = double(ac2.P,t0(1))*double(ac2.rho,t0(1));
            %             Qf = getMaximalQ(funnelI2,Xk2(1+Noverlap,:));
            ac1 = computeAtomicControllerSegmentDubins(u01,x01,sys,ctrloptions,sampSkip,xMssExt);
            t0 = ac1.P.getTimeVec();
            ctrloptions.Qf = double(ac1.P,t0(1))*double(ac1.rho,t0(1));
            %             Qf = getMaximalQ(funnelI1,Xk1(1+Noverlap,:));
            
            %collect into one atomiccontroller
            ac = merge(ac1,ac2); 
        else
            error('Rho infeasible after segmenting the trajecory to its maximal extent. No smaller segments allowed.')
        end
    else
        rethrow(ME)
    end
end
end

function Q = getMaximalQ(funnel,X)

N = 100;
n = length(X);
H = eye(2);
isCyclic = [zeros(n-1,1); 1];

Qi = 1e3*eye(length(X));
Ei = ellipsoid(X',inv(Qi));
E = Ei;
Ep = Ei;

for j = 1:length(funnel.t)
    ellBndInv{1}(j,1).x = funnel.x(j,:)';
    ellBndInv{1}(j,1).P = funnel.P(:,:,j)/funnel.rho(j); 
end

% starting with a tiny rho, iteratively check containment while
% increating rho
for k = 1:100
    [xbar,invQ] = double(E);
    qTest = ellipsoidrand(xbar,invQ,N);
    clear badIndx isect
    [isect,badIndx] = checkIntersection5([],[],ellBndInv,[],qTest,H,n,isCyclic);
    if isect
        break
    end
    Ep = E;
    E = shape(E,1.1);  % increase by 10% and try again..
end
[xbar,invQ] = double(Ep);
Q = inv(invQ);
% if debugFlg
%     plot(projection(Ep,[H; zeros(n-length(H),length(H))]),'g');
%     drawnow
%     %             keyboard
% end

end


