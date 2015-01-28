function [funnelI,rhoMin,rho_d,info,redindx,Qf] = computeFunnelByDrakeRecursive(t,Uk,Xk,Qf,sampSkipFun)

try
    [funnelI,rhoMin,rho_d,info,redindx] = computeFunnelByDrake(t,Uk,Xk,Qf,sampSkipFun);
%     figure(20)
%     hold on
%     Hout = eye(2);
%     n = 3;
%     for j = 1:1:size(funnelI.x,1)
%         tmp = inv(funnelI.P(:,:,j));
%         tmp = (tmp+tmp')/2;
%         E1 = ellipsoid(funnelI.x(j,:)',tmp*funnelI.rho(j));
%         E = projection(E1,[Hout; zeros(n-length(Hout),length(Hout))]);
%         plot(E)
%     end
catch ME
    if strcmp(ME.identifier, 'Drake:PolynomialTrajectorySystem:InfeasibleRho')
        Noverlap = 0;
        if length(t) > 10
            % bisect the (sub) trajectory
            t1 = t(1:floor(end/2)+Noverlap);
            Uk1 = Uk(1:floor(end/2)+Noverlap,:);
            Xk1 = Xk(1:floor(end/2)+Noverlap,:);
            t2 = t(floor(end/2):end);
            Uk2 = Uk(floor(end/2):end,:);
            Xk2 = Xk(floor(end/2):end,:);
            len1 = length(t1);
            
            %TODO: ensure containment of sequenced funnels - reverse ordering 
            [funnelI2,rhoMin2,rho_d2,info2,redindx2,Qf] = computeFunnelByDrakeRecursive(t2,Uk2,Xk2,Qf,sampSkipFun);
            Qf = funnelI2.P(:,:,1)/funnelI2.rho(1);
%             Qf = getMaximalQ(funnelI2,Xk2(1+Noverlap,:));
            [funnelI1,rhoMin1,rho_d1,info1,redindx1,Qf] = computeFunnelByDrakeRecursive(t1,Uk1,Xk1,Qf,sampSkipFun);
            Qf = funnelI1.P(:,:,1)/funnelI1.rho(1);
%             Qf = getMaximalQ(funnelI1,Xk1(1+Noverlap,:));
            
            %collect
            try
                %TODO: two trajectories/funnels instead of one -- below only good for visualization 
            funnelI.t = [funnelI1.t(1:end-1); funnelI2.t]; % discard last Noverlap points of 1st traj b/c funnel *should* be contained
            funnelI.x = [funnelI1.x(1:end-1,:); funnelI2.x];
            funnelI.u = [funnelI1.u(1:end-1,:); funnelI2.u];
            funnelI.P = cat(3,funnelI1.P(:,:,1:end-1),funnelI2.P);
            funnelI.K = cat(3,funnelI1.K(:,:,1:end-1),funnelI2.K);
            funnelI.rho = [funnelI1.rho(1:end-1); funnelI2.rho];
            funnelI.V = [funnelI1.V(1:end-1) funnelI2.V];
            rhoMin = min(rhoMin1,rhoMin2);
            rho_d = [rho_d1(1:end-1); rho_d2];
            redindx = unique([redindx1 redindx2+len1-1]);
            info = [info1; info2];
            catch
                keyboard
            end
            
            % for testing...
%             figure(20)
%             hold on
%             Hout = eye(2);
%             n = 3;
%             for j = 1:1:size(funnelI2.x,1)
%                 tmp = inv(funnelI2.P(:,:,j));
%                 tmp = (tmp+tmp')/2;
%                 E1 = ellipsoid(funnelI2.x(j,:)',tmp*funnelI2.rho(j));
%                 E = projection(E1,[Hout; zeros(n-length(Hout),length(Hout))]);
%                 plot(E)
%             end
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


