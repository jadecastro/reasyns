
global funnelI

debugFlg = 1;

indxSet = 5:6;

clear V ell11

for i = 1:length(funnelI.t), 
    V(i) = (x - [funnelI.x(i,1:2)'; funnelI.x(1,3)])'*funnelI.P(:,:,i)/funnelI.rho(i)*(x - [funnelI.x(i,1:2)'; funnelI.x(1,3)]); 
    tmp = inv(funnelI.P(:,:,i));
    tmp = (tmp+tmp')/2;
    ell11{1,1}(i,1) = ellipsoid(funnelI.x(i,:)',tmp*funnelI.rho(i))';
end
funnelI.x = [funnelI.x(:,1:2) ones(size(funnelI.t))*funnelI.x(1,3)];

%%
figure(4)
clf
hold on
Vbnd{1,1} = ones(size(V)) - V;
VinTmp{1} = [];
ellIn{1} = [];
[funneltmp,rhoMin,rho_d] = computeFunnel(funnelI.t(indxSet),funnelI.x(indxSet,:),funnelI.u(indxSet,:),regBnd,regX{itrans},regAvoid1,regAvoid1,Vbnd,Vbnd,VinTmp,VinTmp,ell11,ell11,ellIn,ellIn);
ellArrayCurr = [];

tmp = inv(funneltmp.P(:,:,end));
tmp = (tmp+tmp')/2;
E1 = ellipsoid(funneltmp.x(end,:)',tmp*rhoMin);
E = projection(E1,[1 0 ;0 1 ;0 0]);
ellArrayCurr = [ellArrayCurr; E1];

% plot(projection(ellArrayCurr,[1 0;0 1; 0 0]))
