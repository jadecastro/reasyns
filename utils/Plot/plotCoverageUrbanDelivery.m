

% figure(3)
% clf
% hold on
% axis equal
% 
options.color = [0.9 0.9 0.9];
% plot(pReg{2},options)
% plot(pReg{1},options)
% plot(pReg{3},options)
figure(5)
% clf
hold on
axis equal
% plot(pReg{2},options)
% plot(pReg{1},options)
% plot(pReg{3},options)

kiter = 1;

% TODO: pick out points in ellFunnel to minimize obvious violations?

% 
options.fill = 0.5;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 1;
for i = 1:size(funnelIn,1)
    if ~isempty(funnelIn{i,j,kiter}),
        for k = 1:length(funnelIn{i,j,kiter}.t)
            tmp = inv(funnelIn{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnelIn{i,j,kiter}.x(k,:)',tmp*funnelIn{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 1;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 4;
for i = 1:size(funnelIn,1)
    if ~isempty(funnelIn{i,j,kiter}),
        for k = 1:length(funnelIn{i,j,kiter}.t)
            tmp = inv(funnelIn{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnelIn{i,j,kiter}.x(k,:)',tmp*funnelIn{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 1;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 3;
for i = 1:size(funnelIn,1)
    if ~isempty(funnelIn{i,j,kiter}),
        for k = 1:length(funnelIn{i,j,kiter}.t)
            tmp = inv(funnelIn{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnelIn{i,j,kiter}.x(k,:)',tmp*funnelIn{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 1;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 2;
for i = 1:size(funnelIn,1)
    if ~isempty(funnelIn{i,j,kiter}),
        for k = 1:length(funnelIn{i,j,kiter}.t)
            tmp = inv(funnelIn{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnelIn{i,j,kiter}.x(k,:)',tmp*funnelIn{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end


kiter = 2;

options.fill = 0.5;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 1;
for i = 1:size(funnel,1)
    if ~isempty(funnel{i,j,kiter}),
        for k = 1:length(funnel{i,j,kiter}.t)
            tmp = inv(funnel{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,kiter}.x(k,:)',tmp*funnel{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 0.5;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 2;
for i = 1:size(funnel,1)
    if ~isempty(funnel{i,j,kiter}),
        for k = 1:length(funnel{i,j,kiter}.t)
            tmp = inv(funnel{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,kiter}.x(k,:)',tmp*funnel{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 0.5;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 3;
for i = 1:size(funnel,1)
    if ~isempty(funnel{i,j,kiter}),
        for k = 1:length(funnel{i,j,kiter}.t)
            tmp = inv(funnel{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,kiter}.x(k,:)',tmp*funnel{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 0.5;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 4;
for i = 1:size(funnel,1)
    if ~isempty(funnel{i,j,kiter}),
        for k = 1:length(funnel{i,j,kiter}.t)
            tmp = inv(funnel{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,kiter}.x(k,:)',tmp*funnel{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 0.5;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 5;
for i = 1:size(funnel,1)
    if ~isempty(funnel{i,j,kiter}),
        for k = 1:length(funnel{i,j,kiter}.t)
            tmp = inv(funnel{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,kiter}.x(k,:)',tmp*funnel{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 0.5;
options.color = [0.5 0.5 0.5];
options.shade = 0.5;
j = 6;
for i = 1:size(funnel,1)
    if ~isempty(funnel{i,j,kiter}),
        for k = 1:length(funnel{i,j,kiter}.t)
            tmp = inv(funnel{i,j,kiter}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,kiter}.x(k,:)',tmp*funnel{i,j,kiter}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
% 
% options.fill = 1;
% options.color = [0.3 0.9 0.3];
% j = 2;
% for i = 1:size(ellFunnel,1)
%     if ~isempty(ellFunnel{i,j,kiter}), plotEllipse(projection(ellFunnel{i,j,kiter},[1 0;0 1;0 0]),options.color); end
% end

% for j = 1:size(ellFunnelC,2)
%     for i = 1:size(ellFunnelC,1)
%         if ~isempty(ellFunnelC{i,j,1}), plot(trajFunnelC{i,j,1}(:,1),trajFunnelC{i,j,1}(:,2),'k','LineWidth',2); end
%     end
% end
%%

j = 4;
for i = 1:size(funnel,1)
    if ~isempty(funnel{i,j,kiter}), plot(trajFunnel{i,j,kiter}(:,1),trajFunnel{i,j,kiter}(:,2),'Color',[0.4 0.4 0.4],'LineWidth',2); end
end
j = 2;
for i = 1:size(funnel,1)
    if ~isempty(funnel{i,j,kiter}), plot(trajFunnel{i,j,kiter}(:,1),trajFunnel{i,j,kiter}(:,2),'Color',[0.4 0.4 0.4],'LineWidth',2); end
end
% 

% j = 2;
% for i = 1:size(ellFunnel,1)
%     if ~isempty(ellFunnel{i,j,kiter}), plot(trajFunnel{i,j,kiter}(:,1),trajFunnel{i,j,kiter}(:,2),'k','LineWidth',2); end
% end

%%
vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});
plot([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],'k--','LineWidth',2)
plot([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],'k--','LineWidth',2)

plot([vReg{1}(2,1) vReg{1}(3,1)],[vReg{1}(2,2) vReg{1}(3,2)],'r--','LineWidth',2)
plot([vReg{1}(3,1) vReg{1}(4,1)],[vReg{1}(3,2) vReg{1}(4,2)],'r--','LineWidth',2)
plot([vRegB{1}(2,1) vRegB{1}(3,1)],[vRegB{1}(2,2) vRegB{1}(3,2)],'b--','LineWidth',1)
plot([vRegB{1}(3,1) vRegB{1}(4,1)],[vRegB{1}(3,2) vRegB{1}(4,2)],'b--','LineWidth',1)
plot([vRegBN{1}(2,1) vRegBN{1}(3,1)],[vRegBN{1}(2,2) vRegBN{1}(3,2)],'b--','LineWidth',1)
plot([vRegBN{1}(3,1) vRegBN{1}(4,1)],[vRegBN{1}(3,2) vRegBN{1}(4,2)],'b--','LineWidth',1)

plot([vReg{3}(1,1) vReg{3}(2,1)],[vReg{3}(1,2) vReg{3}(2,2)],'r--','LineWidth',2)
plot([vReg{3}(4,1) vReg{3}(1,1)],[vReg{3}(4,2) vReg{3}(1,2)],'r--','LineWidth',2)
plot([vRegB{3}(1,1) vRegB{3}(2,1)],[vRegB{3}(1,2) vRegB{3}(2,2)],'b--','LineWidth',1)
plot([vRegB{3}(4,1) vRegB{3}(1,1)],[vRegB{3}(4,2) vRegB{3}(1,2)],'b--','LineWidth',1)
plot([vRegBN{3}(1,1) vRegBN{3}(2,1)],[vRegBN{3}(1,2) vRegBN{3}(2,2)],'b--','LineWidth',1)
plot([vRegBN{3}(4,1) vRegBN{3}(1,1)],[vRegBN{3}(4,2) vRegBN{3}(1,2)],'b--','LineWidth',1)

xlabel('x_r'); ylabel('y_r')

set(gca, 'Box', 'off' )
set(gca, 'TickDir', 'out')

%%
figure(21)
clf
hold on
axis equal

options.color = [0.9 0.9 0.9];

plot(pReg{2},options)
plot(pReg{1},options)
plot(pReg{3},options)

% options.fill = 1;
% options.color = [0.9 0.3 0.3];
% for j = 1:size(ellFunnel,2)
%     for i = 1:size(ellFunnel,1)
%         if ~isempty(ellFunnel{i,j,1}), plotEllipse(projection(ellFunnel{i,j,1},[1 0;0 1;0 0]),options.color); end
%     end
% end
% for j = 1:size(ellFunnel,2)
%     for i = 1:size(ellFunnel,1)
%         if ~isempty(ellFunnel{i,j,1}), plot(trajFunnel{i,j,1}(:,1),trajFunnel{i,j,1}(:,2),'k','LineWidth',2); end
%     end
% end

options.fill = 1;
options.color = [0.3 0.6 0.3];
plotEllipse(projection(ellTmp5,[1 0;0 1;0 0]),options.color);
options.color = [0.3 0.9 0.3];
% plotEllipse(projection(ellTmp1,[1 0;0 1;0 0]),options.color);
% options.color = [0.9 0.3 0.3];
% plotEllipse(projection(ellTmp2,[1 0;0 1;0 0]),options.color);
options.color = [0.9 0.3 0.3];
plotEllipse(projection(ellTmp4,[1 0;0 1;0 0]),options.color);
options.color = [0.3 0.9 0.3];
plotEllipse(projection(ellTmp3,[1 0;0 1;0 0]),options.color);

plot(Xk1(1:19,1),Xk1(1:19,2),'k','LineWidth',2)
plot(Xk2(1:15,1),Xk2(1:15,2),'k','LineWidth',2)
plot(Xk3(1:10,1),Xk3(1:10,2),'k','LineWidth',2)
plot(Xk4(:,1),Xk4(:,2),'k','LineWidth',2)
plot(Xk5(:,1),Xk5(:,2),'k','LineWidth',2)

plot(Xk2(15,1),Xk2(15,2),'b+','LineWidth',3,'MarkerSize',14)
plot(Xk3(10,1),Xk3(10,2),'bx','LineWidth',3,'MarkerSize',14)

vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});
plot([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],'k--','LineWidth',2)
plot([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],'k--','LineWidth',2)

plot([vReg{1}(1,1) vReg{1}(2,1)],[vReg{1}(1,2) vReg{1}(2,2)],'k--','LineWidth',2)
plot([vReg{1}(2,1) vReg{1}(3,1)],[vReg{1}(2,2) vReg{1}(3,2)],'r--','LineWidth',2)
plot([vReg{1}(3,1) vReg{1}(4,1)],[vReg{1}(3,2) vReg{1}(4,2)],'r--','LineWidth',2)
plot([vReg{1}(4,1) vReg{1}(1,1)],[vReg{1}(4,2) vReg{1}(1,2)],'k--','LineWidth',2)

plot([vReg{3}(1,1) vReg{3}(2,1)],[vReg{3}(1,2) vReg{3}(2,2)],'r--','LineWidth',2)
plot([vReg{3}(2,1) vReg{3}(3,1)],[vReg{3}(2,2) vReg{3}(3,2)],'k--','LineWidth',2)
plot([vReg{3}(3,1) vReg{3}(4,1)],[vReg{3}(3,2) vReg{3}(4,2)],'k--','LineWidth',2)
plot([vReg{3}(4,1) vReg{3}(1,1)],[vReg{3}(4,2) vReg{3}(1,2)],'r--','LineWidth',2)

xlabel('x_r'); ylabel('y_r')


