

figure(9)
clf
hold on
axis equal


%%
% TODO: pick out points in ellFunnel to minimize obvious violations?

kindx = 1;

% 
options.fill = 1;
options.color = [0.3 0.3 0.9];
options.shade = 0.5;
j = 5;
clear E
for i = 1:size(ellFunnelC,1)
    if ~isempty(funnelIn{i,j,kindx}),
        for k = 1:length(funnelIn{i,j,kindx}.t)
            tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end

options.fill = 0.5;
options.color = [0.9 0.3 0.3];
options.shade = 0.5;
j = 6;
clear E
for i = 1:size(ellFunnelC,1)
    if ~isempty(funnelIn{i,j,kindx}),
        for k = 1:length(funnelIn{i,j,kindx}.t)
            tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end

options.fill = 1;
options.color = [0.3 0.3 0.9];
options.shade = 0.5;
j = 7;
clear E
for i = 1:size(ellFunnelC,1)
    if ~isempty(funnelIn{i,j,kindx}),
        for k = 1:length(funnelIn{i,j,kindx}.t)
            tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end

% options.fill = 0.5;
% options.color = [0.3 0.3 0.9];
% options.shade = 0.5;
% j = 8;
% clear E
% for i = 1:size(ellFunnelC,1)
%     if ~isempty(funnelIn{i,j,kindx}),
%         for k = 1:length(funnelIn{i,j,kindx}.t)
%             tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
%             tmp = (tmp+tmp')/2;
%             E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
%         end
%         plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
%     end
% end


kindx = 2;

options.fill = 0.5;
options.color = [0.5 0.95 0.5];
options.shade = 0.5;
j = 7;
clear E
for i = 1:size(ellFunnel,1)
    if ~isempty(funnel{i,j,kindx}),
        for k = 1:length(funnel{i,j,kindx}.t)
            tmp = inv(funnel{i,j,kindx}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,kindx}.x(k,:)',tmp*funnel{i,j,kindx}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
for i = 1:size(ellFunnel,1)
    if ~isempty(funnel{i,j,kindx}),
        if ~isempty(ellFunnel{i,j,kindx}), plot(trajFunnel{i,j,kindx}(:,1),trajFunnel{i,j,kindx}(:,2),'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',1); end
    end
end

options.fill = 0.5;
options.color = [0.7 1 0.7];
options.shade = 0.5;
j = 6;
clear E
for i = 1:size(ellFunnel,1)
    if ~isempty(funnel{i,j,kindx}),
        for k = 1:length(funnel{i,j,kindx}.t)
            tmp = inv(funnel{i,j,kindx}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,kindx}.x(k,:)',tmp*funnel{i,j,kindx}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
for i = 1:size(ellFunnel,1)
    if ~isempty(funnel{i,j,kindx}),
        if ~isempty(ellFunnel{i,j,kindx}), plot(trajFunnel{i,j,kindx}(:,1),trajFunnel{i,j,kindx}(:,2),'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',1); end
    end
end

% options.fill = 0.5;
% options.color = [0.5 0.95 0.5];
% options.shade = 0.5;
% j = 1;
% clear E
% for i = 1:size(ellFunnel,1)
%     if ~isempty(funnel{i,j,kindx}),
%         for k = 1:length(funnel{i,j,kindx}.t)
%             tmp = inv(funnel{i,j,kindx}.P(:,:,k));
%             tmp = (tmp+tmp')/2;
%             E(k) = ellipsoid(funnel{i,j,kindx}.x(k,:)',tmp*funnel{i,j,kindx}.rho(k));
%         end
%         plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
%     end
% end


% TODO: optimize color -> B/W conversion!!

%
% options.fill = 1;
% options.color = [0.3 0.9 0.3];
% j = 2;
% for i = 1:size(ellFunnel,1)
%     if ~isempty(ellFunnel{i,j,kindx}), plotEllipse(projection(ellFunnel{i,j,kindx},[1 0;0 1;0 0]),options.color); end
% end

% for j = 1:size(ellFunnelC,2)
%     for i = 1:size(ellFunnelC,1)
%         if ~isempty(ellFunnelC{i,j,1}), plot(trajFunnelC{i,j,1}(:,1),trajFunnelC{i,j,1}(:,2),'k','LineWidth',2); end
%     end
% end
% for i = 1:size(ellFunnel,1)
%     if ~isempty(ellFunnel{i,j,kindx}), plot(trajFunnel{i,j,kindx}(:,1),trajFunnel{i,j,kindx}(:,2),'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',1); end
% end


% j = 2;
% for i = 1:size(ellFunnel,1)
%     if ~isempty(ellFunnel{i,j,kindx}), plot(trajFunnel{i,j,kindx}(:,1),trajFunnel{i,j,kindx}(:,2),'k','LineWidth',2); end
% end

vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});
plot([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],'k--','LineWidth',2)
plot([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],'k--','LineWidth',2)

plot([vReg{1}(2,1) vReg{1}(3,1)],[vReg{1}(2,2) vReg{1}(3,2)],'k--','LineWidth',2)
plot([vReg{1}(3,1) vReg{1}(4,1)],[vReg{1}(3,2) vReg{1}(4,2)],'k--','LineWidth',2)
plot([vRegB{1}(2,1) vRegB{1}(3,1)],[vRegB{1}(2,2) vRegB{1}(3,2)],'b:','LineWidth',2)
plot([vRegB{1}(3,1) vRegB{1}(4,1)],[vRegB{1}(3,2) vRegB{1}(4,2)],'b:','LineWidth',2)
plot([vRegBN{1}(2,1) vRegBN{1}(3,1)],[vRegBN{1}(2,2) vRegBN{1}(3,2)],'b:','LineWidth',2)
plot([vRegBN{1}(3,1) vRegBN{1}(4,1)],[vRegBN{1}(3,2) vRegBN{1}(4,2)],'b:','LineWidth',2)

plot([vReg{2}(2,1) vReg{2}(3,1)],[vReg{2}(2,2) vReg{2}(3,2)],'k--','LineWidth',2)
plot([vReg{2}(3,1) vReg{2}(4,1)],[vReg{2}(3,2) vReg{2}(4,2)],'k--','LineWidth',2)
plot([vRegB{2}(2,1) vRegB{2}(3,1)],[vRegB{2}(2,2) vRegB{2}(3,2)],'b:','LineWidth',2)
plot([vRegB{2}(3,1) vRegB{2}(4,1)],[vRegB{2}(3,2) vRegB{2}(4,2)],'b:','LineWidth',2)
plot([vRegBN{2}(2,1) vRegBN{2}(3,1)],[vRegBN{2}(2,2) vRegBN{2}(3,2)],'b:','LineWidth',2)
plot([vRegBN{2}(3,1) vRegBN{2}(4,1)],[vRegBN{2}(3,2) vRegBN{2}(4,2)],'b:','LineWidth',2)

plot([vReg{3}(2,1) vReg{3}(3,1)],[vReg{3}(2,2) vReg{3}(3,2)],'k--','LineWidth',2)
plot([vReg{3}(3,1) vReg{3}(4,1)],[vReg{3}(3,2) vReg{3}(4,2)],'k--','LineWidth',2)
plot([vRegB{3}(2,1) vRegB{3}(3,1)],[vRegB{3}(2,2) vRegB{3}(3,2)],'b:','LineWidth',2)
plot([vRegB{3}(3,1) vRegB{3}(4,1)],[vRegB{3}(3,2) vRegB{3}(4,2)],'b:','LineWidth',2)
plot([vRegBN{3}(2,1) vRegBN{3}(3,1)],[vRegBN{3}(2,2) vRegBN{3}(3,2)],'b:','LineWidth',2)
plot([vRegBN{3}(3,1) vRegBN{3}(4,1)],[vRegBN{3}(3,2) vRegBN{3}(4,2)],'b:','LineWidth',2)

% plot([vReg{4}(2,1) vReg{4}(3,1)],[vReg{4}(2,2) vReg{4}(3,2)],'k--','LineWidth',2)
% plot([vReg{4}(3,1) vReg{4}(4,1)],[vReg{4}(3,2) vReg{4}(4,2)],'k--','LineWidth',2)
% plot([vRegB{4}(2,1) vRegB{4}(3,1)],[vRegB{4}(2,2) vRegB{4}(3,2)],'b:','LineWidth',2)
% plot([vRegB{4}(3,1) vRegB{4}(4,1)],[vRegB{4}(3,2) vRegB{4}(4,2)],'b:','LineWidth',2)
% plot([vRegBN{4}(2,1) vRegBN{4}(3,1)],[vRegBN{4}(2,2) vRegBN{4}(3,2)],'b:','LineWidth',2)
% plot([vRegBN{4}(3,1) vRegBN{4}(4,1)],[vRegBN{4}(3,2) vRegBN{4}(4,2)],'b:','LineWidth',2)

plot([vReg{5}(3,1) vReg{5}(4,1)],[vReg{5}(3,2) vReg{5}(4,2)],'k--','LineWidth',2)
plot([vReg{5}(4,1) vReg{5}(1,1)],[vReg{5}(4,2) vReg{5}(1,2)],'k--','LineWidth',2)
plot([vReg{5}(1,1) vReg{5}(2,1)],[vReg{5}(1,2) vReg{5}(2,2)],'k--','LineWidth',2)
plot([vRegB{5}(3,1) vRegB{5}(4,1)],[vRegB{5}(3,2) vRegB{5}(4,2)],'b:','LineWidth',2)
plot([vRegB{5}(4,1) vRegB{5}(1,1)],[vRegB{5}(4,2) vRegB{5}(1,2)],'b:','LineWidth',2)
plot([vRegB{5}(1,1) vRegB{5}(2,1)],[vRegB{5}(1,2) vRegB{5}(2,2)],'b:','LineWidth',2)
plot([vRegBN{5}(3,1) vRegBN{5}(4,1)],[vRegBN{5}(3,2) vRegBN{5}(4,2)],'b:','LineWidth',2)
plot([vRegBN{5}(4,1) vRegBN{5}(1,1)],[vRegBN{5}(4,2) vRegBN{5}(1,2)],'b:','LineWidth',2)
plot([vRegBN{5}(1,1) vRegBN{5}(2,1)],[vRegBN{5}(1,2) vRegBN{5}(2,2)],'b:','LineWidth',2)

set(gca,'FontSize',12)
xlabel('x_r'); ylabel('y_r')
text(vReg{1}(1,1)+0.2,vReg{1}(1,2)+0.4,'home','FontSize',16)
text(vReg{2}(2,1)+0.2,vReg{2}(2,2)+0.4,'r_1','FontSize',16)
text(vReg{3}(4,1)+0.2,vReg{3}(4,2)+0.4,'r_2','FontSize',16)
text(vReg{4}(9,1)+0.2,vReg{4}(9,2)+0.4,'r_3','FontSize',16)
text(vReg{5}(1,1)+0.2,vReg{5}(1,2)+0.4,'goal','FontSize',16)

set(gca, 'Box', 'off' )

%%
% figure(21)
% clf
% hold on
% axis equal
% 
% options.color = [0.9 0.9 0.9];
% 
% plot(pReg{2},options)
% plot(pReg{1},options)
% plot(pReg{3},options)
% 
% % options.fill = 1;
% % options.color = [0.9 0.3 0.3];
% % for j = 1:size(ellFunnel,2)
% %     for i = 1:size(ellFunnel,1)
% %         if ~isempty(ellFunnel{i,j,1}), plotEllipse(projection(ellFunnel{i,j,1},[1 0;0 1;0 0]),options.color); end
% %     end
% % end
% % for j = 1:size(ellFunnel,2)
% %     for i = 1:size(ellFunnel,1)
% %         if ~isempty(ellFunnel{i,j,1}), plot(trajFunnel{i,j,1}(:,1),trajFunnel{i,j,1}(:,2),'k','LineWidth',2); end
% %     end
% % end
% 
% options.fill = 1;
% options.color = [0.3 0.6 0.3];
% plotEllipse(projection(ellTmp5,[1 0;0 1;0 0]),options.color);
% options.color = [0.3 0.9 0.3];
% % plotEllipse(projection(ellTmp1,[1 0;0 1;0 0]),options.color);
% % options.color = [0.9 0.3 0.3];
% % plotEllipse(projection(ellTmp2,[1 0;0 1;0 0]),options.color);
% options.color = [0.9 0.3 0.3];
% plotEllipse(projection(ellTmp4,[1 0;0 1;0 0]),options.color);
% options.color = [0.3 0.9 0.3];
% plotEllipse(projection(ellTmp3,[1 0;0 1;0 0]),options.color);
% 
% plot(Xk1(1:19,1),Xk1(1:19,2),'k','LineWidth',2)
% plot(Xk2(1:15,1),Xk2(1:15,2),'k','LineWidth',2)
% plot(Xk3(1:10,1),Xk3(1:10,2),'k','LineWidth',2)
% plot(Xk4(:,1),Xk4(:,2),'k','LineWidth',2)
% plot(Xk5(:,1),Xk5(:,2),'k','LineWidth',2)
% 
% plot(Xk2(15,1),Xk2(15,2),'b+','LineWidth',3,'MarkerSize',14)
% plot(Xk3(10,1),Xk3(10,2),'bx','LineWidth',3,'MarkerSize',14)
% 
% vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});
% plot([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
% plot([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],'k--','LineWidth',2)
% plot([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
% plot([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],'k--','LineWidth',2)
% 
% plot([vReg{1}(1,1) vReg{1}(2,1)],[vReg{1}(1,2) vReg{1}(2,2)],'k--','LineWidth',2)
% plot([vReg{1}(2,1) vReg{1}(3,1)],[vReg{1}(2,2) vReg{1}(3,2)],'r--','LineWidth',2)
% plot([vReg{1}(3,1) vReg{1}(4,1)],[vReg{1}(3,2) vReg{1}(4,2)],'r--','LineWidth',2)
% plot([vReg{1}(4,1) vReg{1}(1,1)],[vReg{1}(4,2) vReg{1}(1,2)],'k--','LineWidth',2)
% 
% plot([vReg{3}(1,1) vReg{3}(2,1)],[vReg{3}(1,2) vReg{3}(2,2)],'r--','LineWidth',2)
% plot([vReg{3}(2,1) vReg{3}(3,1)],[vReg{3}(2,2) vReg{3}(3,2)],'k--','LineWidth',2)
% plot([vReg{3}(3,1) vReg{3}(4,1)],[vReg{3}(3,2) vReg{3}(4,2)],'k--','LineWidth',2)
% plot([vReg{3}(4,1) vReg{3}(1,1)],[vReg{3}(4,2) vReg{3}(1,2)],'r--','LineWidth',2)
% 
% set(gca,'FontSize',14)
% xlabel('x_r'); ylabel('y_r')
% text(vReg{1}(1,1)+0.2,vReg{1}(1,2)+0.4,'r_1','FontSize',16)
% text(vReg{1}(2,1)+0.2,vReg{1}(1,2)+0.4,'r_2','FontSize',16)
% text(vReg{3}(1,1)+0.2,vReg{3}(1,2)+0.4,'r_3','FontSize',16)
% 
% set(gca, 'Box', 'off' )
% % set(gca, 'TickDir', 'out')
% 
