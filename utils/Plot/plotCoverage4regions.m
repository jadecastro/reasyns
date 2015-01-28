

figure(5)
clf
hold on
axis equal

options.color = [0.9 0.9 0.9];
plot(pReg{1},options)
plot(pReg{3},options)
plot(pReg{4},options)

options.fill = 0.5;
options.color = [0.1 0.7 0.1];
options.shade = 0.5;
j = 6;
clear E
for i = 1:size(ellFunnel,1)
    if ~isempty(funnel{i,j,2}),
        for k = 1:length(funnel{i,j,2}.t)
            tmp = inv(funnel{i,j,2}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,2}.x(k,:)',tmp*funnel{i,j,2}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 0.5;
options.color = [0.7 0.1 0.1];
options.shade = 0.5;
j = 2;
clear E
for i = 1:size(ellFunnel,1)
    if ~isempty(funnel{i,j,2}),
        for k = 1:length(funnel{i,j,2}.t)
            tmp = inv(funnel{i,j,2}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,2}.x(k,:)',tmp*funnel{i,j,2}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end
options.fill = 0.5;
options.color = [0.1 0.1 0.7];
options.shade = 0.5;
j = 4;
clear E
for i = setdiff(1:size(ellFunnel,1),[11;15;20]) %1:size(ellFunnel,1)
    if ~isempty(funnel{i,j,2}),
        for k = 1:length(funnel{i,j,2}.t)
            tmp = inv(funnel{i,j,2}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnel{i,j,2}.x(k,:)',tmp*funnel{i,j,2}.rho(k));
        end
        plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
    end
end

% j = 2;
% for i = 1:size(ellFunnel,1)
%     if ~isempty(ellFunnel{i,j,2}), plot(trajFunnel{i,j,2}(:,1),trajFunnel{i,j,2}(:,2),'Color',[0.4 0.4 0.4],'LineWidth',2); end
% end
% j = 4;
% for i = 1:size(ellFunnel,1)
%     if ~isempty(ellFunnel{i,j,2}), plot(trajFunnel{i,j,2}(:,1),trajFunnel{i,j,2}(:,2),'Color',[0.4 0.4 0.4],'LineWidth',2); end
% end
% j = 6;
% for i = 1:size(ellFunnel,1)
%     if ~isempty(ellFunnel{i,j,2}), plot(trajFunnel{i,j,2}(:,1),trajFunnel{i,j,2}(:,2),'Color',[0.4 0.4 0.4],'LineWidth',2); end
% end
%

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

plot([vReg{3}(2,1) vReg{3}(3,1)],[vReg{3}(2,2) vReg{3}(3,2)],'k--','LineWidth',2)
plot([vReg{3}(1,1) vReg{3}(2,1)],[vReg{3}(1,2) vReg{3}(2,2)],'k--','LineWidth',2)
plot([vReg{3}(4,1) vReg{3}(1,1)],[vReg{3}(4,2) vReg{3}(1,2)],'k--','LineWidth',2)
plot([vRegB{3}(2,1) vRegB{3}(3,1)],[vRegB{3}(2,2) vRegB{3}(3,2)],'b:','LineWidth',2)
plot([vRegB{3}(1,1) vRegB{3}(2,1)],[vRegB{3}(1,2) vRegB{3}(2,2)],'b:','LineWidth',2)
plot([vRegB{3}(4,1) vRegB{3}(1,1)],[vRegB{3}(4,2) vRegB{3}(1,2)],'b:','LineWidth',2)
plot([vRegBN{3}(2,1) vRegBN{3}(3,1)],[vRegBN{3}(2,2) vRegBN{3}(3,2)],'b:','LineWidth',2)
plot([vRegBN{3}(1,1) vRegBN{3}(2,1)],[vRegBN{3}(1,2) vRegBN{3}(2,2)],'b:','LineWidth',2)
plot([vRegBN{3}(4,1) vRegBN{3}(1,1)],[vRegBN{3}(4,2) vRegBN{3}(1,2)],'b:','LineWidth',2)

plot([vReg{4}(1,1) vReg{4}(2,1)],[vReg{4}(1,2) vReg{4}(2,2)],'k--','LineWidth',2)
plot([vReg{4}(4,1) vReg{4}(1,1)],[vReg{4}(4,2) vReg{4}(1,2)],'k--','LineWidth',2)
plot([vReg{4}(3,1) vReg{4}(4,1)],[vReg{4}(3,2) vReg{4}(4,2)],'k--','LineWidth',2)
plot([vRegB{4}(1,1) vRegB{4}(2,1)],[vRegB{4}(1,2) vRegB{4}(2,2)],'b:','LineWidth',2)
plot([vRegB{4}(4,1) vRegB{4}(1,1)],[vRegB{4}(4,2) vRegB{4}(1,2)],'b:','LineWidth',2)
plot([vRegB{4}(3,1) vRegB{4}(4,1)],[vRegB{4}(3,2) vRegB{4}(4,2)],'b:','LineWidth',2)
plot([vRegBN{4}(1,1) vRegBN{4}(2,1)],[vRegBN{4}(1,2) vRegBN{4}(2,2)],'b:','LineWidth',2)
plot([vRegBN{4}(4,1) vRegBN{4}(1,1)],[vRegBN{4}(4,2) vRegBN{4}(1,2)],'b:','LineWidth',2)
plot([vRegBN{4}(3,1) vRegBN{4}(4,1)],[vRegBN{4}(3,2) vRegBN{4}(4,2)],'b:','LineWidth',2)

set(gca,'FontSize',14)
xlabel('x_r'); ylabel('y_r')
text(vReg{1}(1,1)+0.2,vReg{1}(1,2)+0.4,'r_1','FontSize',16)
text(vReg{3}(1,1)+0.2,vReg{3}(1,2)+0.4,'r_2','FontSize',16)
text(vReg{4}(1,1)+0.2,vReg{4}(1,2)+0.4,'r_3','FontSize',16)
text(vReg{1}(2,1)+0.2,vReg{1}(1,2)+0.4,'r_4','FontSize',16)

set(gca, 'Box', 'off' )
% set(gca, 'TickDir', 'out')
