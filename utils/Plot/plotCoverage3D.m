

addpath('C:\Cornell\Research\CreatePath\ParametricVerification\ellipsoids')
addpath('C:\Users\Jon\Dropbox\Research\CreatePath\ParametricVerification\ellipsoids')


figure(10)
clf
hold on
axis equal
set(gcf,'RendererMode','manual','Renderer','OpenGL');

vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});

set(gca,'FontSize',12)
xlabel('x_r'); ylabel('y_r'); zlabel('\theta')

%%
% TODO: pick out points in ellFunnel to minimize obvious violations?

kkindx = 1;
jj = 2;
ii = 11;
% k=1: 2, 7, 8, 11, 29, 33, 35
% k=2: 1, 3, 9

options.fill = 0.5;
options.color = [0.2 0.9 0.2];
options.shade = 1;
clear E
for k = 1:length(funnel{ii,jj,kkindx}.t)
    tmp = inv(funnel{ii,jj,kkindx}.P(:,:,k));
    tmp = (tmp+tmp')/2;
    E(k) = ellipsoid(funnel{ii,jj,kkindx}.x(k,:)',tmp*funnel{ii,jj,kkindx}.rho(k));
end
% plotEllipse(E,[]);
plot(E,options);


%%

clear H
H = hyperplane([0 0 1]',funnel{ii,jj,kkindx}.x(end,3));  % hyperplane at final theta slice

for iii = length(H)
    
    kindx = 1;
    j = 2;
    
    options.color = [0.9 0.3 0.3];
    for i = 1:size(funnelIn,1)
        clear E
        if ~isempty(funnelIn{i,j,kindx}),
            for k = 1:length(funnelIn{i,j,kindx}.t)
                tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
                tmp = (tmp+tmp')/2;
                E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
                E1(k) = hpintersection(E(k),H(iii));
                PlotEllipse1(E1(k),options.color,'waterfall');
            end
        end
    end

%     j = 3;
%     
%     options.color = [0.3 0.3 0.9];
%     for i = 1:size(funnelIn,1)
%         clear E
%         if ~isempty(funnelIn{i,j,kindx}),
%             for k = 1:length(funnelIn{i,j,kindx}.t)
%                 tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
%                 tmp = (tmp+tmp')/2;
%                 E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
%                 E1(k) = hpintersection(E(k),H(iii));
%                 PlotEllipse1(E1(k),options.color,'waterfall');
%             end
%         end
%     end
    
%     jj1 = 3;
%     
%     options.fill = 1;
%     options.color = [0.9 0.3 0.3];
%     options.shade = 0.5;
%     for i = 1:size(funnel,1)
%         clear E
%         if ~isempty(funnel{i,jj1,kindx}),
%             for k = 1:length(funnel{i,jj1,kindx}.t)
%                 tmp = inv(funnel{i,jj1,kindx}.P(:,:,k));
%                 tmp = (tmp+tmp')/2;
%                 E(k) = ellipsoid(funnel{i,jj1,kindx}.x(k,:)',tmp*funnel{i,jj1,kindx}.rho(k));
%                 E1(k) = hpintersection(E(k),H);
%                 if ~isempty(E1(k))
%                     PlotEllipse1(E1(k),options.color,'waterfall');
%                 end
%             end
%             
%         end
%     end
    
    %
%     options.fill = 0.5;
%     options.color = [0.5 0.9 0.5];
%     options.shade = 1;
%     clear E
%     for k = 1:length(funnel{ii,jj,kkindx}.t)
%         tmp = inv(funnel{ii,jj,kkindx}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnel{ii,jj,kkindx}.x(k,:)',tmp*funnel{ii,jj,kkindx}.rho(k));
%         E1(k) = hpintersection(E(k),H(iii));
%         PlotEllipse1(E1(k),options.color,'waterfall');
%     end
%     plot3(trajFunnel{ii,jj,kkindx}(:,1),trajFunnel{ii,jj,kkindx}(:,2),funnel{ii,jj,kkindx}.x(end,3)*ones(size(trajFunnel{ii,jj,kkindx}(:,1))),'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',1);
    
end

%%

plot3([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'k--','LineWidth',2)
plot3([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'k--','LineWidth',2)
plot3([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'k--','LineWidth',2)
plot3([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'k--','LineWidth',2)

plot3([vReg{1}(2,1) vReg{1}(3,1)],[vReg{1}(2,2) vReg{1}(3,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'k--','LineWidth',2)
plot3([vReg{1}(3,1) vReg{1}(4,1)],[vReg{1}(3,2) vReg{1}(4,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'k--','LineWidth',2)
% plot3([vRegB{1}(2,1) vRegB{1}(3,1)],[vRegB{1}(2,2) vRegB{1}(3,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegB{1}(3,1) vRegB{1}(4,1)],[vRegB{1}(3,2) vRegB{1}(4,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegBN{1}(2,1) vRegBN{1}(3,1)],[vRegBN{1}(2,2) vRegBN{1}(3,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegBN{1}(3,1) vRegBN{1}(4,1)],[vRegBN{1}(3,2) vRegBN{1}(4,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)

plot3([vReg{3}(1,1) vReg{3}(2,1)],[vReg{3}(1,2) vReg{3}(2,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'k--','LineWidth',2)
plot3([vReg{3}(4,1) vReg{3}(1,1)],[vReg{3}(4,2) vReg{3}(1,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'k--','LineWidth',2)
% plot3([vRegB{3}(1,1) vRegB{3}(2,1)],[vRegB{3}(1,2) vRegB{3}(2,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegB{3}(4,1) vRegB{3}(1,1)],[vRegB{3}(4,2) vRegB{3}(1,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegBN{3}(1,1) vRegBN{3}(2,1)],[vRegBN{3}(1,2) vRegBN{3}(2,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegBN{3}(4,1) vRegBN{3}(1,1)],[vRegBN{3}(4,2) vRegBN{3}(1,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)

text(vReg{1}(1,1)+0.2,vReg{1}(1,2)+0.5,funnel{ii,jj,kkindx}.x(end,3),'r_1','FontSize',16)
text(vReg{1}(2,1)+0.2,vReg{1}(1,2)+0.5,funnel{ii,jj,kkindx}.x(end,3),'r_2','FontSize',16)
text(vReg{3}(1,1)+0.2,vReg{3}(1,2)+0.5,funnel{ii,jj,kkindx}.x(end,3),'r_3','FontSize',16)

%%

figure(11)
clf
hold on
axis equal
set(gcf,'RendererMode','manual','Renderer','OpenGL');

xlabel('x_r'),  ylabel('y_r')

kindx = 1;


j = 2;

options.fill = 1;
options.color = [0.9 0.3 0.3];
options.shade = 0.5;
for i = 1:size(ellFunnelIn,1)
    clear E
    if ~isempty(funnelIn{i,j,kindx}),
        for k = 1:length(funnelIn{i,j,kindx}.t)
            tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
            tmp = (tmp+tmp')/2;
            E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
            E1(k) = hpintersection(E(k),H);
            if ~isempty(E1(k))
                PlotEllipse1(E1(k),options.color,'2d');
                PlotEllipse1(E1(k),options.color,'2d_outline');
            end
        end
    end
end

% j = 3;
% 
% options.fill = 1;
% options.color = [0.3 0.3 0.9];
% options.shade = 0.5;
% for i = 1:size(ellFunnelIn,1)
%     clear E
%     if ~isempty(funnelIn{i,j,kindx}),
%         for k = 1:length(funnelIn{i,j,kindx}.t)
%             tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
%             tmp = (tmp+tmp')/2;
%             E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
%             E1(k) = hpintersection(E(k),H);
%             if ~isempty(E1(k))
%                 PlotEllipse1(E1(k),options.color,'2d');
%                 PlotEllipse1(E1(k),options.color,'2d_outline');
%             end
%         end
%     end
% end

% jj1 = 3;
% 
% options.fill = 1;
% options.color = [0.9 0.3 0.3];
% options.shade = 0.5;
% for i = 1:size(funnel,1)    
%     clear E
%     if ~isempty(funnel{i,jj1,kindx}),
%         for k = 1:length(funnel{i,jj1,kindx}.t)
%             tmp = inv(funnel{i,jj1,kindx}.P(:,:,k));
%             tmp = (tmp+tmp')/2;
%             E(k) = ellipsoid(funnel{i,jj1,kindx}.x(k,:)',tmp*funnel{i,jj1,kindx}.rho(k));
%             E1(k) = hpintersection(E(k),H);
%             if ~isempty(E1(k))
%                 PlotEllipse1(E1(k),options.color,'2d');
%                 PlotEllipse1(E1(k),options.color,'2d_outline');
%             end
%         end
%     end
% end

clear E
options.fill = 1;
options.color = [0.5 0.9 0.5];
options.shade = 0.5;
E2 = [];
for k = 1:length(funnel{ii,jj,kindx}.t)
    tmp = inv(funnel{ii,jj,kindx}.P(:,:,k));
    tmp = (tmp+tmp')/2;
    E(k) = ellipsoid(funnel{ii,jj,kindx}.x(k,:)',tmp*funnel{ii,jj,kindx}.rho(k));
    E1(k) = hpintersection(E(k),H);
    if ~isempty(E1(k))
        PlotEllipse1(E1(k),options.color,'2d',0.25);
    end
end

%%

plot([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
plot([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
plot([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
plot([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])

plot([vReg{1}(2,1) vReg{1}(3,1)],[vReg{1}(2,2) vReg{1}(3,2)],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
plot([vReg{1}(3,1) vReg{1}(4,1)],[vReg{1}(3,2) vReg{1}(4,2)],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
% plot([vRegB{1}(2,1) vRegB{1}(3,1)],[vRegB{1}(2,2) vRegB{1}(3,2)],'b:','LineWidth',2)
% plot([vRegB{1}(3,1) vRegB{1}(4,1)],[vRegB{1}(3,2) vRegB{1}(4,2)],'b:','LineWidth',2)
% plot([vRegBN{1}(2,1) vRegBN{1}(3,1)],[vRegBN{1}(2,2) vRegBN{1}(3,2)],'b:','LineWidth',2)
% plot([vRegBN{1}(3,1) vRegBN{1}(4,1)],[vRegBN{1}(3,2) vRegBN{1}(4,2)],'b:','LineWidth',2)

plot([vReg{3}(1,1) vReg{3}(2,1)],[vReg{3}(1,2) vReg{3}(2,2)],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
plot([vReg{3}(4,1) vReg{3}(1,1)],[vReg{3}(4,2) vReg{3}(1,2)],'k--','LineWidth',2,'Color',[0.5 0.5 0.5])
% plot([vRegB{3}(1,1) vRegB{3}(2,1)],[vRegB{3}(1,2) vRegB{3}(2,2)],'b:','LineWidth',2)
% plot([vRegB{3}(4,1) vRegB{3}(1,1)],[vRegB{3}(4,2) vRegB{3}(1,2)],'b:','LineWidth',2)
% plot([vRegBN{3}(1,1) vRegBN{3}(2,1)],[vRegBN{3}(1,2) vRegBN{3}(2,2)],'b:','LineWidth',2)
% plot([vRegBN{3}(4,1) vRegBN{3}(1,1)],[vRegBN{3}(4,2) vRegBN{3}(1,2)],'b:','LineWidth',2)

text(vReg{1}(1,1)+0.2,vReg{1}(1,2)+0.5,'r_1','FontSize',16)
text(vReg{1}(2,1)+0.2,vReg{1}(1,2)+0.5,'r_2','FontSize',16)
text(vReg{3}(1,1)+0.2,vReg{3}(1,2)+0.5,'r_3','FontSize',16)


%%
% 
% xslice = 3.5:0.1:4.5;
% clear H
% i = 0;
% for x = xslice
%     i = i + 1;
%     H(i) = hyperplane([0 1 0]',x);
% end
% 
% 
% for indx = 1:ceil(length(H)/6)
%     figure(20+indx)
%     clf
% end
% 
% for iii = 1:length(H)
%     
%     figure(20+ceil(iii/6))
%     subplot(3,2,mod(iii-1,6)+1)
%     hold on
%     xlabel('x_r'),  ylabel('\theta')
%     title(['Slice at y_r = ',num2str(xslice(iii))])
% 
%     j = 2;
%     kindx = 1;
%     
%     options.fill = 1;
%     options.color = [0.9 0.3 0.3];
%     options.shade = 0.5;
%     for i = 1:size(ellFunnelC,1)
%         clear E
%         if ~isempty(funnelIn{i,j,kindx}),
%             for k = 1:length(funnelIn{i,j,kindx}.t)
%                 tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
%                 tmp = (tmp+tmp')/2;
%                 E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
%                 E1(k) = hpintersection(E(k),H(iii));
%                 PlotEllipse1(E1(k),options.color,'2d');
%             end
%         end
%     end
%     
%     j = 3;
%     kindx = 1;
%     
%     options.fill = 1;
%     options.color = [0.3 0.3 0.9];
%     options.shade = 0.5;
%     for i = 1:size(ellFunnelC,1)
%         clear E
%         if ~isempty(funnelIn{i,j,kindx}),
%             for k = 1:length(funnelIn{i,j,kindx}.t)
%                 tmp = inv(funnelIn{i,j,kindx}.P(:,:,k));
%                 tmp = (tmp+tmp')/2;
%                 E(k) = ellipsoid(funnelIn{i,j,kindx}.x(k,:)',tmp*funnelIn{i,j,kindx}.rho(k));
%                 E1(k) = hpintersection(E(k),H(iii));
%                 PlotEllipse1(E1(k),options.color,'2d');
%             end
%         end
%     end
%     
%     kindx = 2;
%     
%     clear E
%     options.fill = 1;
%     options.color = [0.5 0.9 0.5];
%     options.shade = 0.5;
%     for k = 1:length(funnel{ii,jj,kindx}.t)
%         tmp = inv(funnel{ii,jj,kindx}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnel{ii,jj,kindx}.x(k,:)',tmp*funnel{ii,jj,kindx}.rho(k));
%         E1(k) = hpintersection(E(k),H(iii));
%         PlotEllipse1(E1(k),options.color,'2d');
%     end
% end
