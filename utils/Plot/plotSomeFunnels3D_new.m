

addpath('C:\Cornell\Research\CreatePath\ParametricVerification\ellipsoids')
addpath('C:\Users\Jon\Dropbox\Research\CreatePath\ParametricVerification\ellipsoids')


figure(101)
clf
hold on
axis equal
set(gcf,'RendererMode','manual','Renderer','OpenGL');

vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});

set(gca,'FontSize',12)
xlabel('x_r'); ylabel('y_r'); zlabel('\theta')

%%
clear Hf
Hf{1} = [];
Hlast = [];
for i = 1:length(ac_trans)
    Hlast = [Hlast; Hf{i}];
    plot(ac_trans{i},sys,101,true)
    Hf{i+1} = setdiff(get(gca,'Children'),Hlast);
end

H0 = get(gca,'Children');

%%
%clear Hf
%Hf{1} = [];
%Hlast = [];
for i = 1:length(ac_inward)
    %Hlast = [Hlast; Hf{i}];
    if ~isempty(ac_inward{i})
        for j = 1:length(ac_inward{i})
            plot(ac_inward{i}(j),sys,101,true)
        end
    end
    %Hf{i+1} = setdiff(get(gca,'Children'),Hlast);
end
H0 = get(gca,'Children');

%%
clear Hf
Hf{1} = [];
Hlast = H0;
for i = 1:length(ac_react)
    if ~isempty(ac_react{i})
        Hlast = [Hlast; Hf{i}];
        i
        for j = 1:length(ac_react{i})
            j
            plot(ac_react{i}(j),sys,101,true)
        end
        Hf{i+1} = setdiff(get(gca,'Children'),Hlast);
    end
end

H0 = get(gca,'Children');

%%
color = {'g','b','y','b','y','b','y','b','y',[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1]};
Options.shade = 0.5;
figure(101)
hold on
for j = 1:9
    Options.color = color{j};
    plot(reg(j),color{j});
end
H1 = get(gca,'Children');
H = setdiff(H1,H0);
for i = 1:length(H)
    set(H(i),'ZData',-10*ones(size(get(H(i),'XData'))));
    set(H(i),'FaceAlpha',Options.shade);
end
for j = 10:12
    plot(reg(j));
end
H2 = get(gca,'Children');
H = setdiff(H2,H1);
for i = 1:length(H)
    set(H(i),'ZData',-10*ones(size(get(H(i),'XData'))));
    set(H(i),'FaceAlpha',Options.shade);
    set(H(i),'FaceColor',[0.1,0.1,0.1]);
end
view(2)
axis([vBndMin(1) vBndMax(1) vBndMin(2) vBndMax(2)])

return

%%
% TODO: pick out points in ellFunnel to minimize obvious violations?

options.fill = 0.5;
options.shade = 1;
clear E
hold on

for itrans = 1:size(trans,1)%setdiff(1:size(transTmp,1)-2,size(transTmp,1)-4)
    options.color = max(0,min(1,[0.2 0.9 0.2] + rand(1,3)*0.4));
    plot(ac_trans{itrans},sys,101,true)
    view(2)
    pause
end
options.color = [0.9 0.2 0.2];
% plotFunnel(funnel{1,16},options)
view(2)
% pause

%%
for imode = 1:NmodesReach
    options.color = [0.2 0.2 0.9];
    for j = 1:size(funnelIn,1)
        if ~isempty(funnelIn{j,imode})
            if funnelIn{j,imode}.rho(1) > 1e-9
                plotFunnel(funnelIn{j,imode},options)
            end
        end
    end
end
view(2)
pause

for itrans = 1:size(transTmp,1)
    options.color = [0.9 0.2 0.2];
    plotFunnel(funnelJoin{1,itrans},options)
end
view(2)
pause

%%
options.color = [0.9 0.2 0.9];

% plotFunnel(funnelReactJoin{8,3},options)

for itrans = 1:size(transTmp,1)
    for j = 1:size(funnelReactJoin,1)
        if ~isempty(funnelReactJoin{j,itrans})
            options.color = [0.9 0.2 0.2];
            plotFunnel(funnelReactJoin{j,itrans},options)
        end
    end
end

%%
ii = 1;
jj = 1;
plot3([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],funnel{ii,jj}.x(end,3)*[1 1],'k','LineWidth',2)
plot3([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],funnel{ii,jj}.x(end,3)*[1 1],'k','LineWidth',2)
plot3([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],funnel{ii,jj}.x(end,3)*[1 1],'k','LineWidth',2)
plot3([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],funnel{ii,jj}.x(end,3)*[1 1],'k','LineWidth',2)

plot3([vReg{1}(3,1) vReg{1}(4,1)],[vReg{1}(3,2) vReg{1}(4,2)],funnel{ii,jj}.x(end,3)*[1 1],'k--','LineWidth',2)
% plot3([vRegB{1}(2,1) vRegB{1}(3,1)],[vRegB{1}(2,2) vRegB{1}(3,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegB{1}(3,1) vRegB{1}(4,1)],[vRegB{1}(3,2) vRegB{1}(4,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegBN{1}(2,1) vRegBN{1}(3,1)],[vRegBN{1}(2,2) vRegBN{1}(3,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegBN{1}(3,1) vRegBN{1}(4,1)],[vRegBN{1}(3,2) vRegBN{1}(4,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)

plot3([vReg{2}(1,1) vReg{2}(2,1)],[vReg{2}(1,2) vReg{2}(2,2)],funnel{ii,jj}.x(end,3)*[1 1],'k--','LineWidth',2)
plot3([vReg{2}(2,1) vReg{2}(3,1)],[vReg{2}(2,2) vReg{2}(3,2)],funnel{ii,jj}.x(end,3)*[1 1],'k--','LineWidth',2)
% plot3([vRegB{3}(1,1) vRegB{3}(2,1)],[vRegB{3}(1,2) vRegB{3}(2,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegB{3}(4,1) vRegB{3}(1,1)],[vRegB{3}(4,2) vRegB{3}(1,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegBN{3}(1,1) vRegBN{3}(2,1)],[vRegBN{3}(1,2) vRegBN{3}(2,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)
% plot3([vRegBN{3}(4,1) vRegBN{3}(1,1)],[vRegBN{3}(4,2) vRegBN{3}(1,2)],funnel{ii,jj,kkindx}.x(end,3)*[1 1],'b:','LineWidth',2)

text(vReg{1}(1,1)+0.2,vReg{1}(1,2)+0.5,funnel{ii,jj}.x(end,3),'r_1','FontSize',16)
text(vReg{1}(2,1)+0.2,vReg{1}(1,2)+0.5,funnel{ii,jj}.x(end,3),'r_2','FontSize',16)
text(vReg{3}(1,1)+0.2,vReg{3}(1,2)+0.5,funnel{ii,jj}.x(end,3),'r_3','FontSize',16)